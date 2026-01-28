
"""
    disturb_clusters(
        remaining_clusters::Vector{Cluster},
        disturbance_df::DataFrame,
        current_location::Point{2,Float64},
        exclusions::POLY_VEC,
        total_tender_capacity::Int;
        tol::Float64=0.01
    )::Vector{Cluster}

Disturb the clusters by randomly removing nodes from them.

# Arguments
- `remaining_clusters`: A vector of `Cluster` objects representing clusters to be disturbed.
- `disturbance_df`: A DataFrame containing disturbance data for each node.
- `current_location`: The current location of the mothership at the time of disturbance.
- `exclusions`: A DataFrame containing exclusion zones that should not be disturbed.
- `total_tender_capacity`: The total capacity available for the tender fleet.
- `tol`: A tolerance value for distance calculations.

# Returns
A vector of `Cluster` objects with nodes removed from them.
"""
function disturb_clusters(
    remaining_clusters::Vector{Cluster},
    disturbance_df::DataFrame,
    current_location::Point{2,Float64},
    exclusions::POLY_VEC,
    total_tender_capacity::Int;
    tol::Float64=0.01
)::Vector{Cluster}
    num_clusters = length(remaining_clusters)
    remaining_nodes = [node for cluster in remaining_clusters for node in cluster.nodes]

    # Filter the disturbance_df so that it only includes rows with nodes in remaining_nodes
    filtered_df = filter(row -> row.node in remaining_nodes, disturbance_df)

    if nrow(filtered_df) != length(remaining_nodes)
        @warn "Warning: Expected $(length(remaining_nodes)) nodes, but filtered df contains $(nrow(filtered_df))"
    end

    disturbed_targets = disturb_remaining_clusters(
        filtered_df,
        num_clusters,
        current_location,
        exclusions,
        total_tender_capacity;
        tol
    )

    # Update the cluster assignments based on previous numbering
    updated_disturbed_targets = update_cluster_assignments(
        disturbed_targets,
        Dict(c.id => c.centroid for c in remaining_clusters)
    )

    disturbed_clusters = calculate_cluster_centroids(
        updated_disturbed_targets;
    )

    return disturbed_clusters
end

"""
    _apply_disturbance_events!(
        cluster_sets::Vector{Vector{Cluster}},
        ms_soln_sets::Vector{MothershipSolution},
        tender_soln_sets::Vector{Vector{TenderSolution}},
        clust_seq::Vector{Int64},
        ordered_disturbances::Vector{Int64},
        problem::Problem,
        total_tender_capacity::Int,
        time_limit::Float64;
        do_improve::Bool=false,
        waypoint_optim_method=nothing,
        wpt_optim_plot_flag::Bool=false,
    )::MSTSolution

Simulates and applies disturbance events to the solution.
Shared disturbance-handling loop used by `initial_solution` and `solve`.
Updates the `cluster_sets`, `ms_soln_sets`, and `tender_soln_sets` at each disturbance event.
- If `do_improve=true`, additionally runs `improve_solution(...)` at each disturbance, using
  vessel_weightings derived from the `problem` instance.
- Assumes index 1 (pre-disturbance state) has already been written into the *_sets vectors.
- `waypoint_optim_method` can be provided to optimize waypoints between disturbance events.
    - NB: Currently, partial optimization perturbs ALL future/unvisited waypoints, rather than
    just those appearing in `candidate_wpt_idxs` and/or between disturbance events.
"""
function _apply_disturbance_events!(
    cluster_sets::Vector{Vector{Cluster}},
    ms_soln_sets::Vector{MothershipSolution},
    tender_soln_sets::Vector{Vector{TenderSolution}},
    clust_seq::Vector{Int64},
    ordered_disturbances::Vector{Int64},
    problem::Problem,
    total_tender_capacity::Int,
    time_limit::Float64;
    do_improve::Bool=false,
    waypoint_optim_method=nothing,
    wpt_optim_plot_flag::Bool=false,
)::MSTSolution
    disturb_idx = 1
    clusters::Vector{Cluster} = cluster_sets[disturb_idx]
    ms_route::MothershipSolution = ms_soln_sets[disturb_idx]
    solution::MSTSolution = MSTSolution(cluster_sets, ms_soln_sets, tender_soln_sets)

    # Iterate through each disturbance event and update solution
    for disturb_clust_idx ∈ ordered_disturbances
        cluster_id = clust_seq[disturb_clust_idx]
        cluster_letter = generate_letter_id(cluster_id)
        @info """Disturbance event #$disturb_idx at
        \t$(ms_route.route.nodes[2*disturb_clust_idx-1])
        \tbefore $(disturb_clust_idx)th cluster_id=$(cluster_letter)=$(cluster_id)"""

        # Update clusters based on the impact of disturbance event on future points/clusters
        clusters = vcat(
            clusters[clust_seq][1:disturb_clust_idx-1],
            disturb_clusters(
                clusters[clust_seq][disturb_clust_idx:end],
                problem.targets.disturbance_gdf,
                ms_route.route.nodes[2*disturb_clust_idx-1],
                problem.tenders.exclusion.geometry,
                total_tender_capacity
            )
        )
        sort!(clusters, by=x -> x.id)

        removed_nodes = setdiff(
            vcat([c.nodes for c in cluster_sets[disturb_idx]]...),
            vcat([c.nodes for c in clusters]...)
        )
        if !isempty(removed_nodes)
            @info """Removed nodes due to disturbance event (since previous cluster):
            \t$(join(removed_nodes, "\n\t"))"""
        end

        # Re-generate the cluster centroids to route mothership
        cluster_centroids_df = generate_cluster_df(clusters, problem.depot)

        # Re-route mothership (respecting pre-existing portion as fixed)
        ms_route = optimize_mothership_route(
            problem,
            cluster_centroids_df,
            disturb_clust_idx,
            ms_route,
            getfield.(clusters[clust_seq][1:disturb_clust_idx-1], :id)
        )
        clust_seq = filter(
            c -> c != 0 && c <= length(clusters),
            ms_route.cluster_sequence.id
        )

        # Update tender solutions (reuse before the disturbance, recompute at/after)
        current_tender_soln = Vector{TenderSolution}(undef, length(clust_seq))

        # Generate tender solutions for the current disturbance cluster set
        for j in 1:length(clust_seq)
            if j < disturb_clust_idx
                current_tender_soln[j] = tender_soln_sets[disturb_idx][j]
            else
                current_tender_soln[j] = tender_sequential_nearest_neighbour(
                    clusters[clust_seq][j],
                    (ms_route.route.nodes[2j], ms_route.route.nodes[2j+1]),
                    problem.tenders.number,
                    problem.tenders.capacity,
                    problem.tenders.exclusion.geometry
                )
            end
        end

        # Increment event index
        disturb_idx += 1

        # Optimize mothership waypoints between disturbance events
        if !isnothing(waypoint_optim_method)
            # Optimize ALL future wpts #! NOTE: not just to next disturbance event
            clust_selection = (ordered_disturbances[disturb_idx-1]:length(clusters)+1)
            candidate_wpt_idxs = 2*clust_selection[1]:2*clust_selection[end]-1

            @info """Optimizing waypoints $candidate_wpt_idxs for clusters $clust_selection
                    using $waypoint_optim_method"""

            solution_tmp::MSTSolution = optimize_waypoints(
                MSTSolution([clusters], [ms_route], [current_tender_soln]),
                problem,
                waypoint_optim_method,
                candidate_wpt_idxs;
                time_limit,
                plot_flag=wpt_optim_plot_flag,
            )

            clusters = solution_tmp.cluster_sets[1]
            ms_route = solution_tmp.mothership_routes[1]
            current_tender_soln = solution_tmp.tenders[1]
        end

        # Solution improvement step (used by `solve`, not by `initial_solution`)
        if do_improve
            next_disturbance_cluster_idx =
                disturb_clust_idx <= length(ordered_disturbances) ?
                ordered_disturbances[disturb_clust_idx] :
                length(clusters) + 1

            vessel_weightings = (problem.mothership.weighting, problem.tenders.weighting)

            optimized_current_solution, _ = improve_solution(
                MSTSolution([clusters], [ms_route], [current_tender_soln]),
                problem.mothership.exclusion.geometry,
                problem.tenders.exclusion.geometry,
                disturb_clust_idx,
                next_disturbance_cluster_idx,
                vessel_weightings,
            )
            # Overwrite with improved
            clusters = optimized_current_solution.cluster_sets[1]
            ms_route = optimized_current_solution.mothership_routes[1]
            current_tender_soln = optimized_current_solution.tenders[1]
            # Update clust_seq in case that it has changed post-improvement
            clust_seq = filter(
                c -> c != 0 && c <= length(clusters),
                ms_route.cluster_sequence.id
            )
        end

        # Update solution sets
        cluster_sets[disturb_idx] = clusters
        ms_soln_sets[disturb_idx] = ms_route
        tender_soln_sets[disturb_idx] = current_tender_soln
    end

    if !isempty(ordered_disturbances)
        solution = MSTSolution(cluster_sets, ms_soln_sets, tender_soln_sets)
    end

    return solution
end
