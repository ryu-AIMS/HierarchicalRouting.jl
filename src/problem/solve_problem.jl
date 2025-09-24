
"""
    initial_solution(
        problem::Problem;
        k::Int=1,
        disturbance_clusters::Set{Int64}=Set{Int64}(),
        waypoint_optim_method=nothing,
        seed::Union{Nothing,Int64}=nothing,
        rng::AbstractRNG=Random.GLOBAL_RNG,
    )::MSTSolution

Generate an initial solution to the problem for:
- the mothership by using the nearest neighbour heuristic
and then improving using the 2-opt heuristic, and
-  for tenders using the sequential nearest neighbour heuristic.
Optionally, optimize waypoints using a provided optimization method.

# Arguments
- `problem`: Problem instance to solve
- `k`: Number of clusters to generate
- `disturbance_clusters`: Set of sequenced clusters to simulate disturbances before.
- `waypoint_optim_method`: Function to use in waypoint optimization.
- `seed`: Optional seed for random number generation
- `rng`: Random number generator

# Returns
Best total MSTSolution found
"""
function initial_solution(
    problem::Problem;
    k::Int=1,
    disturbance_clusters::Set{Int64}=Set{Int64}(),
    waypoint_optim_method=nothing,
    seed::Union{Nothing,Int64}=nothing,
    rng::AbstractRNG=Random.GLOBAL_RNG,
)::MSTSolution
    return solve(
        problem;
        k,
        disturbance_clusters,
        seed,
        rng,
        waypoint_optim_method,
        do_improve=false
    )
end

"""
    solve(
        problem::Problem;
        k::Int=1,
        disturbance_clusters::Set{Int64}=Set{Int64}()
        seed::Union{Nothing,Int64}=nothing,
        rng::AbstractRNG=Random.GLOBAL_RNG,
        waypoint_optim_method=nothing,
    )::MSTSolution

Generate a solution to the problem for:
- the mothership by using the nearest neighbour heuristic to generate an initial solution,
    and then improving using the 2-opt heuristic, and
-  for tenders using the sequential nearest neighbour heuristic.
Optimize solution using simulated annealing as default optimization method.
Optionally, optimize waypoints using a set or provided optimization method.

# Arguments
- `problem`: Problem instance to solve
- `k`: Number of clusters to generate
- `disturbance_clusters`: Set of sequenced clusters to simulate disturbances before.
- `seed`: Optional seed for random number generation
- `rng`: AbstractRNG for random number generation
- `waypoint_optim_method`: Function to use in waypoint optimization.

# Returns
Best total MSTSolution found
"""
function solve(
    problem::Problem;
    k::Int=1,
    disturbance_clusters::Set{Int64}=Set{Int64}(),
    seed::Union{Nothing,Int64}=nothing,
    rng::AbstractRNG=Random.GLOBAL_RNG,
    waypoint_optim_method=nothing,
    do_improve::Bool=true,
)::MSTSolution
    if !isnothing(seed)
        Random.seed!(rng, seed)
    end

    ordered_disturbances = sort(unique(disturbance_clusters))
    n_events = length(ordered_disturbances) + 1
    cluster_sets = Vector{Vector{Cluster}}(undef, n_events)
    ms_soln_sets = Vector{MothershipSolution}(undef, n_events)
    tender_soln_sets = Vector{Vector{TenderSolution}}(undef, n_events)
    total_tender_capacity = Int(problem.tenders.number * problem.tenders.capacity)
    vessel_weightings = (problem.mothership.weighting, problem.tenders.weighting)

    # Cluster the problem data
    clusters::Vector{Cluster} = cluster_problem(problem; k)
    cluster_centroids_df::DataFrame = generate_cluster_df(clusters, problem.depot)

    # Route the mothership using nearest neighbour and 2-opt
    ms_route::MothershipSolution = optimize_mothership_route(problem, cluster_centroids_df)
    clust_seq::Vector{Int64} = filter(
        c -> c != 0 && c <= length(clusters),
        ms_route.cluster_sequence.id
    )

    # Generate initial tender solutions using sequential nearest neighbour
    initial_tenders = [
        tender_sequential_nearest_neighbour(
            clusters[clust_seq][j],
            (ms_route.route.nodes[2j], ms_route.route.nodes[2j+1]),
            problem.tenders.number,
            problem.tenders.capacity,
            problem.tenders.exclusion.geometry
        )
        for j in 1:length(clust_seq)
    ]

    if do_improve
        @info "Improving initial solution using simulated annealing"
        # Optimize the initial tenders solution up to the first disturbance
        next_cluster_idx = !isempty(ordered_disturbances) ?
                           ordered_disturbances[1] :
                           length(clusters)

        optimized_initial, _ = improve_solution(
            MSTSolution([clusters], [ms_route], [initial_tenders]),
            problem.mothership.exclusion.geometry,
            problem.tenders.exclusion.geometry,
            1, next_cluster_idx, vessel_weightings
        )

        # Apply the optimized initial solution to the first set of clusters pre-disturbance
        clusters = optimized_initial.cluster_sets[1]
        ms_route = optimized_initial.mothership_routes[1]
        initial_tenders = optimized_initial.tenders[1]
        clust_seq = filter(
            c -> c != 0 && c <= length(clusters),
            ms_route.cluster_sequence.id
        )
    end

    # Apply solution to the first set of clusters pre-disturbance
    disturb_idx = 1
    cluster_sets[disturb_idx] = clusters
    ms_soln_sets[disturb_idx] = ms_route
    tender_soln_sets[disturb_idx] = initial_tenders

    # Simulate disturbance events
    solution::MSTSolution = _apply_disturbance_events!(
        cluster_sets,
        ms_soln_sets,
        tender_soln_sets,
        clust_seq,
        ordered_disturbances,
        problem,
        total_tender_capacity;
        do_improve,
        waypoint_optim_method
    )

    return solution
end


"""
    _apply_disturbance_events!(
        cluster_sets::Vector{Vector{Cluster}},
        ms_soln_sets::Vector{MothershipSolution},
        tender_soln_sets::Vector{Vector{TenderSolution}},
        clust_seq::Vector{Int64},
        ordered_disturbances::Vector{Int64},
        problem::Problem,
        total_tender_capacity::Int;
        do_improve::Bool=false,
        waypoint_optim_method=nothing,
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
    total_tender_capacity::Int;
    do_improve::Bool=false,
    waypoint_optim_method=nothing,
)::MSTSolution
    disturb_idx = 1
    clusters::Vector{Cluster} = cluster_sets[disturb_idx]
    ms_route::MothershipSolution = ms_soln_sets[disturb_idx]
    solution::MSTSolution = MSTSolution(cluster_sets, ms_soln_sets, tender_soln_sets)

    # Optimize initial waypoints
    if !isnothing(waypoint_optim_method)
        @info "Optimizing pre-disturbed waypoint subset using $(waypoint_optim_method)"

        solution = optimize_waypoints(
            MSTSolution([clusters], [ms_route], [tender_soln_sets[disturb_idx]]),
            problem,
            waypoint_optim_method
        )
        clusters = solution.cluster_sets[disturb_idx]
        ms_route = solution.mothership_routes[disturb_idx]
        tender_soln_sets[disturb_idx] = solution.tenders[disturb_idx]
    end

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
                candidate_wpt_idxs
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

"""
    optimize_mothership_route(
        problem::Problem,
        cluster_centroids_df::DataFrame
    )::MothershipSolution
    optimize_mothership_route(
        problem::Problem,
        cluster_centroids_df::DataFrame,
        cluster_seq_idx::Int64,
        ms_route::MothershipSolution,
        cluster_ids_visited::Vector{Int64}
    )::MothershipSolution

Generate an optimized mothership route using the nearest neighbour heuristic and 2-opt for:
- the whole mothership route to/from the depot, or
- the remaining/partial mothership route to clusters based on current position.

# Arguments
- `problem`: Problem instance to solve
- `cluster_centroids_df`: DataFrame containing cluster centroids
- `cluster_seq_idx`: Index of the cluster sequence
- `ms_route`: Current mothership route
- `cluster_ids_visited`: Vector of cluster IDs that have been visited

# Returns
- The optimized mothership route as a `MothershipSolution` object.
"""
function optimize_mothership_route(
    problem::Problem,
    cluster_centroids_df::DataFrame
)::MothershipSolution
    # Nearest Neighbour to generate initial mothership route & matrix
    ms_soln_NN::MothershipSolution = nearest_neighbour(
        cluster_centroids_df,
        problem.mothership.exclusion.geometry,
        problem.tenders.exclusion.geometry
    )

    # 2-opt to improve the NN soln
    ms_soln_2opt::MothershipSolution = two_opt(
        ms_soln_NN,
        problem.mothership.exclusion.geometry,
        problem.tenders.exclusion.geometry
    )
    return ms_soln_2opt
end
function optimize_mothership_route(
    problem::Problem,
    cluster_centroids_df::DataFrame,
    cluster_seq_idx::Int64,
    ms_route::MothershipSolution,
    cluster_ids_visited::Vector{Int64}
)::MothershipSolution
    start_point::Point{2,Float64} = ms_route.route.nodes[2*cluster_seq_idx-1]

    remaining_clusters_df::DataFrame = filter(
        row -> row.id ∉ cluster_ids_visited,
        cluster_centroids_df
    )

    # Nearest Neighbour to generate initial mothership route & matrix
    ms_soln_NN::MothershipSolution = nearest_neighbour(
        remaining_clusters_df,
        problem.mothership.exclusion.geometry,
        problem.tenders.exclusion.geometry,
        start_point,
        ms_route,
        cluster_seq_idx
    )

    # 2-opt to improve the NN soln
    ms_soln_2opt::MothershipSolution = two_opt(
        ms_soln_NN,
        problem.mothership.exclusion.geometry,
        problem.tenders.exclusion.geometry,
        cluster_seq_idx
    )

    return ms_soln_2opt
end

"""
    improve_solution(
        initial_solution::MSTSolution,
        exclusions_mothership::POLY_VEC,
        exclusions_tender::POLY_VEC,
        current_cluster_idx::Int,
        next_cluster_idx::Int,
        vessel_weightings::NTuple{2,AbstractFloat};
        opt_function::Function=simulated_annealing,
        objective_function::Function=critical_path,
        perturb_function::Function=perturb_swap_solution,
        max_iterations::Int=1_000,
        temp_init::Float64=500.0,
        cooling_rate::Float64=0.95,
        static_limit::Int=20,
    )::Tuple{MSTSolution,Float64}
    improve_solution(
        init_solution::MSTSolution,
        problem::Problem;
        opt_function::Function=simulated_annealing,
        objective_function::Function=critical_path,
        perturb_function::Function=perturb_swap_solution,
        max_iterations::Int=1_000,
        temp_init::Float64=500.0,
        cooling_rate::Float64=0.95,
        static_limit::Int=20,
    )

Improve the solution using the optimization function `opt_function` with the objective
function `objective_function` and the perturbation function `perturb_function`.\n
Multiple dispatch to improve full and partial solutions (respectively).

# Arguments
- `initial_solution`: Initial solution to improve
- `exclusions_mothership`: Exclusion zone polygon polygons for the mothership
- `exclusions_tender`: Exclusion zone polygon polygons for the tenders
- `current_cluster_idx`: Index of the current cluster in the sequence
- `next_cluster_idx`: Index of the next cluster in the sequence
- `opt_function`: Optimization function to improve the solution
- `objective_function`: Objective function to quantify and evaluate the solution
- `perturb_function`: Perturbation function to generate changes in the solution
- `max_iterations`: Maximum number of iterations
- `temp_init`: Initial temperature for simulated annealing
- `cooling_rate`: Cooling rate for simulated annealing
- `static_limit`: Number of iterations to allow stagnation before early exit
- `vessel_weightings`: Weightings for the mothership and tenders in the objective function

# Returns
- `soln_best`: Solution with lowest objective value found
- `z_best`: Objective value associated with the best solution
"""
function improve_solution(
    initial_solution::MSTSolution,
    exclusions_mothership::POLY_VEC,
    exclusions_tender::POLY_VEC,
    current_cluster_idx::Int,
    next_cluster_idx::Int,
    vessel_weightings::NTuple{2,AbstractFloat};
    opt_function::Function=simulated_annealing,
    objective_function::Function=critical_path,
    perturb_function::Function=perturb_swap_solution,
    max_iterations::Int=1_000,
    temp_init::Float64=500.0,
    cooling_rate::Float64=0.95,
    static_limit::Int=20,
)::Tuple{MSTSolution,Float64}
    current_mothership_route::MothershipSolution = initial_solution.mothership_routes[end]

    cluster_seq_ids::Vector{Int64} = current_mothership_route.cluster_sequence.id
    final_cluster_idx = next_cluster_idx == length(cluster_seq_ids) - 1 ?
                        next_cluster_idx - 1 :
                        next_cluster_idx

    clust_seq_current::Vector{Int64} = cluster_seq_ids[
        current_cluster_idx+1:final_cluster_idx+1
    ]
    clust_seq_noncurrent::Vector{Int64} = setdiff(cluster_seq_ids, clust_seq_current)
    filter!(!=(0), clust_seq_noncurrent)

    cluster_set::Vector{Cluster} = initial_solution.cluster_sets[end]
    current_clusters::Vector{Cluster} = cluster_set[clust_seq_current]
    sort!(current_clusters, by=c -> c.id)

    tender_set::Vector{TenderSolution} = initial_solution.tenders[end]
    current_tender_routes::Vector{TenderSolution} = tender_set[current_cluster_idx:final_cluster_idx]
    noncurrent_tender_routes::Vector{TenderSolution} = setdiff(
        tender_set,
        current_tender_routes
    )

    current_solution = MSTSolution(
        [current_clusters],
        [current_mothership_route],
        [current_tender_routes]
    )

    soln_best_partial, z_best = opt_function(
        current_solution,
        objective_function,
        perturb_function,
        exclusions_mothership,
        exclusions_tender,
        max_iterations,
        temp_init,
        cooling_rate,
        static_limit;
        vessel_weightings
    )

    merged_clusters = vcat(
        cluster_set[clust_seq_noncurrent],
        soln_best_partial.cluster_sets[1]
    )
    sort!(merged_clusters, by=c -> c.id)

    merged_tenders = vcat(soln_best_partial.tenders[1], noncurrent_tender_routes)
    sort!(merged_tenders, by=t -> t.id)
    interior_ids = @view cluster_seq_ids[2:end-1]
    ordered_tenders = merged_tenders[interior_ids]

    soln_best = MSTSolution(
        [merged_clusters],
        [soln_best_partial.mothership_routes[1]],
        [ordered_tenders]
    )

    return soln_best, z_best
end
function improve_solution(
    init_solution::MSTSolution,
    problem::Problem;
    opt_function::Function=simulated_annealing,
    objective_function::Function=critical_path,
    perturb_function::Function=perturb_swap_solution,
    max_iterations::Int=1_000,
    temp_init::Float64=500.0,
    cooling_rate::Float64=0.95,
    static_limit::Int=20,
)
    vessel_weightings::NTuple{2,AbstractFloat} = (
        problem.mothership.weighting,
        problem.tenders.weighting
    )
    current_cluster_idx::Int64 = 1
    next_cluster_idx::Int64 = length(init_solution.cluster_sets[end])

    return improve_solution(
        init_solution,
        problem.mothership.exclusion.geometry,
        problem.tenders.exclusion.geometry,
        current_cluster_idx,
        next_cluster_idx,
        vessel_weightings;
        opt_function,
        objective_function,
        perturb_function,
        max_iterations,
        temp_init,
        cooling_rate,
        static_limit
    )
end
