
"""
    Recompute the routes for sorties in a cluster after a perturbation.
"""
function _recompute_sortie_routes(
    original_sorties::Vector{Route},
    modified_sortie_idx::Int,
    modified_nodes::Vector{Point{2,Float64}},
    start_point::Point{2,Float64},
    finish_point::Point{2,Float64},
    exclusions::POLY_VEC
)::Vector{Route}
    new_sorties = Vector{Route}(undef, length(original_sorties))

    for (s_idx, old_sortie) in enumerate(original_sorties)
        nodes_new = s_idx == modified_sortie_idx ? modified_nodes : copy(old_sortie.nodes)
        tender_tour = [[start_point]; nodes_new; [finish_point]]
        dist_mat = get_feasible_matrix(tender_tour, exclusions)[1]
        path_vec = get_feasible_vector(tender_tour, exclusions)[2]
        new_sorties[s_idx] = Route(
            nodes_new,
            dist_mat,
            vcat(path_vec...)
        )
    end
    return new_sorties
end

"""
    _build_sortie_route(
        sortie_nodes::Vector{Point{2,Float64}},
        start_pt::Point{2,Float64},
        finish_pt::Point{2,Float64},
        exclusions::POLY_VEC
    )::Route

Compute a feasible `Route` for a sortie given its interior nodes and start/finish waypoints.
Used when sortie node membership changes (cross-cluster swap), for which `rebuild_sortie`
is insufficient (it only updates the terminal legs, assuming interior nodes are unchanged).
"""
function _build_sortie_route(
    sortie_nodes::Vector{Point{2,Float64}},
    start_pt::Point{2,Float64},
    finish_pt::Point{2,Float64},
    exclusions::POLY_VEC
)::Route
    full_tour::Vector{Point{2,Float64}} = [[start_pt]; sortie_nodes; [finish_pt]]
    dist_vec, path_vec = get_feasible_vector(full_tour, exclusions)
    return Route(sortie_nodes, dist_vec, vcat(path_vec...))
end

""" Resolve new waypoints for tenders affected by a cross-cluster perturbation. """
function _resolve_tender_endpoints(
    updated_waypoints::DataFrame,
    cluster_a_idx::Int,
    cluster_b_idx::Int
)::NTuple{4,Point}
    # Identify new start/finish waypoints for the two modified tenders
    cc::Vector{NTuple{2,Int64}} = updated_waypoints.connecting_clusters
    cc_first::Vector{Int64}, cc_last::Vector{Int64} = first.(cc), last.(cc)

    tender_a_new_start_idx::Int = findlast(cc_last .== cluster_a_idx)
    tender_a_new_finish_idx::Int = findfirst(cc_first .== cluster_a_idx)
    tender_b_new_start_idx::Int = findlast(cc_last .== cluster_b_idx)
    tender_b_new_finish_idx::Int = findfirst(cc_first .== cluster_b_idx)

    tender_a_new_start = updated_waypoints.waypoint[tender_a_new_start_idx]
    tender_a_new_finish = updated_waypoints.waypoint[tender_a_new_finish_idx]
    tender_b_new_start = updated_waypoints.waypoint[tender_b_new_start_idx]
    tender_b_new_finish = updated_waypoints.waypoint[tender_b_new_finish_idx]

    return tender_a_new_start, tender_a_new_finish, tender_b_new_start, tender_b_new_finish
end

""" Rebuild mothership solution after a cross-cluster perturbation."""
function _rebuild_mothership_solution(
    soln,
    new_clusters,
    exclusions_mothership
)::Tuple{DataFrame,MothershipSolution}
    # Update mothership route and waypoints based on updated clusters
    depot::Point{2,Float64} = soln.mothership_routes[end].route.nodes[1]
    cluster_seq_ids::Vector{Int64} = getfield.(soln.tenders[end], :id)
    cluster_centroids::Vector{Point{2,Float64}} = getfield.(new_clusters, :centroid)
    cluster_sequence::DataFrame = get_cluster_sequence_df(
        depot,
        cluster_seq_ids,
        cluster_centroids
    )
    updated_waypoints::DataFrame = get_waypoints(cluster_sequence, exclusions_mothership)

    waypoint_dist_vector, waypoint_path_vector = get_feasible_vector(
        updated_waypoints.waypoint,
        exclusions_mothership
    )

    updated_ms_solution = MothershipSolution(
        cluster_sequence,
        Route(
            updated_waypoints.waypoint,
            waypoint_dist_vector,
            vcat(waypoint_path_vector...)
        )
    )
    return updated_waypoints, updated_ms_solution
end

""" Finalise modified tender solutions for cross-cluster perturbations by rebuilding
    modified sorties with new waypoints and swapped nodes."""
function _finalise_cross_cluster_tenders(
    sorties_a::Vector{Route},
    sorties_b::Vector{Route},
    source_nodes::Vector{Point{2,Float64}},
    dest_nodes::Vector{Point{2,Float64}},
    tender_a_new_start::Point{2,Float64},
    tender_a_new_finish::Point{2,Float64},
    tender_b_new_start::Point{2,Float64},
    tender_b_new_finish::Point{2,Float64},
    exclusions_tender::POLY_VEC,
    cluster_a_idx::Int64,
    cluster_b_idx::Int64,
    source_sortie_idx::Int64,
    dest_sortie_idx::Int64
)::Tuple{TenderSolution,TenderSolution}
    sorties_a[source_sortie_idx], sorties_b[dest_sortie_idx] = _build_sortie_route.(
        [source_nodes, dest_nodes],
        [tender_a_new_start, tender_b_new_start],
        [tender_a_new_finish, tender_b_new_finish],
        Ref(exclusions_tender)
    )

    # Rebuild tender solutions, preserving cluster identity but updating sorties to reflect:
    # (1) swapped node membership, and (2) updated mothership launch/rendezvous waypoints
    tenders_a_new, tenders_b_new = TenderSolution.(
        [cluster_a_idx, cluster_b_idx],
        [tender_a_new_start, tender_b_new_start],
        [tender_a_new_finish, tender_b_new_finish],
        [sorties_a, sorties_b]
    )

    # Waypoints must not appear as interior sortie nodes
    (tender_a_new_start ∈ vcat(getfield.(sorties_a, :nodes)...) ||
     tender_b_new_start ∈ vcat(getfield.(sorties_b, :nodes)...) ||
     tender_a_new_finish ∈ vcat(getfield.(sorties_a, :nodes)...) ||
     tender_b_new_finish ∈ vcat(getfield.(sorties_b, :nodes)...)) &&
        throw(ArgumentError(
            "Tender start/end point wrongly included in sortie nodes after perturbation"
        ))
    return tenders_a_new, tenders_b_new
end

"""
    perturb_swap(
        soln::MSTSolution,
        clust_seq_idx::Int64,
        exclusions_tender::POLY_VEC=POLY_VEC();
        enforce_diff_sortie::Bool=false
    )::MSTSolution
    perturb_swap(
        soln::MSTSolution,
        cluster_pair::Tuple{Int,Int},
        problem::Problem,
    )::MSTSolution

Perturb the solution by swapping two nodes:
- within a cluster, or
- between two clusters.

# Arguments
- `soln`: Solution to perturb.
- `clust_seq_idx`: Sequence index of cluster to perturb.
- `cluster_pair`: Tuple of two cluster sequence indices to swap nodes between.
- `problem`: Problem instance used to access exclusion zones and other problem parameters.
- `exclusions_tender`: Exclusion zone polygons for tender. Default = DataFrame().
- `enforce_diff_sortie`: Boolean flag to enforce different sorties for node swaps.

# Returns
Perturbed full solution.
"""
function perturb_swap(
    soln::MSTSolution,
    clust_seq_idx::Int64,
    exclusions_tender::POLY_VEC=POLY_VEC();
    enforce_diff_sortie::Bool=false
)::MSTSolution
    tender = deepcopy(soln.tenders[end][clust_seq_idx])
    sorties = deepcopy(tender.sorties)

    # If < 2 sorties in cluster, no perturbation possible
    no_sorties = length(sorties)
    no_sorties < 2 && return soln

    # Choose two random sorties to swap nodes between
    # if enforce_diff_sortie is true: ensure different sorties are selected
    sortie_a_idx, sortie_b_idx = enforce_diff_sortie ?
                                 shuffle(1:no_sorties)[1:2] :
                                 rand(1:no_sorties, 2)

    sortie_a, sortie_b = sorties[sortie_a_idx], sorties[sortie_b_idx]

    # No perturbation possible if a sortie has no nodes
    isempty(sortie_a.nodes) || isempty(sortie_b.nodes) && return soln

    # Determine node indices to swap
    if sortie_a_idx == sortie_b_idx
        # Ensure two distinct node indices are selected from same sortie
        sortie_length = length(sortie_a.nodes)
        if sortie_length < 2
            return soln
        end
        node_a_idx, node_b_idx = sortie_length == 2 ? (1, 2) : shuffle(1:sortie_length)[1:2]
    else
        # Chose two random node indices from different sorties
        node_a_idx = rand(1:length(sortie_a.nodes))
        node_b_idx = rand(1:length(sortie_b.nodes))
    end

    # Swap the nodes between the two sorties
    node_a, node_b = sortie_a.nodes[node_a_idx], sortie_b.nodes[node_b_idx]
    sortie_a.nodes[node_a_idx] = node_b
    sortie_b.nodes[node_b_idx] = node_a

    # TODO:Re-run two-opt on the modified sorties
    # Recompute the feasible paths for the modified sorties
    updated_tender_tours::Vector{Vector{Point{2,Float64}}} = [
        [[tender.start]; sortie_a.nodes; [tender.finish]],
        [[tender.start]; sortie_b.nodes; [tender.finish]]
    ]

    # Get new linestrings
    linestrings::Vector{Vector{Vector{LineString{2,Float64}}}} = getindex.(
        get_feasible_vector.(updated_tender_tours, Ref(exclusions_tender)),
        2
    )

    # Get the new feasible distance matrices for the modified sorties
    updated_tender_matrices::Vector{Matrix{Float64}} = getindex.(get_feasible_matrix.(
            updated_tender_tours,
            Ref(exclusions_tender)
        ), 1)

    # Update the modified sorties
    sorties[sortie_a_idx] = Route(
        sortie_a.nodes,
        updated_tender_matrices[1],
        vcat(linestrings[1]...)
    )
    sorties[sortie_b_idx] = Route(
        sortie_b.nodes,
        updated_tender_matrices[2],
        vcat(linestrings[2]...)
    )

    tender_new = TenderSolution(tender, sorties)
    tender_improved = two_opt(
        tender_new,
        exclusions_tender
    )

    tenders_all::Vector{TenderSolution} = copy(soln.tenders[end])
    tenders_all[clust_seq_idx] = tender_improved

    return MSTSolution(soln.cluster_sets, soln.mothership_routes, [tenders_all])
end
function perturb_swap(
    soln::MSTSolution,
    cluster_pair::Tuple{Int,Int},
    problem::Problem,
)::MSTSolution
    exclusions_mothership::POLY_VEC = problem.mothership.exclusion.geometry
    exclusions_tender::POLY_VEC = problem.tenders.exclusion.geometry

    clust_a_seq_idx, clust_b_seq_idx = cluster_pair

    tender_a::TenderSolution = soln.tenders[end][clust_a_seq_idx]
    tender_b::TenderSolution = soln.tenders[end][clust_b_seq_idx]

    # Pick random sorties and ensure both have nodes
    sortie_a_idx = rand(eachindex(tender_a.sorties))
    sortie_b_idx = rand(eachindex(tender_b.sorties))

    # Assert chosen sorties not empty
    @assert(!isempty(tender_a.sorties[sortie_a_idx].nodes),
        "Empty sortie in tender_a during cross-cluster swap")
    @assert(!isempty(tender_b.sorties[sortie_b_idx].nodes),
        "Empty sortie in tender_b during cross-cluster swap")

    sortie_a_nodes = copy(tender_a.sorties[sortie_a_idx].nodes)
    sortie_b_nodes = copy(tender_b.sorties[sortie_b_idx].nodes)

    # Pick two random nodes across the two sorties
    node_a_idx, node_b_idx = rand(1:length(sortie_a_nodes)), rand(1:length(sortie_b_nodes))
    node_a, node_b = sortie_a_nodes[node_a_idx], sortie_b_nodes[node_b_idx]

    # Swap the nodes between the two sorties
    sortie_a_nodes[node_a_idx] = node_b
    sortie_b_nodes[node_b_idx] = node_a

    # Update new clusters
    new_clusters::Vector{Cluster} = copy(soln.cluster_sets[end])

    cluster_a_idx::Int, cluster_b_idx::Int = tender_a.id, tender_b.id

    nodes_a::Vector{Point{2,Float64}} = copy(new_clusters[cluster_a_idx].nodes)
    nodes_b::Vector{Point{2,Float64}} = copy(new_clusters[cluster_b_idx].nodes)

    # Swap nodes between clusters
    node_a_idx_clust::Union{Int,Nothing} = findfirst(isequal(node_a), nodes_a)
    node_b_idx_clust::Union{Int,Nothing} = findfirst(isequal(node_b), nodes_b)
    (node_a_idx_clust === nothing || node_b_idx_clust === nothing) &&
        throw(ArgumentError("Swap node not found in respective cluster nodes"))

    nodes_a[node_a_idx_clust] = node_b
    nodes_b[node_b_idx_clust] = node_a

    centroid_a = Point{2,Float64}(mean(getindex.(nodes_a, 1)), mean(getindex.(nodes_a, 2)))
    centroid_b = Point{2,Float64}(mean(getindex.(nodes_b, 1)), mean(getindex.(nodes_b, 2)))

    new_clusters[cluster_a_idx] = Cluster(
        id=cluster_a_idx,
        centroid=centroid_a,
        nodes=nodes_a
    )
    new_clusters[cluster_b_idx] = Cluster(
        id=cluster_b_idx,
        centroid=centroid_b,
        nodes=nodes_b
    )

    # Update mothership route and waypoints based on updated clusters
    updated_waypoints, updated_ms_solution = _rebuild_mothership_solution(
        soln,
        new_clusters,
        exclusions_mothership
    )

    # Update tender solutions with the new (start/finish) waypoints based on perturbation
    tenders_all::Vector{TenderSolution} = generate_tender_sorties(
        MSTSolution([new_clusters], [updated_ms_solution], [soln.tenders[end]]),
        updated_waypoints.waypoint,
        exclusions_tender
    )

    tender_a_new_start, tender_a_new_finish, tender_b_new_start, tender_b_new_finish =
        _resolve_tender_endpoints(updated_waypoints, cluster_a_idx, cluster_b_idx)

    # Build new routes for modified sorties based on updated waypoints and swapped nodes
    tenders_a_new, tenders_b_new = _finalise_cross_cluster_tenders(
        copy(tenders_all[clust_a_seq_idx].sorties),
        copy(tenders_all[clust_b_seq_idx].sorties),
        sortie_a_nodes,
        sortie_b_nodes,
        tender_a_new_start,
        tender_a_new_finish,
        tender_b_new_start,
        tender_b_new_finish,
        exclusions_tender,
        cluster_a_idx,
        cluster_b_idx,
        sortie_a_idx,
        sortie_b_idx
    )

    # Re-run two-opt on the modified tender solutions
    tenders_all[clust_a_seq_idx], tenders_all[clust_b_seq_idx] = two_opt.(
        [tenders_a_new, tenders_b_new],
        Ref(exclusions_tender)
    )

    return MSTSolution([new_clusters], [updated_ms_solution], [tenders_all])
end

""" Finds longest sortie in a cluster and determines if it's on the critical path."""
function _find_longest_feasible_sortie(
    soln::MSTSolution,
    clust_seq_idx::Int64,
)::Tuple{Int,Vector{Point{2,Float64}}}
    #! Find THE LONGEST sortie with 2+ nodes: so 1 can be moved and remain feasible
    sortie_lengths = tender_clust_dist(soln.tenders[end][clust_seq_idx])
    longest_sortie_idx = argmax(sortie_lengths)
    longest_sortie_nodes = soln.tenders[end][clust_seq_idx].sorties[longest_sortie_idx].nodes

    return longest_sortie_idx, longest_sortie_nodes
end

function _any_feasible_sortie(
    sorties::Vector{Route},
)::Tuple{Int,Vector{Point{2,Float64}}}
    sortie_lengths = length.(getfield.(sorties, :nodes))

    # Filter sorties with length > 2 so a move is feasibile
    feasible_sorties = findall(x -> x >= 2, sortie_lengths)
    isempty(feasible_sorties) && return 0, Vector{Point{2,Float64}}()

    # Return a random feasible sortie
    sortie_idx = rand(feasible_sorties)
    sortie_nodes = copy(sorties[sortie_idx].nodes)
    return sortie_idx, sortie_nodes
end

""" Moves a random node from source to destination node list for perturb_move."""
function _move_node!(
    source_nodes::Vector{Point{2,Float64}},
    dest_nodes::Vector{Point{2,Float64}}
)::Point{2,Float64}
    node_idx = rand(1:length(source_nodes))
    node = source_nodes[node_idx]
    deleteat!(source_nodes, node_idx)
    insert_pos = rand(1:(length(dest_nodes)+1))
    insert!(dest_nodes, insert_pos, node)
    return node
end

"""
    perturb_move(
        soln::MSTSolution,
        clust_seq_idx::Int64,
        problem::Problem,
    )::MSTSolution
    perturb_move(
        soln::MSTSolution,
        cluster_pair::Tuple{Int,Int},
        problem::Problem,
    )::MSTSolution

Perturb the solution by moving a single node from one sortie to another:
- within a cluster, or
- between two clusters.

# Arguments
- `soln`: Solution to perturb.
- `clust_seq_idx`: Sequence index of cluster to perturb.
- `cluster_pair`: Tuple of two cluster sequence indices to move a node between.
- `problem`: Problem instance used to access exclusion zones and other problem parameters.

# Returns
Perturbed full solution.
"""
function perturb_move(
    soln::MSTSolution,
    clust_seq_idx::Int64,
    problem::Problem,
)::MSTSolution
    #! WITHIN CLUSTER MOVE
    exclusions_tender::POLY_VEC = problem.tenders.exclusion.geometry

    tender = deepcopy(soln.tenders[end][clust_seq_idx])
    sorties = deepcopy(tender.sorties)

    # No perturbation possible if < 2 sorties in cluster
    no_sorties::Int = length(sorties)
    no_sorties < 2 && return soln

    # Choose a random feasible (2+ nodes) sortie, and get its nodes
    source_idx, source_nodes = _any_feasible_sortie(
        sorties,
    )

    # If no feasible sortie found, exit
    source_idx == 0 && return soln

    dest_idx::Int = rand(setdiff(1:no_sorties, source_idx))
    dest_nodes = copy(sorties[dest_idx].nodes)

    # Remove a random node from source, insert at random position in destination
    _move_node!(source_nodes, dest_nodes)

    # Recompute routes for both modified sorties
    updated_tours::Vector{Vector{Point{2,Float64}}} = [
        [[tender.start]; source_nodes; [tender.finish]],
        [[tender.start]; dest_nodes; [tender.finish]]
    ]

    linestrings::Vector{Vector{Vector{LineString{2,Float64}}}} = getindex.(
        get_feasible_vector.(updated_tours, Ref(exclusions_tender)),
        2
    )
    updated_matrices::Vector{Matrix{Float64}} = getindex.(
        get_feasible_matrix.(updated_tours, Ref(exclusions_tender)),
        1
    )

    sorties[source_idx] = Route(source_nodes, updated_matrices[1], vcat(linestrings[1]...))
    sorties[dest_idx] = Route(dest_nodes, updated_matrices[2], vcat(linestrings[2]...))

    tender_new = TenderSolution(tender, sorties)
    tender_improved = two_opt(tender_new, exclusions_tender)

    tenders_all::Vector{TenderSolution} = copy(soln.tenders[end])
    tenders_all[clust_seq_idx] = tender_improved

    soln_proposed = MSTSolution(soln.cluster_sets, soln.mothership_routes, [tenders_all])
    return soln_proposed
end
function perturb_move(
    soln::MSTSolution,
    cluster_pair::Tuple{Int,Int},
    problem::Problem,
)::MSTSolution
    exclusions_mothership::POLY_VEC = problem.mothership.exclusion.geometry
    exclusions_tender::POLY_VEC = problem.tenders.exclusion.geometry

    clust_a_seq_idx, clust_b_seq_idx = cluster_pair

    tender_a::TenderSolution = soln.tenders[end][clust_a_seq_idx]
    tender_b::TenderSolution = soln.tenders[end][clust_b_seq_idx]

    # Choose a random feasible (2+ nodes) sortie, and get its nodes
    source_sortie_idx, source_nodes = _any_feasible_sortie(
        tender_a.sorties,
    )

    # If no feasible sortie found, exit
    source_sortie_idx == 0 && return soln

    dest_sortie_idx::Int = rand(eachindex(tender_b.sorties))
    dest_nodes = copy(tender_b.sorties[dest_sortie_idx].nodes)

    # Remove a random node from source, insert at random position in destination
    node::Point{2,Float64} = _move_node!(source_nodes, dest_nodes)

    # Update cluster membership
    new_clusters::Vector{Cluster} = copy(soln.cluster_sets[end])
    cluster_a_idx::Int = tender_a.id
    cluster_b_idx::Int = tender_b.id

    nodes_a::Vector{Point{2,Float64}} = copy(new_clusters[cluster_a_idx].nodes)
    nodes_b::Vector{Point{2,Float64}} = copy(new_clusters[cluster_b_idx].nodes)

    node_a_idx_clust::Union{Int,Nothing} = findfirst(isequal(node), nodes_a)
    node_a_idx_clust === nothing &&
        throw(ArgumentError("Move node not found in source cluster nodes"))

    deleteat!(nodes_a, node_a_idx_clust)
    push!(nodes_b, node)

    centroid_a = Point{2,Float64}(mean(getindex.(nodes_a, 1)), mean(getindex.(nodes_a, 2)))
    centroid_b = Point{2,Float64}(mean(getindex.(nodes_b, 1)), mean(getindex.(nodes_b, 2)))

    new_clusters[cluster_a_idx] = Cluster(id=cluster_a_idx, centroid=centroid_a, nodes=nodes_a)
    new_clusters[cluster_b_idx] = Cluster(id=cluster_b_idx, centroid=centroid_b, nodes=nodes_b)

    # Update mothership route and waypoints based on updated clusters
    updated_waypoints, updated_ms_solution = _rebuild_mothership_solution(
        soln,
        new_clusters,
        exclusions_mothership
    )

    # Update tender solutions with the new (start/finish) waypoints based on perturbation
    tenders_all::Vector{TenderSolution} = generate_tender_sorties(
        MSTSolution([new_clusters], [updated_ms_solution], [soln.tenders[end]]),
        updated_waypoints.waypoint,
        exclusions_tender
    )

    tender_a_new_start, tender_a_new_finish, tender_b_new_start, tender_b_new_finish =
        _resolve_tender_endpoints(updated_waypoints, cluster_a_idx, cluster_b_idx)

    # Build new routes for modified sorties based on updated waypoints and swapped nodes
    tenders_a_new, tenders_b_new = _finalise_cross_cluster_tenders(
        copy(tenders_all[clust_a_seq_idx].sorties),
        copy(tenders_all[clust_b_seq_idx].sorties),
        source_nodes,
        dest_nodes,
        tender_a_new_start,
        tender_a_new_finish,
        tender_b_new_start,
        tender_b_new_finish,
        exclusions_tender,
        cluster_a_idx,
        cluster_b_idx,
        source_sortie_idx,
        dest_sortie_idx
    )

    tenders_all[clust_a_seq_idx], tenders_all[clust_b_seq_idx] = two_opt.(
        [tenders_a_new, tenders_b_new],
        Ref(exclusions_tender)
    )

    soln_proposed = MSTSolution([new_clusters], [updated_ms_solution], [tenders_all])
    return soln_proposed
end

"""
    simulated_annealing(
        problem::Problem,
        soln_init::MSTSolution,
        objective_function::Function,
        max_iterations::Int,
        temp_init::Float64,
        cooling_rate::Float64,
        min_iters::Int,
        static_limit::Int;
        cross_cluster_flag::Bool,
    )::Tuple{MSTSolution,Float64}

Simulated Annealing optimization algorithm to optimize the solution.

# Arguments
- `problem`: Problem instance used to access exclusion zones and other problem parameters.
- `soln_init`: Initial solution.
- `objective_function`: Function to evaluate the solution.
- `max_iterations`: Maximum number of iterations.
- `temp_init`: Initial temperature.
- `cooling_rate`: Rate of cooling to guide acceptance probability for SA algorithm.
- `min_iters`: Minimum number of iterations to perform before allowing early exit.
- `static_limit`: Number of iterations to allow stagnation before early exit.
- `cross_cluster_flag`: Boolean to indicate consideration of cross-cluster perturbations.

# Returns
- `soln_best`: Best solution::MSTSolution found.
- `obj_best`: Objective value of the best solution.
"""
function simulated_annealing(
    problem::Problem,
    soln_init::MSTSolution,
    objective_function::Function,
    max_iterations::Int,
    temp_init::Float64,
    cooling_rate::Float64,
    min_iters::Int,
    static_limit::Int;
    cross_cluster_flag::Bool,
)::Tuple{MSTSolution,Float64}
    exclusions_tender::POLY_VEC = problem.tenders.exclusion.geometry
    vessel_weightings::NTuple{2,AbstractFloat} = (
        problem.mothership.weighting,
        problem.tenders.weighting
    )

    @info "Starting Simulated Annealing Optimization with parameters:"
    @info "  - Max Iterations: $max_iterations"
    @info "  - Initial Temperature: $temp_init"
    @info "  - Cooling Rate: $cooling_rate"
    @info "  - Minimum Iterations: $min_iters"
    @info "  - Static Limit: $static_limit"
    @info "  - Cross-Cluster Flag: $cross_cluster_flag"

    # Initialize best solution as initial
    soln_best = deepcopy(soln_init)
    soln_current = deepcopy(soln_best)
    soln_proposed = MSTSolution([], [], [])
    obj_best = objective_function(soln_init, vessel_weightings)
    obj_init = obj_best
    obj_current = obj_best
    obj_proposed = obj_best
    cluster_set::Vector{Cluster} = soln_init.cluster_sets[end]
    no_clusts::Int = length(cluster_set)

    # Initialize current solution as best, reset temp
    temp = temp_init
    static_ctr = 0

    # Trace vectors for logging
    trace_iters = Int[]
    trace_best = Float64[]
    trace_current = Float64[]
    trace_proposed = Float64[]
    trace_temps = Float64[]

    clust_idx = Int(0)
    clust_alt_idx = Int(0)
    shuffled_clusters = Vector{Int}(undef, no_clusts)
    perturbation_type::Symbol = :none

    @info "Iter\t| Perturbation\t| Best\t\t| Current\t| Proposed\t| Temp\t"

    for iteration in 1:max_iterations
        shuffled_clusters = shuffle(1:no_clusts)
        clust_idx = shuffled_clusters[1]
        clust_alt_idx = 0

        if !cross_cluster_flag || rand() < 0.5
            # SUB-cluster perturbation
            # Swap/move within the same cluster @ 50/50
            if rand() < 0.5
                # Swap 2 nodes across sorties within the same cluster
                soln_proposed = perturb_swap(soln_current, clust_idx, exclusions_tender)
                perturbation_type = :SWAP
            else
                # Move a node across sorties within the same cluster
                soln_proposed = perturb_move(soln_current, clust_idx, problem)
                perturbation_type = :MOVE
            end
        else
            # CROSS-cluster perturbation
            clust_alt_idx = shuffled_clusters[2]

            if rand() < 0.5
                # Swap 2 nodes between a sortie in one cluster to another
                soln_proposed = perturb_swap(
                    soln_current,
                    (clust_idx, clust_alt_idx),
                    problem
                )
                perturbation_type = :SWAP
            else
                # Move a node from a sortie in one cluster to another
                soln_proposed = perturb_move(
                    soln_current,
                    (clust_idx, clust_alt_idx),
                    problem
                )
                perturbation_type = :MOVE
            end
        end

        static_ctr += 1
        obj_proposed = objective_function(soln_proposed, vessel_weightings)
        Δ = obj_proposed - obj_current

        if Δ == 0
            Δ = total_distance(soln_proposed, vessel_weightings) -
                total_distance(soln_current, vessel_weightings)
        end

        # If the new solution is improved OR meets acceptance prob criteria
        if Δ < 0 || rand() < exp(-Δ / temp)
            soln_current = soln_proposed
            obj_current = obj_proposed

            if obj_current < obj_best
                static_ctr = 0
                soln_best = deepcopy(soln_current)
                obj_best = obj_current
            end
        end

        perturbed_clusters = clust_alt_idx == 0 ? " $clust_idx" : "$clust_idx -> $clust_alt_idx"
        @info "$iteration\t| $perturbation_type: $perturbed_clusters\t| $(round(obj_best, digits=4))\t| $(round(obj_current, digits=4))\t| $(round(obj_proposed, digits=4))\t| $(round(temp, digits=4))\t"

        # Record after accept/reject so current/best reflect the updated state
        push!(trace_iters, iteration)
        push!(trace_proposed, obj_proposed)
        push!(trace_current, obj_current)
        push!(trace_best, obj_best)
        push!(trace_temps, temp)

        temp *= cooling_rate

        if iteration >= min_iters && static_ctr >= static_limit
            @info "$iteration\t| $perturbation_type: $perturbed_clusters\t| $(round(obj_best, digits=4))\t| $(round(obj_current, digits=4))\t| $(round(obj_proposed, digits=4))\t| $(round(temp, digits=4))\t"

            # Display trace after each cluster's SA run completes
            display(Plot.sa_trace(
                trace_iters, trace_best, trace_current, trace_proposed, trace_temps
            ))
            break
        end
    end

    @info "Final Value:\t$obj_best\nΔ: $(obj_init - obj_best)"
    return soln_best, obj_best
end
