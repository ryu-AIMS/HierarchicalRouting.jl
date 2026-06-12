
"""
    _build_sortie_route(
        sortie_nodes::Vector{Point{2,Float64}},
        start_pt::Point{2,Float64},
        finish_pt::Point{2,Float64},
        exclusions::POLY_VEC
    )::Route

Compute a feasible `Route` for a sortie given its interior nodes and start/finish waypoints.
Used when sortie node membership changes.
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
    ms_nodes::Vector{Point{2,Float64}},
    connecting_clusters::Vector{NTuple{2,Int64}},
    cluster_a_idx::Int,
    cluster_b_idx::Int
)::NTuple{4,Point}
    cc_first, cc_last = first.(connecting_clusters), last.(connecting_clusters)

    tender_a_new_start = ms_nodes[findlast(cc_last .== cluster_a_idx)]
    tender_a_new_finish = ms_nodes[findfirst(cc_first .== cluster_a_idx)]
    tender_b_new_start = ms_nodes[findlast(cc_last .== cluster_b_idx)]
    tender_b_new_finish = ms_nodes[findfirst(cc_first .== cluster_b_idx)]

    return tender_a_new_start, tender_a_new_finish, tender_b_new_start, tender_b_new_finish
end

""" Rebuild mothership solution after a cross-cluster perturbation."""
function _rebuild_mothership_solution(
    soln,
    new_clusters,
    exclusions,
    perturb_idxs::UnitRange{Int}
)::Tuple{DataFrame,MothershipSolution}
    # Use FULL cluster sequence from existing MS route
    existing_seq = soln.mothership_routes[end].cluster_sequence
    cluster_seq_ids = filter(!=(0), existing_seq.id)

    # Merge existing centroids as base, new_clusters to override modified ones
    existing_centroid_map = Dict(
        row.id => Point{2,Float64}(row.lon, row.lat)
        for row in eachrow(existing_seq) if row.id != 0
    )
    updated_centroid_map = Dict(c.id => c.centroid for c in new_clusters)
    centroid_map = merge(existing_centroid_map, updated_centroid_map)

    ordered_centroids = [centroid_map[id] for id in cluster_seq_ids]
    depot::Point{2,Float64} = soln.mothership_routes[end].route.nodes[1]
    full_pts = [[depot]; ordered_centroids; [depot]]
    cluster_sequence = DataFrame(
        id=[0; cluster_seq_ids; 0],
        lon=getindex.(full_pts, 1),
        lat=getindex.(full_pts, 2),
    )

    updated_waypoints::DataFrame = get_waypoints(cluster_sequence, exclusions)

    n_fixed = 2 * (first(perturb_idxs) - 1) + 1  # depot + 2 per visited cluster

    final_waypoints::Vector{Point{2,Float64}} = updated_waypoints.waypoint
    if n_fixed > 1
        final_waypoints = vcat(
            soln.mothership_routes[end].route.nodes[1:n_fixed],
            updated_waypoints.waypoint[n_fixed+1:end]
        )
    end
    waypoint_dist_vector, waypoint_path_vector = get_feasible_vector(
        final_waypoints,
        exclusions
    )

    updated_ms_solution = MothershipSolution(
        cluster_sequence,
        Route(
            final_waypoints,
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
        rng::AbstractRNG,
        exclusions_tender::POLY_VEC=POLY_VEC();
        enforce_diff_sortie::Bool=false
    )::MSTSolution
    perturb_swap(
        soln::MSTSolution,
        cluster_pair::Tuple{Int,Int},
        problem::Problem,
        perturb_idxs::UnitRange{Int},
        rng::AbstractRNG,
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
- `rng`: Random number generator.
- `enforce_diff_sortie`: Boolean flag to enforce different sorties for node swaps.

# Returns
Perturbed full solution.
"""
function perturb_swap(
    soln::MSTSolution,
    clust_seq_idx::Int64,
    rng::AbstractRNG,
    exclusions_tender::POLY_VEC=POLY_VEC();
    enforce_diff_sortie::Bool=false
)::MSTSolution
    #! WITHIN CLUSTER SWAP
    tender = deepcopy(soln.tenders[end][clust_seq_idx])
    sorties = deepcopy(tender.sorties)

    # If < 2 sorties in cluster, no perturbation possible
    no_sorties = length(sorties)
    no_sorties < 2 && return soln

    # Choose two random sorties to swap nodes between
    # if enforce_diff_sortie is true: ensure different sorties are selected
    sortie_a_idx, sortie_b_idx = enforce_diff_sortie ?
                                 shuffle(rng, 1:no_sorties)[1:2] :
                                 rand(rng, 1:no_sorties, 2)

    sortie_a, sortie_b = sorties[sortie_a_idx], sorties[sortie_b_idx]

    # No perturbation possible if a sortie has no nodes
    (isempty(sortie_a.nodes) || isempty(sortie_b.nodes)) && return soln

    # Determine node indices to swap
    if sortie_a_idx == sortie_b_idx
        # Ensure two distinct node indices are selected from same sortie
        sortie_length = length(sortie_a.nodes)
        if sortie_length < 2
            return soln
        end
        node_a_idx, node_b_idx =
            sortie_length == 2 ?
            (1, 2) :
            shuffle(rng, 1:sortie_length)[1:2]
    else
        # Chose two random node indices from different sorties
        node_a_idx = rand(rng, 1:length(sortie_a.nodes))
        node_b_idx = rand(rng, 1:length(sortie_b.nodes))
    end

    # Swap the nodes between the two sorties
    node_a, node_b = sortie_a.nodes[node_a_idx], sortie_b.nodes[node_b_idx]
    sortie_a.nodes[node_a_idx] = node_b
    sortie_b.nodes[node_b_idx] = node_a

    return _finalise_within_cluster_soln(
        soln,
        clust_seq_idx,
        tender,
        sorties,
        sortie_a_idx,
        sortie_a.nodes,
        sortie_b_idx,
        sortie_b.nodes,
        exclusions_tender
    )
end
function perturb_swap(
    soln::MSTSolution,
    cluster_pair::Tuple{Int,Int},
    problem::Problem,
    perturb_idxs::UnitRange{Int},
    rng::AbstractRNG,
)::MSTSolution
    #! CROSS-CLUSTER SWAP
    clust_a_seq_idx, clust_b_seq_idx = cluster_pair

    tender_a::TenderSolution = soln.tenders[end][clust_a_seq_idx]
    tender_b::TenderSolution = soln.tenders[end][clust_b_seq_idx]

    # Pick random sorties and ensure both have nodes
    sortie_a_idx = rand(rng, eachindex(tender_a.sorties))
    sortie_b_idx = rand(rng, eachindex(tender_b.sorties))

    # Ensure chosen sorties not empty
    (isempty(tender_a.sorties[sortie_a_idx].nodes) ||
     isempty(tender_b.sorties[sortie_b_idx].nodes)) &&
        return soln

    sortie_a_nodes = copy(tender_a.sorties[sortie_a_idx].nodes)
    sortie_b_nodes = copy(tender_b.sorties[sortie_b_idx].nodes)

    # Pick two random nodes across the two sorties
    node_a_idx, node_b_idx = rand(rng, 1:length(sortie_a_nodes)), rand(rng, 1:length(sortie_b_nodes))
    node_a, node_b = sortie_a_nodes[node_a_idx], sortie_b_nodes[node_b_idx]

    # Swap the nodes between the two sorties
    sortie_a_nodes[node_a_idx] = node_b
    sortie_b_nodes[node_b_idx] = node_a

    # Update new clusters
    new_clusters::Vector{Cluster} = copy(soln.cluster_sets[end])

    cluster_a_idx::Int, cluster_b_idx::Int = tender_a.id, tender_b.id

    new_clust_a_pos::Int = findfirst(c -> c.id == cluster_a_idx, new_clusters)
    new_clust_b_pos::Int = findfirst(c -> c.id == cluster_b_idx, new_clusters)
    nodes_a::Vector{Point{2,Float64}} = copy(new_clusters[new_clust_a_pos].nodes)
    nodes_b::Vector{Point{2,Float64}} = copy(new_clusters[new_clust_b_pos].nodes)

    # Swap nodes between clusters
    node_a_idx_clust::Union{Int,Nothing} = findfirst(isequal(node_a), nodes_a)
    node_b_idx_clust::Union{Int,Nothing} = findfirst(isequal(node_b), nodes_b)
    (node_a_idx_clust === nothing || node_b_idx_clust === nothing) &&
        throw(ArgumentError("Swap node not found in respective cluster nodes"))

    nodes_a[node_a_idx_clust] = node_b
    nodes_b[node_b_idx_clust] = node_a

    centroid_a = Point{2,Float64}(mean(getindex.(nodes_a, 1)), mean(getindex.(nodes_a, 2)))
    centroid_b = Point{2,Float64}(mean(getindex.(nodes_b, 1)), mean(getindex.(nodes_b, 2)))

    new_clusters[new_clust_a_pos] = Cluster(
        id=cluster_a_idx,
        centroid=centroid_a,
        nodes=nodes_a
    )
    new_clusters[new_clust_b_pos] = Cluster(
        id=cluster_b_idx,
        centroid=centroid_b,
        nodes=nodes_b
    )

    return _apply_cross_cluster_perturbation(
        soln,
        new_clusters,
        problem,
        clust_a_seq_idx,
        clust_b_seq_idx,
        cluster_a_idx,
        cluster_b_idx,
        sortie_a_nodes,
        sortie_b_nodes,
        sortie_a_idx,
        sortie_b_idx,
        perturb_idxs,
    )
end

""" Finds longest sortie in a cluster and determines if it's on the critical path."""
function _find_longest_feasible_sortie(
    soln::MSTSolution,
    clust_seq_idx::Int64,
)::Tuple{Int,Vector{Point{2,Float64}}}
    sortie_lengths = tender_clust_dist(soln.tenders[end][clust_seq_idx])
    longest_sortie_idx = argmax(sortie_lengths)
    longest_sortie_nodes =
        soln.tenders[end][clust_seq_idx].sorties[longest_sortie_idx].nodes

    return longest_sortie_idx, longest_sortie_nodes
end

function _any_feasible_sortie(
    sorties::Vector{Route},
    rng::AbstractRNG,
)::Tuple{Int,Vector{Point{2,Float64}}}
    sortie_lengths = length.(getfield.(sorties, :nodes))

    # Filter sorties with length > 2 so a move is feasibile
    feasible_sorties = findall(x -> x >= 2, sortie_lengths)
    isempty(feasible_sorties) && return 0, Vector{Point{2,Float64}}()

    # Return a random feasible sortie
    sortie_idx = rand(rng, feasible_sorties)
    sortie_nodes = copy(sorties[sortie_idx].nodes)
    return sortie_idx, sortie_nodes
end

""" Moves a random node from source to destination node list for perturb_move."""
function _move_node!(
    source_nodes::Vector{Point{2,Float64}},
    dest_nodes::Vector{Point{2,Float64}},
    rng::AbstractRNG,
)::Point{2,Float64}
    node_idx = rand(rng, 1:length(source_nodes))
    node = source_nodes[node_idx]
    deleteat!(source_nodes, node_idx)
    insert_pos = rand(rng, 1:(length(dest_nodes)+1))
    insert!(dest_nodes, insert_pos, node)
    return node
end

"""
    perturb_move(
        soln::MSTSolution,
        clust_seq_idx::Int64,
        problem::Problem,
        rng::AbstractRNG,
    )::MSTSolution
    perturb_move(
        soln::MSTSolution,
        cluster_pair::Tuple{Int,Int},
        problem::Problem,
        rng::AbstractRNG,
    )::MSTSolution

Perturb the solution by moving a single node from one sortie to another:
- within a cluster, or
- between two clusters.

# Arguments
- `soln`: Solution to perturb.
- `clust_seq_idx`: Sequence index of cluster to perturb.
- `cluster_pair`: Tuple of two cluster sequence indices to move a node between.
- `problem`: Problem instance used to access exclusion zones and other problem parameters.
- `rng`: Random number generator.

# Returns
Perturbed full solution.
"""
function perturb_move(
    soln::MSTSolution,
    clust_seq_idx::Int64,
    problem::Problem,
    rng::AbstractRNG,
)::MSTSolution
    #! WITHIN CLUSTER MOVE
    exclusions_tender::POLY_VEC = problem.tenders.exclusion.geometry

    tender = deepcopy(soln.tenders[end][clust_seq_idx])
    sorties = deepcopy(tender.sorties)

    # No perturbation possible if < 2 sorties in cluster
    no_sorties::Int = length(sorties)
    no_sorties < 2 && return soln

    # Choose a random feasible (2+ nodes) sortie, and get its nodes
    source_idx, source_nodes = _any_feasible_sortie(sorties, rng)

    # If no feasible sortie found, exit
    source_idx == 0 && return soln

    dest_idx::Int = rand(rng, setdiff(1:no_sorties, source_idx))
    dest_nodes = copy(sorties[dest_idx].nodes)

    # Reject move if destination sortie is already at capacity
    length(dest_nodes) >= problem.tenders.capacity && return soln

    # Remove a random node from source, insert at random position in destination
    _move_node!(source_nodes, dest_nodes, rng)

    return _finalise_within_cluster_soln(
        soln,
        clust_seq_idx,
        tender,
        sorties,
        source_idx,
        source_nodes,
        dest_idx,
        dest_nodes,
        exclusions_tender
    )
end
function perturb_move(
    soln::MSTSolution,
    cluster_pair::Tuple{Int,Int},
    problem::Problem,
    perturb_idxs::UnitRange{Int},
    rng::AbstractRNG,
)::MSTSolution
    #! CROSS-CLUSTER MOVE
    clust_a_seq_idx, clust_b_seq_idx = cluster_pair

    tender_a::TenderSolution = soln.tenders[end][clust_a_seq_idx]
    tender_b::TenderSolution = soln.tenders[end][clust_b_seq_idx]

    # Choose a random feasible (2+ nodes) sortie, and get its nodes
    source_sortie_idx, source_nodes = _any_feasible_sortie(
        tender_a.sorties, rng
    )

    # If no feasible sortie found, exit
    source_sortie_idx == 0 && return soln

    dest_sortie_idx::Int = rand(rng, eachindex(tender_b.sorties))
    dest_nodes = copy(tender_b.sorties[dest_sortie_idx].nodes)

    # Reject move if destination sortie is already at capacity
    length(dest_nodes) >= problem.tenders.capacity && return soln

    # Remove a random node from source, insert at random position in destination
    node::Point{2,Float64} = _move_node!(source_nodes, dest_nodes, rng)

    # Update cluster membership
    new_clusters::Vector{Cluster} = copy(soln.cluster_sets[end])
    cluster_a_idx::Int = tender_a.id
    cluster_b_idx::Int = tender_b.id

    new_clust_a_pos::Int = findfirst(c -> c.id == cluster_a_idx, new_clusters)
    new_clust_b_pos::Int = findfirst(c -> c.id == cluster_b_idx, new_clusters)
    nodes_a::Vector{Point{2,Float64}} = copy(new_clusters[new_clust_a_pos].nodes)
    nodes_b::Vector{Point{2,Float64}} = copy(new_clusters[new_clust_b_pos].nodes)

    node_a_idx_clust::Union{Int,Nothing} = findfirst(isequal(node), nodes_a)
    node_a_idx_clust === nothing &&
        throw(ArgumentError("Move node not found in source cluster nodes"))

    deleteat!(nodes_a, node_a_idx_clust)
    push!(nodes_b, node)

    centroid_a = Point{2,Float64}(mean(getindex.(nodes_a, 1)), mean(getindex.(nodes_a, 2)))
    centroid_b = Point{2,Float64}(mean(getindex.(nodes_b, 1)), mean(getindex.(nodes_b, 2)))

    new_clusters[new_clust_a_pos], new_clusters[new_clust_b_pos] =
        Cluster(id=cluster_a_idx, centroid=centroid_a, nodes=nodes_a),
        Cluster(id=cluster_b_idx, centroid=centroid_b, nodes=nodes_b)

    return _apply_cross_cluster_perturbation(
        soln,
        new_clusters,
        problem,
        clust_a_seq_idx,
        clust_b_seq_idx,
        cluster_a_idx,
        cluster_b_idx,
        source_nodes,
        dest_nodes,
        source_sortie_idx,
        dest_sortie_idx,
        perturb_idxs,
    )
end

"""
Rebuild mothership, waypoints, and tender sorties after any cross-cluster perturbation.
Shared body for cross-cluster perturb_move and perturb_swap.
"""
function _apply_cross_cluster_perturbation(
    soln::MSTSolution,
    new_clusters::Vector{Cluster},
    problem::Problem,
    clust_a_seq_idx::Int,
    clust_b_seq_idx::Int,
    cluster_a_idx::Int,
    cluster_b_idx::Int,
    nodes_a::Vector{Point{2,Float64}},
    nodes_b::Vector{Point{2,Float64}},
    sortie_a_idx::Int,
    sortie_b_idx::Int,
    perturb_idxs::UnitRange{Int}
)::MSTSolution
    exclusions_tender = problem.tenders.exclusion.geometry
    exclusions_all = vcat(problem.mothership.exclusion.geometry, exclusions_tender)

    updated_waypoints, updated_ms_solution = _rebuild_mothership_solution(
        soln, new_clusters, exclusions_all, perturb_idxs
    )

    tender_ids = getfield.(soln.tenders[end], :id)
    cc = updated_waypoints.connecting_clusters
    ms_nodes = updated_ms_solution.route.nodes
    partial_waypoints = vcat(
        [ms_nodes[1]],
        [[ms_nodes[findlast(last.(cc) .== id)],
            ms_nodes[findfirst(first.(cc) .== id)]]
         for id in tender_ids]...,
        [ms_nodes[end]]
    )

    tenders_all = generate_tender_sorties(
        MSTSolution([new_clusters], [updated_ms_solution], [soln.tenders[end]]),
        partial_waypoints,
        exclusions_tender
    )

    tender_a_new_start, tender_a_new_finish, tender_b_new_start, tender_b_new_finish =
        _resolve_tender_endpoints(ms_nodes, updated_waypoints.connecting_clusters, cluster_a_idx, cluster_b_idx)

    tenders_a_new, tenders_b_new = _finalise_cross_cluster_tenders(
        copy(tenders_all[clust_a_seq_idx].sorties),
        copy(tenders_all[clust_b_seq_idx].sorties),
        nodes_a,
        nodes_b,
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

    tenders_all[clust_a_seq_idx], tenders_all[clust_b_seq_idx] = two_opt.(
        [tenders_a_new, tenders_b_new], Ref(exclusions_tender)
    )

    return MSTSolution([new_clusters], [updated_ms_solution], [tenders_all])
end

"""
Rebuild two modified sorties, apply two-opt, and return an updated solution.
"""
function _finalise_within_cluster_soln(
    soln::MSTSolution,
    clust_seq_idx::Int,
    tender::TenderSolution,
    sorties::Vector{Route},
    idx_a::Int,
    nodes_a::Vector{Point{2,Float64}},
    idx_b::Int,
    nodes_b::Vector{Point{2,Float64}},
    exclusions::POLY_VEC,
)::MSTSolution
    # Update modified sorties
    sorties[idx_a], sorties[idx_b] = _build_sortie_route.(
        [nodes_a, nodes_b],
        Ref(tender.start),
        Ref(tender.finish),
        Ref(exclusions),
    )

    tender_improved = two_opt(TenderSolution(tender, sorties), exclusions)

    tenders_all = copy(soln.tenders[end])
    tenders_all[clust_seq_idx] = tender_improved
    return MSTSolution(soln.cluster_sets, soln.mothership_routes, [tenders_all])
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
        static_limit::Int,
        rng::AbstractRNG;
        perturb_idxs::UnitRange{Int}=1:length(soln_init.cluster_sets[end]),
        output_dir::String="",
        info_log::Bool,
        plot_flag::Bool,
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
- `rng`: Random number generator for reproducibility.
- `perturb_idxs`: Range of cluster sequence indices to consider for perturbations.
- `output_dir::String`: Path to output directory. If empty, do not save outputs.
- `info_log::Bool`: Flag to switch info statement logging
- `plot_flag`: Flag to plot solution progress for debugging/visualization

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
    static_limit::Int,
    rng::AbstractRNG;
    perturb_idxs::UnitRange{Int}=1:length(soln_init.cluster_sets[end]),
    output_dir::String="",
    info_log::Bool,
    plot_flag::Bool,
)::Tuple{MSTSolution,Float64}
    exclusions_tender::POLY_VEC = problem.tenders.exclusion.geometry
    vessel_weightings::NTuple{2,AbstractFloat} = (
        problem.mothership.weighting,
        problem.tenders.weighting
    )

    if info_log
        @info "Starting Simulated Annealing Optimization with parameters:"
        max_iterations == typemax(Int) || @info "  - Max Iterations: $max_iterations"
        @info "  - Minimum Iterations: $min_iters"
        @info "  - Initial Temperature: $temp_init"
        @info "  - Cooling Rate: $cooling_rate"
        @info "  - Static Limit: $static_limit"
    end

    # Initialize best solution as initial
    soln_best::MSTSolution = deepcopy(soln_init)
    soln_current::MSTSolution = deepcopy(soln_best)
    soln_proposed::MSTSolution = MSTSolution([], [], [])
    obj_best::Float64 = objective_function(soln_init, vessel_weightings)
    obj_init::Float64 = obj_best
    obj_current::Float64 = obj_best
    obj_proposed::Float64 = obj_best

    total_dist_current::Float64 = total_distance(soln_current, vessel_weightings)
    total_dist_proposed::Float64 = total_dist_current
    total_dist_best::Float64 = total_dist_current

    no_clusts::Int = length(perturb_idxs)

    # Initialize current solution as best, reset temp
    temp::Float64 = temp_init
    static_ctr::Int = 0

    # Trace vectors for logging
    trace_iters::Vector{Int} = Int[]
    trace_best::Vector{Float64} = Float64[]
    trace_current::Vector{Float64} = Float64[]
    trace_proposed::Vector{Float64} = Float64[]
    trace_temps::Vector{Float64} = Float64[]

    clust_idx::Int = Int(0)
    clust_alt_idx::Int = Int(0)
    shuffled_clusters::Vector{Int} = Vector{Int}(undef, no_clusts)
    perturbation_type::Symbol = :none

    info_log && @info "Iter\t| Perturbation\t| Best\t\t| Current\t| Proposed\t| Temp\t"

    for iteration in 1:max_iterations
        shuffled_clusters = shuffle(rng, perturb_idxs)
        clust_idx = shuffled_clusters[1]
        clust_alt_idx = 0

        if rand(rng) < 0.5
            # SUB-cluster perturbation
            # Swap/move within the same cluster @ 50/50
            if rand(rng) < 0.5
                # Swap 2 nodes across sorties within the same cluster
                soln_proposed = perturb_swap(
                    soln_current,
                    clust_idx,
                    rng,
                    exclusions_tender,
                )
                perturbation_type = :SWAP
            else
                # Move a node across sorties within the same cluster
                soln_proposed = perturb_move(soln_current, clust_idx, problem, rng)
                perturbation_type = :MOVE
            end
        elseif no_clusts >= 2
            # CROSS-cluster perturbation
            clust_alt_idx = shuffled_clusters[2]

            if rand(rng) < 0.5
                # Swap 2 nodes between a sortie in one cluster to another
                soln_proposed = perturb_swap(
                    soln_current,
                    (clust_idx, clust_alt_idx),
                    problem,
                    perturb_idxs,
                    rng,
                )
                perturbation_type = :SWAP
            else
                # Move a node from a sortie in one cluster to another
                soln_proposed = perturb_move(
                    soln_current,
                    (clust_idx, clust_alt_idx),
                    problem,
                    perturb_idxs,
                    rng
                )
                perturbation_type = :MOVE
            end
        else
            # Only 1 cluster — fall back to within-cluster swap
            soln_proposed = perturb_swap(soln_current, clust_idx, rng, exclusions_tender)
            perturbation_type = :SWAP
        end

        obj_proposed = objective_function(soln_proposed, vessel_weightings)
        Δ = obj_proposed - obj_current

        if Δ == 0
            total_dist_proposed = total_distance(soln_proposed, vessel_weightings)
            Δ = total_dist_proposed - total_dist_current
        end

        perturbed_clusters = clust_alt_idx == 0 ?
                             " $clust_idx" :
                             "$clust_idx -> $clust_alt_idx"

        # If the new solution is improved OR meets acceptance prob criteria
        if Δ < 0 || rand(rng) < exp(-Δ / temp)
            soln_current = soln_proposed
            total_dist_current = obj_proposed == obj_current ?
                                 total_dist_proposed :
                                 total_distance(soln_proposed, vessel_weightings)
            obj_current = obj_proposed

            if obj_current < obj_best ||
               (obj_current == obj_best && total_dist_current < total_dist_best)
                static_ctr = 0
                soln_best = deepcopy(soln_current)
                obj_best = obj_current
                total_dist_best = total_dist_current
                info_log && @info "$iteration\t| $perturbation_type: $perturbed_clusters" *
                                  "\t| $(round(obj_best, digits=4))\t| " *
                                  "$(round(obj_current, digits=4))\t| " *
                                  "$(round(obj_proposed, digits=4))\t| $temp"
                plot_flag && display(Plot.solution(
                    problem, soln_best;
                    title="New Best Solution at Iteration $iteration",
                    highlight_critical_path_flag=true,
                ))
            end
        end

        # Record after accept/reject so current/best reflect the updated state
        push!(trace_iters, iteration)
        push!(trace_proposed, obj_proposed)
        push!(trace_current, obj_current)
        push!(trace_best, obj_best)
        push!(trace_temps, temp)

        (iteration >= min_iters && static_ctr >= static_limit) && break
        temp *= cooling_rate
        static_ctr += 1
    end

    # Display trace after each cluster's SA run completes
    fig_trace = Plot.sa_trace(
        trace_iters, trace_best, trace_current, trace_proposed, trace_temps
    )
    plot_flag && display(fig_trace)
    !isempty(output_dir) && CairoMakie.save("$output_dir/2_sa_trace.png", fig_trace)
    info_log && @info ("Final Value:\t$obj_best\nΔ: $(obj_init - obj_best)\n" *
                       "Iterations: $(trace_iters[end])")
    return soln_best, obj_best
end
