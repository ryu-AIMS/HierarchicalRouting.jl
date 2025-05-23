
"""
    perturb_swap_solution(
        soln::MSTSolution,
        clust_seq_idx_target::Int64=-1,
        exclusions_tender::DataFrame = DataFrame();
        enforce_diff_sortie::Bool = false
    )::MSTSolution
    perturb_swap_solution(
        soln::MSTSolution,
        cluster_pair::Tuple{Int, Int},
        exclusions_mothership::DataFrame = DataFrame(),
        exclusions_tender::DataFrame = DataFrame()
    )::MSTSolution

Perturb the solution by swapping two nodes:
- within a cluster, or
- between two clusters.

# Arguments
- `soln`: Solution to perturb.
- `clust_seq_idx_target`: Sequence index of cluster to perturb. Default = -1: randomly
    selects cluster.
- `cluster_pair`: Tuple of two cluster sequence indices to swap nodes between.
- `exclusions_tender`: DataFrame of exclusions for tender. Default = DataFrame().
- `exclusions_mothership`: DataFrame of exclusions for mothership. Default = DataFrame().
- `enforce_diff_sortie`: Boolean flag to enforce different sorties for node swaps.

# Returns
Perturbed full solution.
"""
function perturb_swap_solution(
    soln::MSTSolution,
    clust_seq_idx_target::Int64=-1,
    exclusions_tender::DataFrame = DataFrame();
    enforce_diff_sortie::Bool = false
)::MSTSolution
    clust_seq_idx = clust_seq_idx_target == -1 ?
        rand(1:length(soln.tenders[end])) :
        clust_seq_idx_target

    tender = deepcopy(soln.tenders[end][clust_seq_idx])
    sorties = deepcopy(tender.sorties)

    # If < 2 sorties in cluster, no perturbation possible
    if length(sorties) < 2
        return soln
    end

    # Choose two random sorties to swap nodes between
    # if enforce_diff_sortie is true: ensure different sorties are selected
    sortie_a_idx, sortie_b_idx = enforce_diff_sortie ?
        shuffle(1:length(sorties))[1:2] :
        rand(1:length(sorties), 2)

    sortie_a, sortie_b = sorties[sortie_a_idx], sorties[sortie_b_idx]

    # No perturbation possible if a sortie has no nodes
    if isempty(sortie_a.nodes) || isempty(sortie_b.nodes)
        return soln
    end

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
    updated_tender_tours::Vector{Vector{Point{2, Float64}}} = [
        [[tender.start]; sortie_a.nodes; [tender.finish]],
        [[tender.start]; sortie_b.nodes; [tender.finish]]
    ]

    # Get new linestrings
    linestrings::Vector{Vector{Vector{LineString{2, Float64}}}} = getindex.(
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

    tenders_all::Vector{TenderSolution} = copy(soln.tenders[end])
    tenders_all[clust_seq_idx] = TenderSolution(
        tender.id,
        tender.start,
        tender.finish,
        sorties,
        tender.dist_matrix  #? recompute
    )
    tenders_full_updated::Vector{Vector{TenderSolution}} = [
        copy(soln.tenders[end]),
        tenders_all
    ]
    return MSTSolution(soln.cluster_sets, soln.mothership_routes, tenders_full_updated)
end
function perturb_swap_solution(
    soln::MSTSolution,
    cluster_pair::Tuple{Int, Int},
    exclusions_mothership::DataFrame = DataFrame(),
    exclusions_tender::DataFrame = DataFrame()
)::MSTSolution
    clust_a_seq_idx, clust_b_seq_idx = cluster_pair

    tender_a = soln.tenders[end][clust_a_seq_idx]
    tender_b = soln.tenders[end][clust_b_seq_idx]

    # Pick random sorties and ensure both have nodes
    sortie_a_idx = rand(1:length(tender_a.sorties))
    sortie_b_idx = rand(1:length(tender_b.sorties))
    sortie_a = deepcopy(tender_a.sorties[sortie_a_idx])
    sortie_b = deepcopy(tender_b.sorties[sortie_b_idx])

    if isempty(sortie_a.nodes) || isempty(sortie_b.nodes)
        return soln
    end

    # Swap two random nodes across the two sorties
    node_a_idx, node_b_idx = rand(1:length(sortie_a.nodes)), rand(1:length(sortie_b.nodes))
    node_a, node_b = sortie_a.nodes[node_a_idx], sortie_b.nodes[node_b_idx]
    sortie_a.nodes[node_a_idx] = node_b
    sortie_b.nodes[node_b_idx] = node_a

    # Update new clusters
    new_clusters::Vector{Cluster} = deepcopy(soln.cluster_sets[end])

    cluster_a_idx = tender_a.id
    cluster_b_idx = tender_b.id
    nodes_a = new_clusters[cluster_a_idx].nodes
    nodes_b = new_clusters[cluster_b_idx].nodes

    # Swap nodes between clusters
    node_a_idx_clust = findfirst(isequal(node_a), nodes_a)
    node_b_idx_clust = findfirst(isequal(node_b), nodes_b)
    nodes_a[node_a_idx_clust] = node_b
    nodes_b[node_b_idx_clust] = node_a

    centroid_a, centroid_b = Point{2, Float64}.([
        (mean(getindex.(nodes_a, 1)), mean(getindex.(nodes_a, 2))),
        (mean(getindex.(nodes_b, 1)), mean(getindex.(nodes_b, 2)))
    ])

    new_clusters[cluster_a_idx], new_clusters[cluster_b_idx] = Cluster.(
        [cluster_a_idx, cluster_b_idx],
        [centroid_a, centroid_b],
        [nodes_a, nodes_b],
    )

    # Update mothership route and waypoints based on updated clusters
    cluster_seq_ids = getfield.(soln.tenders[end], :id)
    cluster_centroids = getfield.(new_clusters, :centroid)
    ordered_cluster_centroids = cluster_centroids[cluster_seq_ids]
    depot = soln.mothership_routes[end].route.nodes[1]
    full_ms_route_pts = [[depot]; ordered_cluster_centroids; [depot]]

    cluster_sequence = DataFrame(
        id = [0; cluster_seq_ids; 0],
        lon = getindex.(full_ms_route_pts, 1),
        lat = getindex.(full_ms_route_pts, 2)
    )

    updated_waypoints = get_waypoints(cluster_sequence, exclusions_mothership)
    waypoint_path_vector = get_feasible_vector(
        updated_waypoints.waypoint,
        exclusions_mothership
    )[2]

    ordered_clusters = sort(cluster_sequence[1:end-1, :], :id)
    ordered_centroid_pts = Point{2, Float64}.(ordered_clusters.lon, ordered_clusters.lat)
    full_ms_route_matrix = get_feasible_matrix(
        ordered_centroid_pts,
        exclusions_mothership
    )[1]

    updated_ms_solution = MothershipSolution(
        cluster_sequence,
        Route(
            updated_waypoints.waypoint,
            full_ms_route_matrix,
            vcat(waypoint_path_vector...)
        )
    )

    # Update routes for modified sorties
    tours = [
        [[tender_a.start]; sortie_a.nodes; [tender_a.finish]],
        [[tender_b.start]; sortie_b.nodes; [tender_b.finish]]
    ]
    updated_linestrings = getindex.(
        get_feasible_vector.(tours, Ref(exclusions_tender)),
        2
    )
    updated_sortie_matrices = getindex.(
        get_feasible_matrix.(tours, Ref(exclusions_tender)),
        1
    )

    sortie_a = Route(
        sortie_a.nodes,
        updated_sortie_matrices[1],
        vcat(updated_linestrings[1]...)
    )
    sortie_b = Route(
        sortie_b.nodes,
        updated_sortie_matrices[2],
        vcat(updated_linestrings[2]...)
    )

    # Re-compute sorties for the modified clusters
    sorties_a = [
        i == sortie_a_idx ? sortie_a : tender_a.sorties[i]
        for i in 1:length(tender_a.sorties)
    ]
    sorties_b = [
        i == sortie_b_idx ? sortie_b : tender_b.sorties[i]
        for i in 1:length(tender_b.sorties)
    ]

    # Update tenders with existing start/finish (not yet adjusted)
    tenders_all::Vector{TenderSolution} = copy(soln.tenders[end])
    tenders_all[clust_a_seq_idx] = TenderSolution(
        tender_a.id,
        tender_a.start,
        tender_a.finish,
        sorties_a,
        tender_a.dist_matrix
    )
    tenders_all[clust_b_seq_idx] = TenderSolution(
        tender_b.id,
        tender_b.start,
        tender_b.finish,
        sorties_b,
        tender_b.dist_matrix
    )
    tenders_full_updated::Vector{Vector{TenderSolution}} = [
        copy(soln.tenders[end]),
        tenders_all
    ]
    # Create new perturbed solution
    soln_perturbed = MSTSolution(
        [new_clusters],
        [updated_ms_solution],
        tenders_full_updated
    )
    return soln_perturbed
end

"""
    find_unallocated_nodes(
        soln::MSTSolution
    )::Set{Point{2, Float64}}
    find_unallocated_nodes(
        clusters::Vector{Cluster},
        tenders::Vector{TenderSolution}
    )::Set{Point{2, Float64}}

Find nodes that are not allocated to any sortie in the solution.

# Arguments
- `soln`: Solution to find unallocated nodes in.
- `clusters`: Vector of clusters.
- `tenders`: Vector of tender solutions.

# Returns
 A set of unallocated nodes.
"""
function find_unallocated_nodes(
    soln::MSTSolution
)::Set{Point{2, Float64}}
    clusters = soln.cluster_sets[end]
    tenders = soln.tenders[end]
    return find_unallocated_nodes(clusters, tenders)
end
function find_unallocated_nodes(
    clusters::Vector{Cluster},
    tenders::Vector{TenderSolution}
)::Set{Point{2, Float64}}
    cluster_sorties = getfield.(tenders, :sorties)

    all_nodes::Set{Point{2, Float64}} = Set(vcat(getfield.(clusters, :nodes)...))
    allocated_nodes::Set{Point{2, Float64}} = Set(
        collect(Base.Iterators.flatten(
            getfield.(Base.Iterators.flatten(cluster_sorties), :nodes))
        )
    )
    return setdiff(all_nodes, allocated_nodes)
end

"""
    insert_unallocated_node(
        soln::MSTSolution,
        exclusions::DataFrame = DataFrame()
        max_dist::Int64 = 15000;
        t_cap::Int16,
        current_cluster_seq_idx::Int=-1,
    )::MSTSolution
    insert_unallocated_node(
        clusters_ex::Vector{Cluster},
        tenders_ex::Vector{TenderSolution},
        exclusions::DataFrame = DataFrame(),
        max_dist::Float64 = 15000;
        t_cap::Int16,
        current_cluster_seq_idx::Int=-1,
    )::Vector{Cluster}

Insert unallocated nodes into the solution.

# Arguments
- `soln`: Solution to insert unallocated nodes into.
- `clusters_ex`: Vector of existing clusters.
- `tenders_ex`: Vector of existing tender solutions.
- `exclusions`: DataFrame of exclusion zones to avoid. Default = DataFrame().
- `max_dist`: Maximum distance to consider for cluster assignment - feasible distance (m)
    from start and finish waypoints to the unallocated node. Default = 15000.
- `t_cap`: Tender capacity.
- `current_cluster_seq_idx`: Sequence index of the current cluster.
    Default = -1: randomly selects cluster.
"""
function insert_unallocated_node(
    soln::MSTSolution,
    exclusions::DataFrame = DataFrame(),
    max_dist::Int64 = 15000;
    t_cap::Int16,
    current_cluster_seq_idx::Int=-1,
)::MSTSolution
    clusters = deepcopy(soln.cluster_sets[end])
    tenders = deepcopy(soln.tenders[end])
    unallocated_nodes = find_unallocated_nodes(soln)

    waypoints = soln.mothership_routes[end].route.nodes[2:end-1]
    min_tender_sorties = [minimum(length.(getfield.(t.sorties, :nodes))) for t in tenders]
    updated_tenders = copy(tenders)

    for node in unallocated_nodes
        centroids = getfield.(clusters, :centroid)

        # Find the closest cluster by shortest distance to ms path
        waypoint_distances = getindex.(shortest_feasible_path.(
            Ref(node),
            waypoints,
            Ref(exclusions)
        ), 1)
        @views cluster_waypoint_distances =
            waypoint_distances[1:2:end] .+
            waypoint_distances[2:2:end]

        cluster_mask = (min_tender_sorties .< t_cap) .&
            (cluster_waypoint_distances .< max_dist)

        if isempty(cluster_mask)
            break
        end

        # Find the closest cluster with available tender capacity to the unallocated node
        valid_idxs = findall(cluster_mask)
        cluster_seq_idx = valid_idxs[argmin(cluster_waypoint_distances[cluster_mask])]
        cluster_idx = tenders[cluster_seq_idx].id

        # Update cluster: Insert unallocated node
        cluster_nodes = vcat(clusters[cluster_idx].nodes, [node])
        cluster_centroid = Point{2, Float64}(
            mean(getindex.(cluster_nodes, 1)),
            mean(getindex.(cluster_nodes, 2))
        )
        clusters[cluster_idx] = Cluster(
            id = clusters[cluster_idx].id,
            centroid = cluster_centroid,
            nodes = cluster_nodes
        )

        # Update tender sortie: Add unallocated node to an available tender sortie
        sorties = updated_tenders[cluster_seq_idx].sorties
        sortie_mask = length.(getfield.(sorties, :nodes)) .< t_cap
        if !any(sortie_mask)
            break
        end

        # Find the closest available sortie to the unallocated node
        allocated_sortie_nodes = collect(Base.Iterators.flatten(
            getfield.(sorties[sortie_mask], :nodes)
        ))
        sortie_distances = getindex.(shortest_feasible_path.(
            Ref(node),
            allocated_sortie_nodes,
            Ref(exclusions)
        ), 1)
        closest_available_sortie_node = allocated_sortie_nodes[argmin(sortie_distances)]

        closest_available_sortie_idx = findfirst(
            x -> closest_available_sortie_node in x,
            getfield.(sorties, :nodes)
        )

        # Create updated sortie::Route
        updated_sortie_nodes = vcat(sorties[closest_available_sortie_idx].nodes, [node])
        updated_full_tour = [
            [tenders[cluster_seq_idx].start];
            updated_sortie_nodes;
            [tenders[cluster_seq_idx].finish]
        ]
        updated_sortie_matrix = get_feasible_matrix(updated_full_tour, exclusions)[1]
        #? Check these linestrings from get_feasible_vector() vs get_feasible_matrix()[2]
        updated_sortie_linestrings = get_feasible_vector(updated_full_tour, exclusions)[2]
        updated_sortie_linestrings_flat = vcat(updated_sortie_linestrings...)

        updated_sortie = Route(
            updated_sortie_nodes,
            updated_sortie_matrix,
            updated_sortie_linestrings_flat
        )

        # Update the sortie in the tender solution
        updated_sorties = [
            i == closest_available_sortie_idx ? updated_sortie : sorties[i]
            for i in 1:length(sorties)
        ]

        #? Does this modify the original solution in-place?
        updated_tenders[cluster_seq_idx] = TenderSolution(
            tenders[cluster_seq_idx].id,
            tenders[cluster_seq_idx].start,
            tenders[cluster_seq_idx].finish,
            updated_sorties,
            tenders[cluster_seq_idx].dist_matrix #! recompute or un-used??
        )

        #! update min_tender_sorties
        min_tender_sorties[cluster_seq_idx] = minimum(
            length.(getfield.(sorties, :nodes))
        )
    end
    tenders_full_updated::Vector{Vector{TenderSolution}} = [
        copy(soln.tenders[end]),
        updated_tenders
    ]
    soln_new = MSTSolution(
        [copy(soln.cluster_sets[end]), clusters],
        [
            deepcopy(soln.mothership_routes[end]),
            deepcopy(soln.mothership_routes[end])  # TODO: recompute if needed
        ],
        tenders_full_updated
    )
    return soln_new
end
function insert_unallocated_node(
    clusters_ex::Vector{Cluster},
    tenders_ex::Vector{TenderSolution},
    exclusions::DataFrame = DataFrame(),
    max_dist::Int64 = 15000;
    t_cap::Int16,
    current_cluster_seq_idx::Int=-1,
)::Vector{Cluster}
    clusters = deepcopy(clusters_ex)
    tenders = deepcopy(tenders_ex)
    unallocated_nodes = find_unallocated_nodes(clusters_ex, tenders)

    min_tender_sorties = [
        minimum(length.(getfield.(t.sorties, :nodes)))
        for t in tenders
    ]

    for node in unallocated_nodes
        waypoint_distances = [
            sum(getindex.(shortest_feasible_path.(
                Ref(node),
                [t.start, t.finish],
                Ref(exclusions)
            ),
            1))
            for t in tenders
        ]

        cluster_mask = (min_tender_sorties .< t_cap) .& (waypoint_distances .< max_dist)

        if isempty(cluster_mask)
            continue
        end
        # Select cluster closest by centroid distance (not waypoint this time)
        valid_idxs = getfield.(tenders[findall(cluster_mask)], :id)
        cluster_idx = valid_idxs[argmin(waypoint_distances[cluster_mask])]
        cluster_nodes = vcat(clusters[cluster_idx].nodes, [node])
        cluster_centroid = Point{2, Float64}(
            mean(getindex.(cluster_nodes, 1)),
            mean(getindex.(cluster_nodes, 2))
        )

        clusters[cluster_idx] = Cluster(
            id = clusters[cluster_idx].id,
            centroid = cluster_centroid,
            nodes = cluster_nodes
        )

        # Update tender sortie: Add unallocated node to an available tender sortie
        cluster_seq_idx = findfirst(==(cluster_idx), getfield.(tenders, :id))
        sorties = tenders[cluster_seq_idx].sorties
        sortie_mask = length.(getfield.(sorties, :nodes)) .< t_cap

        # Find the closest available sortie to the unallocated node
        allocated_sortie_nodes = collect(Base.Iterators.flatten(
            getfield.(sorties[sortie_mask], :nodes)
        ))
        sortie_distances = getindex.(shortest_feasible_path.(
            Ref(node),
            allocated_sortie_nodes,
            Ref(exclusions)
        ), 1)
        closest_available_sortie_node = allocated_sortie_nodes[argmin(sortie_distances)]

        closest_available_sortie_idx = findfirst(
            x -> closest_available_sortie_node in x,
            getfield.(sorties, :nodes)
        )

        updated_full_tour = [
            [tenders[cluster_seq_idx].start];
            vcat(sorties[closest_available_sortie_idx].nodes, [node]);
            [tenders[cluster_seq_idx].finish]
        ]

        updated_sortie_matrix = get_feasible_matrix(updated_full_tour, exclusions)[1]
        updated_sortie_linestrings = get_feasible_vector(updated_full_tour, exclusions)[2]
        updated_sortie_linestrings_flat = vcat(updated_sortie_linestrings...)
        updated_sortie = Route(
            updated_full_tour[2:end-1],
            updated_sortie_matrix,
            updated_sortie_linestrings_flat
        )

        # Update the sortie in the tender solution
        tenders[cluster_seq_idx].sorties[closest_available_sortie_idx] = updated_sortie

        # Update capacity tracker
        min_tender_sorties[cluster_seq_idx] = minimum(
            length.(getfield.(tenders[cluster_seq_idx].sorties, :nodes))
        )
    end

    return clusters
end

"""
    simulated_annealing(
        soln_init::MSTSolution,
        objective_function::Function,
        perturb_function::Function,
        exclusions_mothership::DataFrame = DataFrame(),
        exclusions_tender::DataFrame = DataFrame(),
        max_iterations::Int = 5_000,
        temp_init::Float64 = 500.0,
        cooling_rate::Float64 = 0.95,
        static_limit::Int = 150;
        cross_cluster_flag::Bool = false,
    )

Simulated Annealing optimization algorithm to optimize the solution.

# Arguments
- `soln_init`: Initial solution.
- `objective_function`: Function to evaluate the solution.
- `perturb_function`: Function to perturb the solution.
- `exclusions_mothership`: DataFrame of exclusions for mothership. Default = DataFrame().
- `exclusions_tender`: DataFrame of exclusions for tender. Default = DataFrame().
- `max_iterations`: Maximum number of iterations. Default = 5_000.
- `temp_init`: Initial temperature. Default = 500.0.
- `cooling_rate`: Rate of cooling to guide acceptance probability for SA algorithm.
    Default = 0.95 = 95%.
- `static_limit`: Number of iterations to allow stagnation before early exit. Default = 150.
- `cross_cluster_flag`: Boolean flag to indicate if perturbation across clusters should be
    considered. Default = false.

# Returns
- `soln_best`: Best solution::MSTSolution found.
- `obj_best`: Objective value of the best solution.
"""
function simulated_annealing(
    soln_init::MSTSolution,
    objective_function::Function,
    perturb_function::Function,
    exclusions_mothership::DataFrame = DataFrame(),
    exclusions_tender::DataFrame = DataFrame(),
    max_iterations::Int = 5_000,
    temp_init::Float64 = 500.0,
    cooling_rate::Float64 = 0.95,
    static_limit::Int = 150;
    cross_cluster_flag::Bool = false,
)::Tuple{MSTSolution, Float64}
    # Initialize best solution as initial
    soln_best = deepcopy(soln_init)
    obj_best = objective_function(soln_init)
    cluster_set::Vector{Cluster} = soln_init.cluster_sets[end]
    # TODO: Display cost values for each cluster, rather than total solution cost
    for clust_idx in 1:length(cluster_set)
        # Initialize current solution as best, reset temp
        temp = temp_init
        soln_current = deepcopy(soln_best)
        obj_current = obj_best
        static_ctr = 0

        @info "\n\tCluster\t\t$(cluster_set[clust_idx].id)\n\t" *
            "Iteration \tBest Value \t\tTemp\n\t0\t\t$obj_best\t$temp"

        for iteration in 1:max_iterations
            if !cross_cluster_flag
                soln_proposed = perturb_function(soln_current, clust_idx, exclusions_tender)
            else
                if rand() < 0.5
                # swap two nodes within the same cluster
                soln_proposed = perturb_function(soln_current, clust_idx, exclusions_tender)
                else
                    # swap two nodes between two different random clusters
                    clust_swap_idx = shuffle(setdiff(1:length(cluster_set), clust_idx))[1]
                    soln_proposed = perturb_function(
                        soln_current,
                        (clust_idx, clust_swap_idx),
                        exclusions_mothership,
                        exclusions_tender
                    )
                end
            end
            obj_proposed = objective_function(soln_proposed)
            improvement = obj_current - obj_proposed
            static_ctr += 1

            # If the new solution is improved OR meets acceptance prob criteria
            if improvement > 0 || exp(improvement / temp) > rand()
                soln_current = soln_proposed
                obj_current = obj_proposed

                if obj_current < obj_best
                    static_ctr = 0
                    soln_best = deepcopy(soln_current)
                    obj_best = obj_current
                    @info "$iteration\t\t$obj_best\t$temp"
                end
            end

            temp *= cooling_rate

            if iteration % 100 == 0
                @info "$iteration\t\t$obj_best\t$temp"
            end

            if static_ctr >= static_limit
                @info "$iteration\t\t$obj_best\t$temp\n\t" *
                    "Early exit at iteration $iteration due to stagnation."
                break
            end
        end
    end

    @info "\nFinal Value: $obj_best\nΔ: $(objective_function(soln_init) - obj_best)"
    return soln_best, obj_best
end
