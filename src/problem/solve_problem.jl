
"""
    initial_solution(
        problem::Problem;
        disturbance_clusters::Set{Int64}=Set{Int64}()
    )::MSTSolution

Generate a solution to the problem for:
- the mothership by using the nearest neighbour heuristic to generate an initial solution,
    and then improving using the 2-opt heuristic, and
-  for tenders using the sequential nearest neighbour heuristic.

# Arguments
- `problem`: Problem instance to solve
- `disturbance_clusters`: Set of sequenced clusters to simulate disturbances before.

# Returns
Best total MSTSolution found
"""
function initial_solution(
    problem::Problem;
    k::Int=1,
    disturbance_clusters::Set{Int64}=Set{Int64}()
)::MSTSolution
    n_events = length(disturbance_clusters) + 1
    cluster_sets = Vector{Vector{Cluster}}(undef, n_events)
    ms_soln_sets = Vector{MothershipSolution}(undef, n_events)
    tender_soln_sets = Vector{Vector{TenderSolution}}(undef, n_events)
    total_tender_capacity = Int(problem.tenders.number * problem.tenders.capacity)

    # Load problem data
    clusters::Vector{Cluster} = cluster_problem(
        problem; k
    )
    cluster_centroids_df::DataFrame = generate_cluster_df(clusters, problem.depot)

    ms_route::MothershipSolution = optimize_mothership_route(problem, cluster_centroids_df)
    clust_seq::Vector{Int64} = filter(
        i -> i != 0 && i <= length(clusters),
        ms_route.cluster_sequence.id
    )
    initial_tenders = [
        tender_sequential_nearest_neighbour(
            clusters[clust_seq][j],
            (ms_route.route.nodes[2j], ms_route.route.nodes[2j+1]),
            problem.tenders.number,
            problem.tenders.capacity,
            problem.tenders.exclusion
        )
        for j in 1:length(clust_seq)
    ]

    cluster_sets[1] = deepcopy(clusters)
    ms_soln_sets[1] = ms_route
    tender_soln_sets[1] = initial_tenders

    disturbance_index = 1
    for i ∈ sort(collect(disturbance_clusters))
        @info "Disturbance event #$disturbance_index at $(ms_route.route.nodes[2i-1]) " *
              "before $(i)th cluster_id=$(clust_seq[i])"
        disturbance_index += 1
        clusters = vcat(
            clusters[clust_seq][1:i-1],
            disturb_clusters(
                clusters[clust_seq][i:end],
                problem.targets.disturbance_gdf,
                ms_route.route.nodes[2i-1],
                problem.tenders.exclusion,
                total_tender_capacity
            )
        )

        removed_nodes = setdiff(
            vcat([c.nodes for c in cluster_sets[disturbance_index-1]]...),
            vcat([c.nodes for c in clusters]...)
        )
        if !isempty(removed_nodes)
            @info "Removed nodes due to disturbance event (since previous cluster):\n" *
                  "\t$(join(removed_nodes, "\n\t"))"
        end

        sort!(clusters, by=x -> x.id)

        cluster_centroids_df = generate_cluster_df(clusters, problem.depot)

        ms_route = optimize_mothership_route(
            problem,
            cluster_centroids_df,
            i,
            ms_route,
            getfield.(clusters[clust_seq][1:i-1], :id)
        )
        clust_seq = filter(
            i -> i != 0 && i <= length(clusters),
            ms_route.cluster_sequence.id
        )

        cluster_sets[disturbance_index] = clusters
        ms_soln_sets[disturbance_index] = ms_route
        current_tender_soln = Vector{TenderSolution}(undef, length(clust_seq))
        for j in 1:length(clust_seq)
            if j < i
                current_tender_soln[j] = tender_soln_sets[disturbance_index-1][j]
            else
                current_tender_soln[j] = tender_sequential_nearest_neighbour(
                    clusters[clust_seq][j],
                    (ms_route.route.nodes[2j], ms_route.route.nodes[2j+1]),
                    problem.tenders.number,
                    problem.tenders.capacity,
                    problem.tenders.exclusion
                )
            end
            tender_soln_sets[disturbance_index] = current_tender_soln
        end
    end

    return MSTSolution(cluster_sets, ms_soln_sets, tender_soln_sets)
end

"""
    solve(
        problem::Problem;
        k::Int=1,
        disturbance_clusters::Set{Int64}=Set{Int64}()
        seed::Union{Nothing,Int64}=nothing,
        rng::AbstractRNG=Random.GLOBAL_RNG
    )::MSTSolution

Generate a solution to the problem for:
- the mothership by using the nearest neighbour heuristic to generate an initial solution,
    and then improving using the 2-opt heuristic, and
-  for tenders using the sequential nearest neighbour heuristic.

# Arguments
- `problem`: Problem instance to solve
- `k`: Number of clusters to generate
- `disturbance_clusters`: Set of sequenced clusters to simulate disturbances before.
- `seed`: Optional seed for random number generation
- `rng`: AbstractRNG for random number generation

# Returns
Best total MSTSolution found
"""
function solve(
    problem::Problem;
    k::Int=1,
    disturbance_clusters::Set{Int64}=Set{Int64}(),
    seed::Union{Nothing,Int64}=nothing,
    rng::AbstractRNG=Random.GLOBAL_RNG
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

    # Cluster the problem data
    clusters::Vector{Cluster} = cluster_problem(problem; k);
    cluster_centroids_df::DataFrame = generate_cluster_df(clusters, problem.depot)

    # Route the mothership using nearest neighbour and 2-opt
    ms_route::MothershipSolution = optimize_mothership_route(problem, cluster_centroids_df)
    clust_seq::Vector{Int64} = filter(
        i -> i != 0 && i <= length(clusters),
        ms_route.cluster_sequence.id
    )

    # Generate initial tender solutions using sequential nearest neighbour
    initial_tenders = [
        tender_sequential_nearest_neighbour(
            clusters[clust_seq][j],
            (ms_route.route.nodes[2j], ms_route.route.nodes[2j+1]),
            problem.tenders.number,
            problem.tenders.capacity,
            problem.tenders.exclusion
        )
        for j in 1:length(clust_seq)
    ]

    # Optimize the initial tenders solution up to the first disturbance
    next_cluster_idx = !isempty(ordered_disturbances) ?
                       ordered_disturbances[1] :
                       length(clusters)

    optimized_initial, _ = improve_solution(
        MSTSolution([clusters], [ms_route], [initial_tenders]),
        problem.mothership.exclusion,
        problem.tenders.exclusion,
        1, next_cluster_idx
    )

    # Apply the optimized initial solution to the first set of clusters pre-disturbance
    cluster_sets[1], ms_soln_sets[1], tender_soln_sets[1] = apply_improved!(
        optimized_initial,
        clusters,
        initial_tenders,
    )

    # Iterate through each disturbance event and update solution
    disturbance_index_count = 1
    for disturbance_cluster_idx ∈ ordered_disturbances
        @info "Disturbance event #$disturbance_cluster_idx at " *
            "$(ms_route.route.nodes[2*disturbance_cluster_idx-1]) before " *
            "$(disturbance_cluster_idx)th cluster_id=$(clust_seq[disturbance_cluster_idx])"
        disturbance_index_count += 1
        # Update clusters based on the impact of disturbance event on future points/clusters
        clusters = sort!(
            vcat(
                clusters[clust_seq][1:disturbance_cluster_idx-1],
                disturb_clusters(
                    clusters[clust_seq][disturbance_cluster_idx:end],
                    problem.targets.disturbance_gdf,
                    ms_route.route.nodes[2*disturbance_cluster_idx-1],
                    problem.tenders.exclusion,
                    total_tender_capacity
                )
            ), by=x -> x.id
        )

        removed_nodes = setdiff(
            vcat([c.nodes for c in cluster_sets[disturbance_cluster_idx-1]]...),
            vcat([c.nodes for c in clusters]...)
        )
        if !isempty(removed_nodes)
            @info "Removed nodes due to disturbance event (since previous cluster):\n" *
                  "\t$(join(removed_nodes, "\n\t"))"
        end
        # Re-generate the cluster centroids to route mothership
        cluster_centroids_df = generate_cluster_df(clusters, problem.depot)

        ms_route = optimize_mothership_route(
            problem,
            cluster_centroids_df,
            disturbance_cluster_idx,
            ms_route,
            getfield.(clusters[clust_seq][1:disturbance_cluster_idx-1], :id)
        )
        clust_seq = filter(
            i -> i != 0 && i <= length(clusters),
            ms_route.cluster_sequence.id
        )

        # Update cluster sets, mothership solution, and tender solutions
        cluster_sets[disturbance_cluster_idx] = clusters
        ms_soln_sets[disturbance_cluster_idx] = ms_route
        current_tender_soln = Vector{TenderSolution}(undef, length(clust_seq))

        # Generate tender solutions for the current disturbance cluster
        for j in 1:length(clust_seq)
            if j < disturbance_cluster_idx
                current_tender_soln[j] = tender_soln_sets[disturbance_cluster_idx-1][j]
            else
                current_tender_soln[j] = tender_sequential_nearest_neighbour(
                    clusters[clust_seq][j],
                    (ms_route.route.nodes[2j], ms_route.route.nodes[2j+1]),
                    problem.tenders.number,
                    problem.tenders.capacity,
                    problem.tenders.exclusion
                )
            end
        end

        next_disturbance_cluster_idx = disturbance_cluster_idx <= length(ordered_disturbances) ?
            ordered_disturbances[disturbance_cluster_idx] :
            length(clusters) + 1

        # Improve the current solution using the optimization function
        optimized_current_solution, _ = improve_solution(
            MSTSolution([clusters], [ms_route], [current_tender_soln]),
            problem.mothership.exclusion,
            problem.tenders.exclusion,
            disturbance_cluster_idx, next_disturbance_cluster_idx
        )

        # Update with improved solution to the current cluster set and tender solutions
        cluster_sets[disturbance_cluster_idx],
        ms_soln_sets[disturbance_cluster_idx],
        tender_soln_sets[disturbance_cluster_idx] = apply_improved!(
            optimized_current_solution,
            clusters,
            current_tender_soln
        )
    end

    return MSTSolution(cluster_sets, ms_soln_sets, tender_soln_sets)
end

"""
    apply_improved!(
        improved_soln::MSTSolution,
        existing_clusters::Vector{Cluster},
        existing_tenders::Vector{TenderSolution},
    )::Tuple{Vector{Cluster}, MothershipSolution, Vector{TenderSolution}}

Update the clusters and tenders with the improved solution.

# Arguments
- `improved_soln`: The improved solution containing updated clusters and tenders.
- `existing_clusters`: Vector of clusters to update.
- `existing_tenders`: Vector of tender solutions for each event.

# Returns
- A tuple containing the updated clusters, the mothership route, and the updated tenders for the specified event index.
"""
function apply_improved!(
    improved_soln::MSTSolution,
    existing_clusters::Vector{Cluster},
    existing_tenders::Vector{TenderSolution},
)
    # Overwrite updated clusters
    updated_clusters = improved_soln.cluster_sets[1]
    ids = getfield.(updated_clusters, :id)
    existing_clusters[ids] .= updated_clusters

    #Overwrite updated tenders
    updated_tenders = improved_soln.tenders[1]
    tender_ids = getfield.(updated_tenders, :id)

    tender_idxs = findfirst.(.==(tender_ids), Ref(getfield.(existing_tenders, :id)))
    existing_tenders[tender_idxs] .= updated_tenders

    # Record the updated clusters and tenders for the current index
    return (
      deepcopy(existing_clusters),
      improved_soln.mothership_routes[1],
      existing_tenders
    )
end
    )
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
        cluster_centroids_df, problem.mothership.exclusion, problem.tenders.exclusion
    )

    # 2-opt to improve the NN soln
    ms_soln_2opt::MothershipSolution = two_opt(
        ms_soln_NN, problem.mothership.exclusion, problem.tenders.exclusion
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
        problem.mothership.exclusion,
        problem.tenders.exclusion,
        start_point,
        ms_route,
        cluster_seq_idx
    )

    # 2-opt to improve the NN soln
    ms_soln_2opt::MothershipSolution = two_opt(
        ms_soln_NN,
        problem.mothership.exclusion,
        problem.tenders.exclusion,
        cluster_seq_idx
    )

    return ms_soln_2opt
end

"""
    improve_solution(
        initial_solution::MSTSolution,
        exclusions_mothership::DataFrame=DataFrame(),
        exclusions_tender::DataFrame=DataFrame(),
        current_cluster_idx::Int=1,
        next_cluster_idx::Int=length(initial_solution.cluster_sets[end]);
        opt_function::Function=HierarchicalRouting.simulated_annealing,
        objective_function::Function=HierarchicalRouting.critical_path,
        perturb_function::Function=HierarchicalRouting.perturb_swap_solution,
        max_iterations::Int=1_000,
        temp_init::Float64=500.0,
        cooling_rate::Float64=0.95,
        static_limit::Int=20,
        vessel_weightings::NTuple{2,AbstractFloat}=(1.0, 1.0)
    )::Tuple{MSTSolution,Float64}

Improve the solution using the optimization function `opt_function` with the objective \n
function `objective_function` and the perturbation function `perturb_function`.

# Arguments
- `initial_solution`: Initial solution to improve
- `exclusions_mothership`: DataFrame of exclusion polygons for the mothership
- `exclusions_tender`: DataFrame of exclusion polygons for the tenders
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
    exclusions_mothership::DataFrame=DataFrame(),
    exclusions_tender::DataFrame=DataFrame(),
    current_cluster_idx::Int=1,
    next_cluster_idx::Int=length(initial_solution.cluster_sets[end]);
    opt_function::Function=HierarchicalRouting.simulated_annealing,
    objective_function::Function=HierarchicalRouting.critical_path,
    perturb_function::Function=HierarchicalRouting.perturb_swap_solution,
    max_iterations::Int=1_000,
    temp_init::Float64=500.0,
    cooling_rate::Float64=0.95,
    static_limit::Int=20,
    vessel_weightings::NTuple{2,AbstractFloat}=(1.0, 1.0)
)::Tuple{MSTSolution,Float64}
    current_mothership_route = initial_solution.mothership_routes[end]
    clust_seq_current = current_mothership_route.cluster_sequence.id[current_cluster_idx+1:next_cluster_idx]
    current_clusters = initial_solution.cluster_sets[end][clust_seq_current]
    current_tender_routes = initial_solution.tenders[end][current_cluster_idx:next_cluster_idx-1]
    current_solution = MSTSolution(
        [current_clusters],
        [current_mothership_route],
        [current_tender_routes]
    )

    soln_best, z_best = opt_function(
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

    return soln_best, z_best
end
