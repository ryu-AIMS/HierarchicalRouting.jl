
"""
    initial_solution(
        problem::Problem,
        num_clusters::Int;
        cluster_tolerance::Real = Float64(5E-5),
        disturbance_clusters::Set{Int64} = Set{Int64}()
    )::MSTSolution

Generate a solution to the problem for:
- the mothership by using the nearest neighbour heuristic to generate an initial solution,
    and then improving using the 2-opt heuristic, and
-  for tenders using the sequential nearest neighbour heuristic.

# Arguments
- `problem`: Problem instance to solve
- `num_clusters`: Number of clusters to create
- `cluster_tolerance`: Tolerance for clustering. Default is 5E-5.
- `disturbance_clusters`: Set of sequenced clusters to simulate disturbances before.

# Returns
Best total MSTSolution found
"""
function initial_solution(
    problem::Problem,
    num_clusters::Int;
    cluster_tolerance::Real = Float64(5E-5),
    disturbance_clusters::Set{Int64} = Set{Int64}()
)::MSTSolution
    num_clusters = Int8(num_clusters)
    cluster_tolerance = Float64(cluster_tolerance)

    n_events = length(disturbance_clusters) + 1
    cluster_sets = Vector{Vector{Cluster}}(undef, n_events)
    ms_soln_sets = Vector{MothershipSolution}(undef, n_events)
    tender_soln_sets = Vector{Vector{TenderSolution}}(undef, n_events)

    # Load problem data
    clusters::Vector{Cluster} = cluster_problem(
        problem,
    );
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
        @info "Disturbance event at $(ms_route.route.nodes[2i-1]) " *
            "before $(i)th cluster_id=$(clust_seq[i])"
        disturbance_index += 1
        clusters = vcat(
            clusters[clust_seq][1:i-1],
            disturb_clusters(
                clusters[clust_seq][i:end],
                problem.targets.disturbance_gdf,
                ms_route.route.nodes[2i-1],
                problem.tenders.exclusion
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

        sort!(clusters, by = x -> x.id)

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
        );

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
        );

        # 2-opt to improve the NN soln
        ms_soln_2opt::MothershipSolution = two_opt(
            ms_soln_NN, problem.mothership.exclusion, problem.tenders.exclusion
        );
    return ms_soln_2opt
end
function optimize_mothership_route(
    problem::Problem,
    cluster_centroids_df::DataFrame,
    cluster_seq_idx::Int64,
    ms_route::MothershipSolution,
    cluster_ids_visited::Vector{Int64}
)::MothershipSolution
    start_point::Point{2, Float64} =  ms_route.route.nodes[2 * cluster_seq_idx - 1]

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
    );

    # 2-opt to improve the NN soln
    ms_soln_2opt::MothershipSolution = two_opt(
        ms_soln_NN,
        problem.mothership.exclusion,
        problem.tenders.exclusion,
        cluster_seq_idx
    );

    return ms_soln_2opt
end

"""
    improve_solution(
        initial_solution::MSTSolution,
        exclusions_mothership::DataFrame = DataFrame(),
        exclusions_tender::DataFrame = DataFrame();
        opt_function::Function = HierarchicalRouting.simulated_annealing,
        objective_function::Function = HierarchicalRouting.critical_path,
        perturb_function::Function = HierarchicalRouting.perturb_swap_solution,
        max_iterations::Int = 1_000,
        temp_init::Float64 = 500.0,
        cooling_rate::Float64 = 0.95,
        static_limit::Int = 20
    )::Tuple{MSTSolution, Float64}

Improve the solution using the optimization function `opt_function` with the objective \n
function `objective_function` and the perturbation function `perturb_function`.

# Arguments
- `initial_solution`: Initial solution to improve
- `exclusions_mothership`: DataFrame of exclusion polygons for the mothership
- `exclusions_tender`: DataFrame of exclusion polygons for the tenders
- `opt_function`: Optimization function to improve the solution
- `objective_function`: Objective function to quantify and evaluate the solution
- `perturb_function`: Perturbation function to generate changes in the solution
- `max_iterations`: Maximum number of iterations
- `temp_init`: Initial temperature for simulated annealing
- `cooling_rate`: Cooling rate for simulated annealing
- `static_limit`: Number of iterations to allow stagnation before early exit

# Returns
- `soln_best`: Solution with lowest objective value found
- `z_best`: Objective value associated with the best solution
"""
function improve_solution(
    initial_solution::MSTSolution,
    exclusions_mothership::DataFrame = DataFrame(),
    exclusions_tender::DataFrame = DataFrame();
    opt_function::Function = HierarchicalRouting.simulated_annealing,
    objective_function::Function = HierarchicalRouting.critical_path,
    perturb_function::Function = HierarchicalRouting.perturb_swap_solution,
    max_iterations::Int = 1_000,
    temp_init::Float64 = 500.0,
    cooling_rate::Float64 = 0.95,
    static_limit::Int = 20
)::Tuple{MSTSolution, Float64}
    soln_best, z_best = opt_function(
        initial_solution,
        objective_function,
        perturb_function,
        exclusions_mothership,
        exclusions_tender,
        max_iterations,
        temp_init,
        cooling_rate,
        static_limit
    )
    return soln_best, z_best
end
