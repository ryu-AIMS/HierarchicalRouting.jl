
"""
    initial_solution(
        problem::Problem,
        disturbance_clusters::Set{Int64} = Set{Int64}()
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
    problem::Problem,
    num_clusters::Int8,
    cluster_tolerance::Float64,
    disturbance_clusters::Set{Int64} = Set{Int64}()
)::MSTSolution
    # Load problem data
    clusters::Vector{Cluster} = cluster_problem(
        problem,
        num_clusters,
        cluster_tolerance,
    );
    cluster_centroids_df::DataFrame = generate_cluster_df(clusters, problem.depot)

    ms_route::MothershipSolution = optimize_mothership_route(problem, cluster_centroids_df)
    clust_seq::Vector{Int64} = filter(
        i -> i != 0 && i <= length(clusters),
        ms_route.cluster_sequence.id
    )

    tender_soln = Vector{TenderSolution}(undef, length(clust_seq))
    cluster_set = Vector{Vector{Cluster}}(undef, length(disturbance_clusters)+1)
    ms_soln_sets = Vector{MothershipSolution}(undef, length(disturbance_clusters)+1)

    cluster_set[1] = clusters
    ms_soln_sets[1] = ms_route

    i = 1
    disturbance_index = 1
    while i <= length(clust_seq)
        cluster_id::Int64 = clust_seq[i];
        if i ∈ disturbance_clusters
            @info "Disturbance event at $(ms_route.route.nodes[2i-1]): before $(i)th cluster_id=$(cluster_id)"
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
                vcat([c.nodes for c in cluster_set[disturbance_index-1]]...),
                vcat([c.nodes for c in clusters]...)
            )
            if !isempty(removed_nodes)
                @info "Removed nodes due to disturbance event (since previous cluster):"
                @info join(removed_nodes, "\n\t")
            end

            sort!(clusters, by = x -> x.id)
            cluster_set[disturbance_index] = clusters

            cluster_centroids_df = generate_cluster_df(clusters, problem.depot)

            ms_soln_sets[disturbance_index] = ms_route = optimize_mothership_route(
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
        end

        start_waypoint::Point{2, Float64} = ms_route.route.nodes[2i]
        end_waypoint::Point{2, Float64} = ms_route.route.nodes[2i + 1]
        cluster::Cluster = clusters[clust_seq][i]
        @info "$(i): Clust $(cluster.id) from $(start_waypoint) to $(end_waypoint)"

        tender_soln[i] = tender_sequential_nearest_neighbour(
            cluster,
            (start_waypoint, end_waypoint),
            problem.tenders.number, problem.tenders.capacity, problem.tenders.exclusion
        )
        i+=1
    end

    return MSTSolution(cluster_set, ms_soln_sets, tender_soln)
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
        soln::MSTSolution,
        opt_function::Function,
        objective_function::Function,
        perturb_function::Function,
        exclusions::DataFrame = DataFrame();
        max_iterations::Int = 5_000,
        temp_init::Float64 = 500.0,
        cooling_rate::Float64 = 0.95,
        static_limit::Int = 150
    )::Tuple{MSTSolution, Float64}

Improve the solution using the optimization function `opt_function` with the objective \n
function `objective_function` and the perturbation function `perturb_function`.

# Arguments
- `soln`: Initial solution to improve
- `opt_function`: Optimization function to use
- `objective_function`: Objective function to use
- `perturb_function`: Perturbation function to use
- `exclusions`: Exclusions DataFrame
- `max_iterations`: Maximum number of iterations
- `temp_init`: Initial temperature for simulated annealing
- `cooling_rate`: Cooling rate for simulated annealing
- `static_limit`: Number of iterations to allow stagnation before early exit

# Returns
- `soln_best`: Solution with lowest objective value found
- `z_best`: Objective value associated with the best solution
"""
function improve_solution(
    soln::MSTSolution,
    opt_function::Function,
    objective_function::Function,
    perturb_function::Function,
    exclusions::DataFrame = DataFrame();
    max_iterations::Int = 5_000,
    temp_init::Float64 = 500.0,
    cooling_rate::Float64 = 0.95,
    static_limit::Int = 150
)::Tuple{MSTSolution, Float64}
    soln_best, z_best = opt_function(
        soln,
        objective_function,
        perturb_function,
        exclusions,
        max_iterations,
        temp_init,
        cooling_rate,
        static_limit
    )
    return soln_best, z_best
end
