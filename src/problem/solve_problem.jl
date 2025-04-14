
"""
    initial_solution(problem::Problem)::MSTSolution

Generate a solution to the problem for:
- the mothership by using the nearest neighbour heuristic to generate an initial solution,
    and then improving using the 2-opt heuristic, and
-  for tenders using the sequential nearest neighbour heuristic.

# Arguments
- `problem`: Problem instance to solve

# Returns
Best total MSTSolution found
"""
function initial_solution(problem::Problem)::MSTSolution
    # Load problem data
    clusters::Vector{Cluster} = cluster_problem(problem);
    cluster_centroids_df::DataFrame = generate_cluster_df(clusters, problem.depot)

    ms_route::MothershipSolution = optimize_mothership_route(problem, cluster_centroids_df)
    clust_seq::Vector{Int64} = filter(
        i -> i != 0 && i <= length(clusters),
        ms_route.cluster_sequence.id
    )

    disturbance_clusters::Set{Int64} = Set((2, 4)) #! disturbance events are hardcoded for now at/before 2 and 4
    tender_soln = Vector{TenderSolution}(undef, length(clust_seq))
    cluster_set = Vector{Vector{Cluster}}(undef, length(disturbance_clusters)+1)

    ms_soln_sets = Vector{MothershipSolution}(undef, length(clust_seq))

    disturbance_index = 1
    cluster_set[1] = clusters
    ms_soln_sets[1] = ms_route

    for (i, cluster_id) in enumerate(clust_seq)
        if i ∈ disturbance_clusters
            @info "Disturbance event at $(ms_route.route.nodes[2*i-1]): before $(i)th cluster_id=$(cluster_id)"
            disturbance_index += 1
            cluster_set[disturbance_index] = clusters = sort(
                vcat(
                    clusters[clust_seq][1:i-1],
                    disturb_clusters(
                        clusters[clust_seq][i:end],
                        problem.targets.disturbance_gdf
                    )
                ),
                by = x -> x.id
            )
        end

        start_waypoint::Point{2, Float64} =  ms_route.route.nodes[2 * i]
        end_waypoint::Point{2, Float64} =  ms_route.route.nodes[2 * i + 1]
        @info "$(i): Clust $(cluster_id) from $(start_waypoint) to $(end_waypoint)"

        #? order by deployment sequence, rather than ID
        tender_soln[i] = tender_sequential_nearest_neighbour(
            clusters[cluster_id],
            (start_waypoint, end_waypoint),
            problem.tenders.number, problem.tenders.capacity, problem.tenders.exclusion
        )
    end

    return MSTSolution(cluster_set, ms_soln_sets, tender_soln)
end

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
