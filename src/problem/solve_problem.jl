
"""
    initial_solution(problem::Problem)

Generate an initial solution to the problem for mothership using the nearest neighbour
heuristic, and then improving the solution using the 2-opt heuristic; for tenders using
the sequential nearest neighbour heuristic.

# Arguments
- `problem::Problem`: Problem instance to solve

# Returns
- `soln_best::MSTSolution`: Best solution found
- `z_best::Float64`: Best objective value found
"""
function initial_solution(problem::Problem)
    # Load problem data
    clusters::Vector{Cluster} = process_problem(problem)

    cluster_centroids_df::DataFrame = DataFrame(
        id  = [0; 1:length(clusters)],
        lon = [problem.depot[1]; [clust.centroid[1] for clust in clusters]],
        lat = [problem.depot[2]; [clust.centroid[2] for clust in clusters]]
    )

    # Nearest Neighbour to generate initial mothership route & matrix
    ms_soln_NN = nearest_neighbour(
        cluster_centroids_df, problem.mothership.exclusion, problem.tenders.exclusion
    )

    # 2-opt to improve the NN soln
    ms_soln_2opt = two_opt(
        ms_soln_NN, problem.mothership.exclusion, problem.tenders.exclusion
    )

    clust_seq = filter(
        i -> i != 0 && i <= length(clusters),
        ms_soln_2opt.cluster_sequence.id
    )
    tender_soln = HierarchicalRouting.TenderSolution[]

    for (i, cluster_id) in enumerate(clust_seq)
        start_waypoint =  ms_soln_2opt.route.nodes[2 * i]
        end_waypoint =  ms_soln_2opt.route.nodes[2 * i + 1]
        @info "$(i): Clust $(cluster_id) from $(start_waypoint) to $(end_waypoint)"

        t_solution = tender_sequential_nearest_neighbour(
            clusters[cluster_id],
            (start_waypoint, end_waypoint),
            problem.tenders.number, problem.tenders.capacity, problem.tenders.exclusion)

        push!(tender_soln, t_solution)
    end

    return MSTSolution(clusters, ms_soln_2opt, tender_soln)
end

"""
    improve_solution(
        soln::MSTSolution,
        opt_function::Function,
        objective_function::Function,
        perturb_function::Function,
        exclusions::DataFrame = DataFrame(),
        max_iterations::Int = 100_000,
        temp_init::Float64 = 500.0,
        cooling_rate::Float64 = 0.95
    )

Improve the solution using the optimization function `opt_function`
    with the objective function `objective_function` and
    the perturbation function `perturb_function`.

# Arguments
- `soln::MSTSolution`: Initial solution to improve
- `opt_function::Function`: Optimization function to use
- `objective_function::Function`: Objective function to use
- `perturb_function::Function`: Perturbation function to use
- `exclusions::DataFrame = DataFrame()`: Exclusions DataFrame
- `max_iterations::Int = 5_000`: Maximum number of iterations
- `temp_init::Float64 = 500.0`: Initial temperature for simulated annealing
- `cooling_rate::Float64 = 0.95`: Cooling rate for simulated annealing

# Returns
- `soln_best::MSTSolution`: Best solution found
- `z_best::Float64`: Best objective value found
"""
function improve_solution(
    soln::MSTSolution,
    opt_function::Function,
    objective_function::Function,
    perturb_function::Function,
    exclusions::DataFrame = DataFrame(),
    max_iterations::Int = 5_000,
    temp_init::Float64 = 500.0,
    cooling_rate::Float64 = 0.95
)
    soln_best, z_best = opt_function(soln, objective_function, perturb_function, exclusions, max_iterations, temp_init, cooling_rate)
    return soln_best, z_best
end
