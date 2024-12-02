
function initial_solution(problem::Problem)

    # Load problem data
    clusters, cluster_centroids_df = process_problem(problem)

    # Nearest Neighbour to generate initial mothership route & matrix
    ms_soln_NN, ms_feasible_matrix_centroids, ms_feasible_path_centroids = nearest_neighbour(cluster_centroids_df, problem.mothership.exclusion)

    # 2-opt to improve the NN soln
    ms_soln_2opt = two_opt(ms_soln_NN, ms_feasible_matrix_centroids, ms_feasible_path_centroids, problem.mothership.exclusion)

    clust_seq = [i for i in ms_soln_2opt.cluster_sequence.id if i!==0 && i <= length(clusters)]
    tender_soln = HierarchicalRouting.TenderSolution[]

    for (i, cluster_id) in enumerate(clust_seq)
        start_waypoint =  ms_soln_2opt.route.waypoint[2 * i]
        end_waypoint =  ms_soln_2opt.route.waypoint[2 * i + 1]
        @info "$(i): Clust $(cluster_id) from $(start_waypoint) to $(end_waypoint)"

        t_solution = tender_sequential_nearest_neighbour(
            clusters[cluster_id],
            (start_waypoint, end_waypoint),
            problem.tenders.number, problem.tenders.capacity, problem.tenders.exclusion)

        push!(tender_soln, t_solution[1])
    end

    return MSTSolution(clusters, ms_soln_2opt, tender_soln)#, feasible_path, [Point{2, Float64}(row.lat, row.lon) for row in eachrow(cluster_centroids_df)]
end

function improve_solution(soln::MSTSolution, opt_function::Function, objective_function::Function, perturb_function::Function, max_iterations::Int = 100_000, temp_init::Float64 = 500.0, cooling_rate::Float64 = 0.99_99)
    soln_best, z_best = opt_function(soln, objective_function, perturb_function, max_iterations, temp_init, cooling_rate)
    return soln_best, z_best
end
