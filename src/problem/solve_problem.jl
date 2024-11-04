
function initial_solution(problem::Problem)

    config = TOML.parsefile(joinpath("src",".config.toml"))

    suitable_threshold = config["parameters"]["suitable_threshold"]
    k = config["parameters"]["k"]
    cluster_tolerance = config["parameters"]["cluster_tolerance"]
    EPSG_code = config["parameters"]["EPSG_code"]

    target_scenario = problem.target_scenario
    suitable_targets_prefix = target_scenario[1:findlast(".",target_scenario)[1]-1]

    site_dir = config["data_dir"]["site"]
    output_dir = config["output_dir"]["path"]

    clustered_targets_path = joinpath(output_dir, "clustered_$(suitable_targets_prefix)_targets_k=$(k).tif")
    target_subset_threshold_path = joinpath(output_dir, "target_subset_$(suitable_targets_prefix)_threshold=$(suitable_threshold).tif")
    target_subset_path = joinpath(output_dir, "target_subset_$(suitable_targets_prefix).tif")

    subset = GDF.read(first(glob("*.gpkg", site_dir)))

    target_scenario_dir = config["data_dir"]["target_scenarios"]
    suitable_targets_all_path = joinpath(target_scenario_dir, target_scenario)

    clusters = cluster_targets(
        clustered_targets_path,
        target_subset_threshold_path,
        k, cluster_tolerance,
        suitable_targets_all_path,
        suitable_threshold,
        target_subset_path,
        subset,
        EPSG_code
    )

    cluster_centroids_df = DataFrame(id = Int[0], lon = Float64[problem.depot[1]], lat = Float64[problem.depot[2]])
    [push!(cluster_centroids_df, (i, clust.centroid[1], clust.centroid[2])) for (i, clust) in enumerate(clusters)]

    # Nearest Neighbour to generate initial mothership route & matrix
    ms_soln_NN, ms_feasible_matrix = nearest_neighbour(cluster_centroids_df, problem.mothership.exclusion)

    # 2-opt to improve the NN soln
    ms_soln_2opt = two_opt(ms_soln_NN.cluster_sequence, ms_feasible_matrix)

    clust_seq = [i for i in ms_soln_2opt.cluster_sequence.id if i!==0 && i <= length(clusters)]
    tender_soln_NN = HierarchicalRouting.TenderSolution[]

    for (i, cluster_id) in enumerate(clust_seq)
        start_waypoint =  ms_soln_2opt.route.waypoint[2 * i]
        end_waypoint =  ms_soln_2opt.route.waypoint[2 * i + 1]
        @info "$(i): Clust $(cluster_id) from $(start_waypoint) to $(end_waypoint)"

        t_solution = tender_sequential_nearest_neighbour(
            clusters[cluster_id],
            (start_waypoint, end_waypoint),
            problem.tenders.number, problem.tenders.capacity, problem.tenders.exclusion)

        push!(tender_soln_NN, t_solution[1])
    end

    return MSTSolution(clusters, ms_soln_2opt, tender_soln_NN)
end

function improve_solution(soln::MSTSolution, opt_function::Function, objective_function::Function, perturb_function::Function, max_iterations::Int = 100_000, temp_init::Float64 = 500.0, cooling_rate::Float64 = 0.99_99)
    soln_best, z_best = opt_function(soln, objective_function, perturb_function, max_iterations, temp_init, cooling_rate)
    return soln_best, z_best
end
