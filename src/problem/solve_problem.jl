
function initial_solution(clusts::Vector{Cluster}, ms_exclusions::DataFrame, t_exclusions::DataFrame)

    config = TOML.parsefile(joinpath("src",".config.toml"))

    depot = Point{2, Float64}(config["parameters"]["depot_x"], config["parameters"]["depot_y"]) #(0.0, 0.0)
    n_tenders = config["parameters"]["n_tenders"]
    t_cap = config["parameters"]["t_cap"]

    cluster_centroids_df = DataFrame(id = Int[0], lon = Float64[depot[1]], lat = Float64[depot[2]]) # DataFrame(id = Int[], lon = Float64[], lat = Float64[]) # push!(cluster_centroids_df, (0, depot[1], depot[2]))
    [push!(cluster_centroids_df, (i, clust.attributes.centroid[1], clust.attributes.centroid[2])) for (i, clust) in enumerate(clusts)] #clusters)]

    # Nearest Neighbour to generate initial mothership route & matrix
    ms_soln_NN, ms_feasible_matrix = nearest_neighbour(cluster_centroids_df, ms_exclusions)

    # 2-opt to improve the NN soln
    ms_soln_2opt = two_opt(ms_soln_NN.cluster_sequence, ms_feasible_matrix)

    #Refs to `clusts` rather than `clusters` because `clusts` have sub-cluster nodes
    clust_seq = [i for i in ms_soln_2opt.cluster_sequence.id if i!==0 && i <= length(clusts)]
    tender_soln_NN = HierarchicalRouting.ClusterSolution[]

    for (i, cluster_id) in enumerate(clust_seq)
        start_waypoint =  ms_soln_2opt.route.waypoint[2 * i]
        end_waypoint =  ms_soln_2opt.route.waypoint[2 * i + 1]
        @info "$(i): Clust $(cluster_id) from $(start_waypoint) to $(end_waypoint)"

        t_solution = tender_sequential_nearest_neighbour(
            clusts[cluster_id],
            (start_waypoint, end_waypoint),
            n_tenders, t_cap, t_exclusions)

        push!(tender_soln_NN, t_solution[1])
    end

    # MSTSolution(mothership, clusters)
    return MSTSolution(ms_soln_2opt, tender_soln_NN)
end

function improve_solution(soln::MSTSolution, opt_function::Function, objective_function::Function, perturb_function::Function, max_iterations::Int = 100_000, temp_init::Float64 = 50.0, cooling_rate::Float64 = 0.99_99)
    soln_best, z_best = opt_function(soln, objective_function, perturb_function, max_iterations, temp_init, cooling_rate)
    return soln_best, z_best
end
