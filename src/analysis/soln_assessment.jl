
function initial_solution(clusts::Vector{Cluster}, ms_exclusions::DataFrame, t_exclusions::DataFrame)

    config = TOML.parsefile(joinpath("src",".config.toml"))

    depot = Point{2, Float64}(config["parameters"]["depot_x"], config["parameters"]["depot_y"]) #(0.0, 0.0)
    n_tenders = config["parameters"]["n_tenders"]
    t_cap = config["parameters"]["t_cap"]


    cluster_centroids_df = DataFrame(id = Int[0], lon = Float64[depot[1]], lat = Float64[depot[2]]) # DataFrame(id = Int[], lon = Float64[], lat = Float64[]) # push!(cluster_centroids_df, (0, depot[1], depot[2]))
    [push!(cluster_centroids_df, (i, clust.attributes.centroid[1], clust.attributes.centroid[2])) for (i, clust) in enumerate(clusts)] #clusters)]
    #####################
    plot_centroids_and_exclusions(cluster_centroids_df, ms_exclusions, 10)

    # Nearest Neighbour to generate initial mothership route & matrix
    ms_soln_NN, ms_feasible_matrix = nearest_neighbour(cluster_centroids_df, ms_exclusions)
    # plot_waypoints_and_exclusions(ms_soln_NN.route, ms_soln_NN.cluster_sequence, ms_exclusions, 10)

    # 2-opt to improve the NN soln
    ms_soln_2opt = two_opt(ms_soln_NN.cluster_sequence, ms_feasible_matrix)
    plot_waypoints_and_exclusions(ms_soln_2opt.route, ms_soln_2opt.cluster_sequence, ms_exclusions, 10)

    #Refs to `clusts` rather than `clusters` because `clusts` have sub-cluster nodes
    clust_seq = [i for i in ms_soln_2opt.cluster_sequence.id if i!==0 && i <= length(clusts)]
    tender_soln_NN = HierarchicalRouting.ClusterSolution[]

    for (i, cluster_id) in enumerate(clust_seq)
        start_waypoint =  ms_soln_2opt.route.waypoint[2 * i]
        end_waypoint =  ms_soln_2opt.route.waypoint[2 * i + 1]

        println("$(i): Clust $(cluster_id) from $(start_waypoint) to $(end_waypoint)")

        t_solution = tender_sequential_nearest_neighbour(
            clusts[cluster_id],
            (
                start_waypoint,
                end_waypoint
                ),
            n_tenders, t_cap, t_exclusions)

        push!(tender_soln_NN, t_solution[1])
    end

    # plot_tender_routes(tender_soln_NN, ms_soln_2opt.route.waypoint)


    # mothership = ClusterSolution(0, [clusters[1]], 0.0)
    return MSTSolution(ms_soln_2opt, tender_soln_NN) # MSTSolution(mothership, clusters)
end
function initial_solution()
    clusters_raster, ms_exclusions, t_exclusions = load_problem()
    clusts = create_clusters(clusters_raster)

    return initial_solution(clusts, ms_exclusions, t_exclusions)
end

function total_dist(soln::MSTSolution)
    return soln.mothership.cost + sum([cluster.cost for cluster in soln.clusters])
end

function critical_path(soln::MSTSolution)
    longest_sortie = [maximum([sortie.cost for sortie in cluster.sorties]) for cluster in soln.clusters]
    tender_cost = sum(longest_sortie)
    return tender_cost + soln.mothership.cost
end