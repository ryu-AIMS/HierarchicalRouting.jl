
function initial_solution(clusts::Vector{Cluster}, ms_exclusions::DataFrame, t_exclusions::DataFrame)

    config = TOML.parsefile(joinpath("src",".config.toml"))

    depot = Point{2, Float64}(config["parameters"]["depot_x"], config["parameters"]["depot_y"]) #(0.0, 0.0)
    n_tenders = config["parameters"]["n_tenders"]
    t_cap = config["parameters"]["t_cap"]


    cluster_centroids_df = DataFrame(id = Int[0], lon = Float64[depot[1]], lat = Float64[depot[2]]) # DataFrame(id = Int[], lon = Float64[], lat = Float64[]) # push!(cluster_centroids_df, (0, depot[1], depot[2]))
    [push!(cluster_centroids_df, (i, clust.attributes.centroid[1], clust.attributes.centroid[2])) for (i, clust) in enumerate(clusts)] #clusters)]
    #####################
    # plot_centroids_and_exclusions(cluster_centroids_df, ms_exclusions, 10)

    # Nearest Neighbour to generate initial mothership route & matrix
    ms_soln_NN, ms_feasible_matrix = nearest_neighbour(cluster_centroids_df, ms_exclusions)
    # plot_waypoints_and_exclusions(ms_soln_NN.route, ms_soln_NN.cluster_sequence, ms_exclusions) #, 10)

    # 2-opt to improve the NN soln
    ms_soln_2opt = two_opt(ms_soln_NN.cluster_sequence, ms_feasible_matrix)
    # plot_waypoints_and_exclusions(ms_soln_2opt.route, ms_soln_2opt.cluster_sequence, ms_exclusions)#, 10)

    #Refs to `clusts` rather than `clusters` because `clusts` have sub-cluster nodes
    clust_seq = [i for i in ms_soln_2opt.cluster_sequence.id if i!==0 && i <= length(clusts)]
    tender_soln_NN = HierarchicalRouting.ClusterSolution[]

    for (i, cluster_id) in enumerate(clust_seq)
        start_waypoint =  ms_soln_2opt.route.waypoint[2 * i]
        end_waypoint =  ms_soln_2opt.route.waypoint[2 * i + 1]
        @info "$(i): Clust $(cluster_id) from $(start_waypoint) to $(end_waypoint)"

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
function initial_solution(target_scenario::String="")
    clusters_raster, ms_exclusions, t_exclusions = load_problem(target_scenario)
    clusts = create_clusters(clusters_raster)

    return initial_solution(clusts, ms_exclusions, t_exclusions), clusts, ms_exclusions, t_exclusions
end

"""
    tender_clust_cost(cluster::ClusterSolution)::Vector{Float64}

Compute the cost of each sortie in a cluster.

# Arguments
- `cluster::ClusterSolution` : The cluster solution.

# Returns
- `sortie_dist` : The cost of each sortie in the cluster.
"""
function tender_clust_cost(cluster::ClusterSolution)::Vector{Float64}
    sortie_dist = [euclidean(cluster.start, sortie.nodes[1]) for sortie in cluster.sorties]

    for sortie in cluster.sorties
        if length(sortie.nodes) > 1
            sortie_dist .+= sum([euclidean(sortie.nodes[s], sortie.nodes[s+1]) for s in 1:(length(sortie.nodes)-1)])
        end
    end
    # TODO: Why does this get different result??
    # sortie_dist .+= [
    #     length(sortie.nodes) > 1 ? sum([euclidean(sortie.nodes[s], sortie.nodes[s+1]) for s in 1:(length(sortie.nodes)-1)]) : 0.0
    #     for sortie in cluster.sorties
    # ]

    sortie_dist .+= [euclidean(cluster.sorties[i].nodes[end], cluster.finish) for i in 1:length(cluster.sorties)]
    return sortie_dist
end

function total_dist(soln::MSTSolution)
    cluster_dist = sum([sum([euclidean(sortie.nodes[i], sortie.nodes[i+1]) for i in 1:(length(sortie.nodes)-1)]) for sortie in soln.mothership.sorties])
    return soln.mothership.cost + cluster_dist # sum([cluster.cost for cluster in soln.clusters])
end

"""
    mothership_cost_between_clusts(soln::MSTSolution)

Compute the cost of the mothership between clusters, not including across each cluster.

# Arguments
- `soln::MSTSolution` : Full MST solution.

# Returns
- The sum of (euclidean) mothership distances between clusters.
"""
function mothership_cost_between_clusts(soln::MSTSolution)
    return sum(euclidean(soln.mothership.route.waypoint[i], soln.mothership.route.waypoint[i+1]) for i in 1:2:(nrow(soln.mothership.route)))
end

"""
    mothership_cost_within_clusts(soln::MSTSolution)

Compute the cost of the mothership within each cluster, not including between clusters.

# Arguments
- `soln::MSTSolution` : Full MST solution.

# Returns
- The sum of (euclidean) mothership distances within each cluster.
"""
function mothership_cost_within_clusts(soln::MSTSolution)
    return [euclidean(soln.mothership.route.waypoint[i], soln.mothership.route.waypoint[i+1]) for i in 2:2:(nrow(soln.mothership.route)-1)]
end

function critical_path(soln::MSTSolution)

    longest_sortie = [maximum(tender_clust_cost(cluster)) for cluster in soln.clusters]
    mothership_sub_clust = mothership_cost_within_clusts(soln)

    cluster_cost = [max(longest_sortie[i], mothership_sub_clust[i]) for i in 1:length(longest_sortie)]

    return sum(cluster_cost) + mothership_cost_between_clusts(soln)
end

function perturb_solution(soln::MSTSolution, clust_idx::Int=-1)
    new_soln = deepcopy(soln)

    if clust_idx == -1
        # Choose ONE random cluster - assume fixed clustering
        clust_idx = rand(1:length(new_soln.clusters))
    end
    cluster = new_soln.clusters[clust_idx]

    if length(cluster.sorties) < 2
        return new_soln # No perturbation possible if a cluster has less than 2 sorties
    end

    # Choose TWO random sorties from the cluster
    sortie_a_idx, sortie_b_idx = rand(1:length(cluster.sorties), 2)
    # TODO: Allow same sortie if not applying two-opt
    while sortie_a_idx == sortie_b_idx
        sortie_b_idx = rand(1:length(cluster.sorties))
    end
    sortie_a, sortie_b = cluster.sorties[sortie_a_idx], cluster.sorties[sortie_b_idx]

    if isempty(sortie_a.nodes) || isempty(sortie_b.nodes)
        return new_soln # No perturbation possible if a sortie has no nodes
    end

    node_a_idx, node_b_idx = rand(1:length(sortie_a.nodes)), rand(1:length(sortie_b.nodes))
    node_a, node_b = sortie_a.nodes[node_a_idx], sortie_b.nodes[node_b_idx]

    # Swap the nodes between the two sorties
    new_soln.clusters[clust_idx].sorties[sortie_a_idx].nodes[node_a_idx] = node_b
    new_soln.clusters[clust_idx].sorties[sortie_b_idx].nodes[node_b_idx] = node_a

    # TODO:Re-run two-opt on the modified sorties

    # TODO: Recompute the cost for each modified sortie
    # This assumes a function `compute_sortie_cost` exists
    # sortie1.cost = compute_sortie_cost(sortie1)
    # sortie2.cost = compute_sortie_cost(sortie2)

    # Update the tender's total cost if necessary
    # tender.cost = sum(sortie.cost for sortie in tender.sorties)

    return new_soln
end

"""
    simulated_annealing(
        soln_init::MSTSolution,
        objective_function::Function,
        perturb_function::Function,
        max_iterations::Int = 100_000,
        temp_init::Float64 = 50.0,
        cooling_rate::Float64 = 0.99_99
    )

Simulated Annealing optimization algorithm to optimize the solution.

# Arguments
- `soln_init` : Initial solution.
- `objective_function` : Function to evaluate the solution.
- `perturb_function` : Function to perturb the solution.
- `max_iterations` : Maximum number of iterations. Default = 10_000.
- `temp_init` : Initial temperature. Default = 500.0.
- `cooling_rate` : Rate of cooling to guide acceptance probability for SA algorithm. Default = 0.995 = 99.95%.
- `static_limit` : Number of iterations to allow stagnation before early exit. Default = 200.

# Returns
- `soln_best` : Best solution::MSTSolution found.
- `obj_best` : Objective value of the best solution.
"""
function simulated_annealing(
    soln_init::MSTSolution,
    objective_function::Function,
    perturb_function::Function,
    max_iterations::Int = 10_000,
    temp_init::Float64 = 500.0,
    cooling_rate::Float64 = 0.995,
    static_limit::Int = 200
)

    # Initialize best solution as initial
    soln_best = soln_init
    obj_best = objective_function(soln_init)

    for clust_idx in 1:length(soln_init.clusters)
        # Initialize current solution as best, reset temp
        temp = temp_init
        soln_current = soln_best
        obj_current = obj_best
        static_ctr = 0

        @info "\nCluster $(soln_init.clusters[clust_idx].id)"
        @info "Iteration \tBest Value \t\tTemp"
        @info "0\t\t$obj_best\t$temp"

        for iteration in 1:max_iterations
            soln_proposed = perturb_function(soln_current, clust_idx)
            obj_proposed = objective_function(soln_proposed)

            static_ctr += 1
            improvement = obj_current - obj_proposed

            # If the new solution is improved OR meets acceptance prob criteria
            if improvement > 0 || exp(improvement / temp) > rand()
                soln_current = soln_proposed
                obj_current = obj_proposed

                if obj_current < obj_best
                    static_ctr = 0
                    soln_best = soln_current
                    obj_best = obj_current
                    @info "$iteration\t\t$obj_best\t$temp"
                end
            end

            temp *= cooling_rate

            if static_ctr > static_limit
                @info "Early exit at iteration $iteration due to stagnation"
                break
            end
        end
    end

    @info "Final Value: $obj_best"
    @info "Δ: $(objective_function(soln_init) - obj_best)"
    return soln_best, obj_best
end
