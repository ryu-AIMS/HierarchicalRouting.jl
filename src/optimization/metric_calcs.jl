
"""
    tender_clust_cost(cluster::ClusterSolution)::Vector{Float64}

Compute the cost of each sortie in a cluster.

# Arguments
- `cluster::ClusterSolution` : The cluster solution.

# Returns
- `sortie_dist` : The cost of each sortie in the cluster.
"""
function tender_clust_cost(tenders::TenderSolution)::Vector{Float64}
    sortie_dist = [euclidean(tenders.start, sortie.nodes[1]) for sortie in tenders.sorties] # haversine

    sortie_dist .+= sum([
        length(sortie.nodes) > 1 ?
        sum(euclidean.(sortie.nodes[1:end-1], sortie.nodes[2:end])) :
        0.0
        for sortie in tenders.sorties
    ])

    sortie_dist .+= euclidean.([sortie.nodes[end] for sortie in tenders.sorties], tenders.finish)
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

    longest_sortie = [maximum(tender_clust_cost(cluster)) for cluster in soln.tenders]
    mothership_sub_clust = mothership_cost_within_clusts(soln)

    cluster_cost = [max(longest_sortie[i], mothership_sub_clust[i]) for i in 1:length(longest_sortie)]

    return sum(cluster_cost) + mothership_cost_between_clusts(soln)
end
