
"""
    return_route_distance(route::Vector{Int64}, dist_matrix::Matrix{Float64})

Calculate the total distance of a route starting from index 1, and returning to index 1.

# Arguments
- `route` : Vector of cluster indices.
- `dist_matrix` : Distance matrix between clusters.

# Returns
Total distance of the return route.
"""
function return_route_distance(route::Vector{Int64}, dist_matrix::Matrix{Float64})
    total_dist = 0.0
    n = length(route)  # Adjust for duplicate last point as the start point)

    total_dist = sum([dist_matrix[route[i], route[i + 1]] for i in 1:n-1])

    total_dist += dist_matrix[route[n], route[1]]
    return total_dist
end

"""
    sortie_cost(sortie::Vector{Vector{Int64}}, dist_matrix::Matrix{Float64})

Calculate the total distance of all sorties.

# Arguments
- `sorties` : Vector of Vector of node indices.
- `dist_matrix` : Distance matrix between nodes.

# Returns
Vector of total distance of each sortie.
"""
function sortie_cost(sorties::Vector{Vector{Int64}}, dist_matrix::Matrix{Float64})::Vector{Float64}
    # TODO: Handling for empty tender tours
    sortie_dist = [dist_matrix[1, tour[1]] for tour in sorties]
    sortie_dist .+= [length(tour) > 1 ? sum([dist_matrix[tour[i], tour[i+1]] for i in 1:(length(tour)-1)]) : 0 for tour in sorties]
    sortie_dist .+= [dist_matrix[tour[end], 1] for tour in sorties]
    return sortie_dist
end

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
