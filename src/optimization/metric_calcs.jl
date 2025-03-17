
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
    sortie_cost(
        sortie::Vector{Vector{Int64}},
        dist_matrix::Matrix{Float64}
    )::Vector{Float64}

Calculate the total distance of all sorties.

# Arguments
- `sorties` : Vector of Vector of node indices.
- `dist_matrix` : Distance matrix between nodes.

# Returns
Vector of total distance of each sortie.
"""
function sortie_cost(
    sorties::Vector{Vector{Int64}},
    dist_matrix::Matrix{Float64}
)::Vector{Float64}
    # TODO: Handling for empty tender tours
    sortie_dist = [dist_matrix[1, tour[1]] for tour in sorties]
    sortie_dist .+= [
        length(tour) > 1 ?
            sum([dist_matrix[tour[i], tour[i+1]] for i in 1:(length(tour)-1)]) :
            0
        for tour in sorties
    ]
    sortie_dist .+= [dist_matrix[tour[end], 1] for tour in sorties]
    return sortie_dist
end

"""
    tender_clust_dist(tenders::TenderSolution)::Vector{Float64}

Compute the cost of each sortie in a cluster.

# Arguments
- `tenders::TenderSolution` : Tender solution.

# Returns
The cost of each sortie in the cluster.
"""
function tender_clust_dist(tenders::TenderSolution)::Vector{Float64}
    sortie_dist = [euclidean(tenders.start, sortie.nodes[1]) for sortie in tenders.sorties]

    sortie_dist .+= sum([
        length(sortie.nodes) > 1 ?
            sum(euclidean.(sortie.nodes[1:end-1], sortie.nodes[2:end])) :
            0.0
        for sortie in tenders.sorties
    ])

    sortie_dist .+= euclidean.(
        [sortie.nodes[end] for sortie in tenders.sorties],
        tenders.finish
    )
    return sortie_dist
end

"""
    mothership_dist_between_clusts(route::Route)

Compute the cost of the mothership route between clusters, not including across each cluster.

# Arguments
- `route::Route` : Full mothership route between waypoints.

# Returns
The sum of (euclidean) mothership distances between clusters.
"""
function mothership_dist_between_clusts(route::Route)
    return sum(euclidean(route.nodes[i], route.nodes[i+1]) for i in 1:2:(length(route.nodes)))
end

"""
    mothership_dist_within_clusts(route::Route)

Compute the cost of the mothership within each cluster, not including between clusters.

# Arguments
- `route::Route` : Full mothership route between waypoints.

# Returns
The sum of (euclidean) mothership distances across/within each cluster.
"""
function mothership_dist_within_clusts(route::Route)
    return [euclidean(route.nodes[i], route.nodes[i+1]) for i in 2:2:(length(route.nodes)-1)]
end

"""
    critical_path(soln::MSTSolution, vessel_weightings::NTuple{2, Float64}=(1.0, 1.0))

Compute the critical path cost of the solution.
This is the longest path for a return trip, quantified as the sum of:
    - the longest path within each cluster - `cluster_cost_each`,
        i.e., the maximum(sortie, mothership) cost within each cluster
        summed for `cluster_cost`, and
    - the sum of the mothership cost between clusters - `tow_cost`.

# Arguments
- `soln::MSTSolution` : The solution to evaluate.
- `vessel_weightings::NTuple{2, Float64}` : The weightings for mothership and sortie costs.

# Returns
The total (critical path) cost of the solution.
"""
function critical_path(soln::MSTSolution, vessel_weightings::NTuple{2, Float64}=(1.0, 1.0))

    # Within clusters
    longest_sortie_cost = maximum.(tender_clust_dist.(soln.tenders)) .* vessel_weightings[2]
    mothership_sub_clust_cost = mothership_dist_within_clusts(soln.mothership.route) * vessel_weightings[1]

    cluster_cost_each = max.(longest_sortie_cost, mothership_sub_clust_cost)
    cluster_cost = sum(cluster_cost_each)

    # Between clusters
    tow_cost = (mothership_dist_between_clusts(soln.mothership.route) * vessel_weightings[1])

    return cluster_cost + tow_cost
end
