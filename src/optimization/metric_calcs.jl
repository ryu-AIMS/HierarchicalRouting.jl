
"""
    return_route_distance(route::Vector{Int64}, dist_matrix::Matrix{Float64})::Float64

Calculate the total distance of a route starting from index 1, and returning to index 1.

# Arguments
- `route`: Vector of cluster indices.
- `dist_matrix`: Distance matrix between clusters.

# Returns
Total distance of the return route.
"""
function return_route_distance(route::Vector{Int64}, dist_matrix::Matrix{Float64})::Float64
    return sum(getindex.(Ref(dist_matrix), route[1:end-1], route[2:end])) +
           dist_matrix[route[end], route[1]]  # Return to start
end

"""
    sortie_dist(
        sortie::Vector{Vector{Int64}},
        dist_matrix::Matrix{Float64}
    )::Vector{Float64}

Calculate the total distance of all sorties.

# Arguments
- `sorties`: Vector of Vector of node indices.
- `dist_matrix`: Distance matrix between nodes.

# Returns
Vector of total distance of each sortie.
"""
function sortie_dist(
    sorties::Vector{Vector{Int64}},
    dist_matrix::Matrix{Float64}
)::Vector{Float64}
    sorties = filter(!isempty, sorties)
    n = length(sorties)
    total_dist = Vector{Float64}(undef, n)

    for i in 1:n
        tour = sorties[i]
        m = length(tour)

        dist = dist_matrix[1, tour[1]] # dist from depot (node 1) to first node
        dist += m > 1 ? sum(getindex.(Ref(dist_matrix), tour[1:m-1], tour[2:m])) : 0.0
        dist += dist_matrix[tour[end], 1] # dist from last node back to depot

        total_dist[i] = dist
    end
    return total_dist
end

"""
    tender_clust_dist(tenders::TenderSolution)::Vector{Float64}

Compute the cost of each sortie in a cluster.

# Arguments
- `tenders`: Tender solution.

# Returns
The cost of each sortie in the cluster.
"""
function tender_clust_dist(tenders::TenderSolution)::Vector{Float64}
    tender_matrices = getfield.(tenders.sorties, :dist_matrix)
    sortie_dists = sum.(tender_matrices)
    return sortie_dists
end

"""
    tender_sortie_dist(sortie::Route)::Float64
    tender_sortie_dist(
        node_order::Vector{Int64},
        dist_matrix::Matrix{Float64}
    )::Float64

Compute the cost of a sortie, which starts and ends at given points, but does not return to
start.

# Arguments
- `sortie`: Route object containing nodes and distance matrix.
- `node_order`: Vector of node indices representing the order of the sortie.
- `dist_matrix`: Distance matrix between nodes, ordered by node index.

# Returns
The total distance of the sortie.
"""
function tender_sortie_dist(sortie::Route)::Float64
    return sum(@view sortie.dist_matrix[1:end-1, 2:end])
end
function tender_sortie_dist(
    node_order::Vector{Int64},
    dist_matrix::Matrix{Float64}
)::Float64
    return sum(
        dist_matrix[node_order[i], node_order[i+1]]
        for i in 1:(length(node_order)-1)
    )
end

"""
    mothership_dist_between_clusts(route::Route, num_clusters::Int=0)::Float64

Compute the cost of the mothership route between clusters, not including across each cluster.
Optionally, the number of clusters can be specified to limit the calculation to the first
`num_clusters` clusters.

# Arguments
- `route`: Full mothership route between waypoints.

# Returns
- The sum of mothership distances between clusters.
"""
function mothership_dist_between_clusts(route::Route, num_clusters::Int=0)::Float64
    full_route = route.dist_matrix
    n = iszero(num_clusters) ? length(full_route) : 2 * num_clusters + 1
    return sum(full_route[1:2:n])
end

"""
    mothership_dist_within_clusts(route::Route)::Vector{Float64}

Compute the cost of the mothership within each cluster, not including between clusters.

# Arguments
- `route`: Full mothership route between waypoints.

# Returns
- The mothership distance within each cluster.
"""
function mothership_dist_within_clusts(route::Route, num_clusters::Int=0)::Vector{Float64}
    full_route = route.dist_matrix
    n = iszero(num_clusters) ? length(full_route) : 2 * num_clusters + 1
    return full_route[2:2:n]
end

"""
    critical_path(
        soln::MSTSolution,
        vessel_weightings::NTuple{2,AbstractFloat}=(1.0, 1.0)
    )::Float64

Compute the critical path cost of the solution.
This is the longest path for a return trip, quantified as the sum of:
- the sum (`cluster_cost_total`) of the longest path within each cluster (`cluster_cost_each`),
    i.e., the sum of the maximum(sortie, mothership) cost within each cluster, and
- the sum of the mothership cost between clusters: `tow_cost`.

# Arguments
- `soln`: The solution to evaluate.
- `vessel_weightings`: The weightings for mothership and sortie costs.

# Returns
The total (critical path) cost of the solution.
"""
function critical_path(
    soln::MSTSolution,
    vessel_weightings::NTuple{2,AbstractFloat}=(1.0, 1.0)
)::Float64
    tenders = soln.tenders[end]
    ms_route = soln.mothership_routes[end].route
    num_clusters = length(tenders)

    # Within clusters
    cluster_sorties = tender_clust_dist.(tenders)
    cluster_sorties = map(x -> isempty(x) ? [0.0] : x, cluster_sorties)
    longest_sortie_cost = maximum.(cluster_sorties) .* vessel_weightings[2]
    mothership_within_clusts = mothership_dist_within_clusts(ms_route)[1:num_clusters]
    mothership_sub_clust_cost = vessel_weightings[1] * mothership_within_clusts

    cluster_cost_each = max.(longest_sortie_cost, mothership_sub_clust_cost)
    cluster_cost_total = sum(cluster_cost_each)

    # Between clusters
    tow_cost = vessel_weightings[1] * mothership_dist_between_clusts(ms_route, num_clusters)

    return cluster_cost_total + tow_cost
end

function critical_distance_path(
    soln::MSTSolution,
    vessel_weightings::NTuple{2,AbstractFloat}=(1.0, 1.0)
)::Float64
    tenders = soln.tenders[end]
    ms_route = soln.mothership_routes[end].route
    num_clusters = length(tenders)

    # Within clusters
    cluster_sorties = tender_clust_dist.(tenders)
    cluster_sorties = map(x -> isempty(x) ? [0.0] : x, cluster_sorties)
    total_sortie_cost = sum.(cluster_sorties) .* vessel_weightings[2]
    mothership_within_clusts = mothership_dist_within_clusts(ms_route)[1:num_clusters]
    mothership_sub_clust_cost = vessel_weightings[1] * mothership_within_clusts

    cluster_cost_total = sum(vcat(total_sortie_cost, mothership_sub_clust_cost))

    # Between clusters
    tow_cost = vessel_weightings[1] * mothership_dist_between_clusts(ms_route, num_clusters)

    return cluster_cost_total + tow_cost
end
