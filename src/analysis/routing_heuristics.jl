
"""
    create_exclusion_zones(target_bathy::Raster, ms_depth)

Create exclusion zones based on environmental raster data and vessel threshold.

# Arguments
- `env_constraint` : Environmental constraint raster.
- `threshold` : Threshold for given vessel's environmental constraint.

# Returns
Exclusion zones for environmental constraint and vessel threshold provided.
"""
function create_exclusion_zones(env_constraint::Raster, threshold::Float64)
    exclusion_zones = env_constraint .>= threshold
    return exclusion_zones
end

"""
    nearest_neighbour(dist_matrix::Matrix{Float64})

Apply the nearest neighbor algorithm starting from the depot (1st row/col) and returning to the depot.

# Arguments
- `cluster_centroids` : DataFrame containing cluster_id, lat, lon. Depot is cluster 0 in row 1.
- `exclusions` : DataFrame containing exclusion zones.
# Returns
- `ordered_centroids` : Centroid sequence DataFrame (containing cluster_id, lat, lon).
- `total_distance` : Total distance of the route.
- `dist_matrix` : Distance matrix between centroids.
"""
function nearest_neighbour(cluster_centroids::DataFrame, exclusions::DataFrame)
    dist_matrix = get_feasible_matrix([Point(row.lat, row.lon) for row in eachrow(cluster_centroids)], exclusions)

    num_clusters = size(dist_matrix, 1) - 1  # excludes the depot
    visited = falses(num_clusters + 1)
    tour = Int[]
    total_distance = 0.0

    current_location = 1
    push!(tour, current_location)
    visited[current_location] = true

    while length(tour) <= num_clusters
        distances = dist_matrix[current_location, :]
        distances[visited] .= Inf
        nearest_idx = argmin(distances)

        push!(tour, nearest_idx)
        total_distance += dist_matrix[current_location, nearest_idx]
        visited[nearest_idx] = true
        current_location = nearest_idx
    end

    # Return to the depot
    push!(tour, 1)
    total_distance += dist_matrix[current_location, 1]

    # Adjust cluster_sequence to zero-based indexing
    cluster_sequence = tour .- 1

    ordered_centroids = cluster_centroids[[findfirst(==(id), cluster_centroids.cluster_id) for id in cluster_sequence], :]

    return ordered_centroids, total_distance, dist_matrix
end

"""
    get_waypoints(centroid_sequence::DataFrame)::Vector{Point{2, Float64}

Calculate mothership waypoints between sequential clusters.
Based on calc: midpoint of current and next cluster.

# Arguments
- `sequence` :  centroid lat long coordinates in sequence; including depot as the first and last cluster.

# Returns
Route as vector of lat long points. Depot is first and last waypoints.
"""
#TODO: Implement convex hull exclusion zones
function get_waypoints(sequence::DataFrame)::DataFrame
    n_cluster_seqs = nrow(sequence)

    waypoints = Vector{Point{2, Float64}}(undef, n_cluster_seqs+1)
    connecting_clusters = Vector{NTuple{2, Int64}}(undef, n_cluster_seqs+1)

    waypoints[1] = (sequence.lat[1], sequence.lon[1])
    connecting_clusters[1] = (sequence.cluster_id[1], sequence.cluster_id[1])

    for i in 1:(n_cluster_seqs - 1)
        current_lat, current_lon, current_clust = sequence.lat[i], sequence.lon[i], sequence.cluster_id[i] #first(sequence[sequence.cluster_id .== cluster_seq_ids[i], :])
        next_lat, next_lon, next_clust = sequence.lat[i+1], sequence.lon[i+1], sequence.cluster_id[i+1] #first(sequence[sequence.cluster_id .== cluster_seq_ids[i+1], :])

        centroid = (
            (current_lat + next_lat) / 2,
            (current_lon + next_lon) / 2
        )
        waypoints[i+1] = centroid
        connecting_clusters[i+1] = (current_clust, next_clust)
    end

    waypoints[n_cluster_seqs+1] = (sequence.lat[end], sequence.lon[end])
    connecting_clusters[n_cluster_seqs+1] = (sequence.cluster_id[end], sequence.cluster_id[end])

    centroid_df = DataFrame(
        waypoint = waypoints,
        connecting_clusters = connecting_clusters
    )
    return centroid_df
end

"""
    two_opt(waypoints::Vector{Point{2, Float64}}, dist_matrix::Matrix{Float64})

Apply 2-opt heuristic to improve the current route (by uncrossing crossed links).

# Arguments
- `cluster_centroids` : DataFrame (cluster_id, lat, lon) incl depot as cluster 0 in row 1. NOTE: Remove depot as final pt
- `dist_matrix` : Distance matrix between waypoints. Depot is the first, but not last point.

# Returns
- `ordered_centroids` : Improved return route to/from depot.
- `best_distance` : Total distance of the best route.
"""
function two_opt(cluster_centroids::DataFrame, dist_matrix::Matrix{Float64})
    points = [Point(row.lat, row.lon) for row in eachrow(cluster_centroids)]
    # Initialize route as ordered waypoints
    # Remove the last point (depot) from route
    best_route = [row.cluster_id+1 for row in eachrow(cluster_centroids)]
    best_distance = return_route_distance(best_route, dist_matrix)
    improved = true

    while improved
        improved = false
        for i in 1:(length(best_route) - 1)
            for j in (i + 1):length(best_route)
                new_route = two_opt_swap(best_route, i, j)
                new_distance = return_route_distance(new_route, dist_matrix)

                if new_distance < best_distance
                    best_route = new_route
                    best_distance = new_distance
                    improved = true
                end
            end
        end
    end

    # Re-orient route to start from the depot (1), and add the depot as final point
    best_route = orient_route(best_route)
    push!(best_route, best_route[1])

    best_route .-= 1

    ordered_centroids = cluster_centroids[[findfirst(==(id), cluster_centroids.cluster_id) for id in best_route], :]

    return ordered_centroids, best_distance
end

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
    two_opt_swap(route::Vector{Int64}, i::Int, j::Int)

Swap two points in a route.

# Arguments
- `route` : Vector of cluster indices.
- `i` : Index of the first point to swap.
- `j` : Index of the second point to swap.

# Returns
Route with swapped points.
"""
function two_opt_swap(route::Vector{Int64}, i::Int, j::Int)
    return vcat(route[1:(i-1)], reverse(route[i:j]), route[(j+1):end])
end

"""
    orient_route(route::Vector{Int64})

Orient the route such that it starts with the depot (1).

# Arguments
- `route` : Vector of cluster indices.

# Returns
Route starting with the depot.
"""
function orient_route(route::Vector{Int64})
    idx = findfirst(==(1), route)
    return vcat(route[idx:end], route[1:idx-1])
end
