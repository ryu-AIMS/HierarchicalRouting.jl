
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

# Returns
- `ordered_centroids` : Centroid sequence DataFrame (containing cluster_id, lat, lon).
- `total_distance` : Total distance of the route.
- `dist_matrix` : Distance matrix between centroids.
"""
function nearest_neighbour(cluster_centroids::DataFrame)
    dist_matrix = distance_matrix(cluster_centroids)

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

function nearest_neighbour(dist_matrix::Matrix{Float64})
    num_pts = size(dist_matrix, 1) - 1  # excludes the depot
    visited = falses(num_pts + 1)
    total_distance = 0.0
    current_location = 1
    tour = Vector{Int}(undef, num_pts + 2)

    tour[1] = current_location
    visited[current_location] = true
    tour_idx = 2

    while tour_idx <= num_pts + 1 #length(tour) <= num_pts
        distances = dist_matrix[current_location, :]
        distances[visited] .= Inf
        nearest_idx = argmin(distances)

        total_distance += dist_matrix[current_location, nearest_idx]

        tour[tour_idx] = nearest_idx
        visited[nearest_idx] = true
        current_location = nearest_idx
        tour_idx += 1
    end

    # Return to the depot
    tour[tour_idx] = 1 # push!(tour, 1)
    total_distance += dist_matrix[current_location, 1]

    # Adjust waypt_sequence to zero-based indexing
    waypt_sequence = tour .- 1

    return waypt_sequence, total_distance
end

"""
    distance_matrix(cluster_centroids::DataFrame)

Calculate the haversine distance matrix between cluster centroids.

# Arguments
- `cluster_centroids` : DataFrame containing cluster_id, lat, lon.

# Returns
Distance matrix between cluster centroids.
"""
function distance_matrix(cluster_centroids::DataFrame)
    num_centroids = nrow(cluster_centroids)
    dist_matrix = Matrix{Float64}(undef, num_centroids, num_centroids)

    centroid_coords = [(row.lat, row.lon) for row in eachrow(cluster_centroids)]

    for i in 1:num_centroids
        for j in 1:num_centroids
            dist_matrix[i, j] = euclidean(centroid_coords[i], centroid_coords[j])
        end
    end

    return dist_matrix
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
function get_waypoints(sequence::DataFrame)::Vector{Point{2, Float64}}
    cluster_seq_ids = sequence.cluster_id
    n_cluster_seqs = length(cluster_seq_ids)

    waypoints = Vector{Tuple{Float64, Float64}}(undef, n_cluster_seqs+1)

    waypoints[1] = (sequence.lat[1], sequence.lon[1])

    for i in 1:(n_cluster_seqs - 1)
        current_centroid = first(sequence[sequence.cluster_id .== cluster_seq_ids[i], :])
        next_centroid = first(sequence[sequence.cluster_id .== cluster_seq_ids[i+1], :])

        centroid = (
            (current_centroid.lat + next_centroid.lat) / 2,
            (current_centroid.lon + next_centroid.lon) / 2
        )
        waypoints[i+1] = centroid
    end

    waypoints[n_cluster_seqs+1] = (sequence.lat[n_cluster_seqs], sequence.lon[n_cluster_seqs])

    return [GeometryBasics.Point{2, Float64}(wp...) for wp in waypoints]
end

"""
    two_opt(waypoints::Vector{Point{2, Float64}}, dist_matrix::Matrix{Float64})

Apply 2-opt heuristic to improve the current route (by uncrossing crossed links).

# Arguments
- `waypoints` : Vector of lat long tuples to be visited. Depot is the first and last point.
- `dist_matrix` : Distance matrix between waypoints. Depot is the first, but not last point.

# Returns
- `best_route` : Improved return route to/from depot.
- `best_distance` : Total distance of the best route.
"""
function two_opt(waypoints::Vector{Point{2, Float64}}, dist_matrix::Matrix{Float64})
    # Initialize route as ordered waypoints
    # Remove the last point (depot) from route
    best_route = collect(1:length(waypoints)-1)
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

    return best_route, best_distance
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
