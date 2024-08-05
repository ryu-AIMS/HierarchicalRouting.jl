
function create_exclusion_zones(target_bathy::Raster, ms_depth)
    # Create exclusion zones based on the bathymetry data
    exclusion_zones = target_bathy .<= ms_depth
    return exclusion_zones
end

function is_feasible_path(start_pt::Tuple{Float64, Float64}, end_pt::Tuple{Float64, Float64}, env_constraint::DataFrame)
    # Create a line from start_pt to end_pt
    line = GI.LineString([start_pt, end_pt])

    # Check if the line intersects with any polygon in the env_constraint
    for row in eachrow(env_constraint)
        polygon = row.geometry
        if GI.intersects(line, polygon)
            return false
        end
    end
    return true
end

"""
    shortest_path(waypoints::Vector{Tuple{Float64, Float64}}, exclusions::Raster{Int16, 2})

Determine the shortest path between waypoints accounting for exclusion zones.

# Arguments
- `waypoints` : Path of waypoints to visit.
- `exclusions` : Exclusion zones of all environmental constraints.

# Returns
The identified shortest path.
"""
function shortest_path(waypoints::Vector{Tuple{Float64, Float64}}, exclusions::Raster{Int16, 2})
    # Initialize the total distance
    total_distance = 0.0

    # Iterate over the waypoints and calculate the distance between them
    for i in 1:length(waypoints)-1
        # Check if the path between the waypoints is feasible
        if !is_feasible_path(waypoints[i], waypoints[i+1], exclusions)
            # find shortest path around constraint
        else
            # Calculate the distance between the waypoints
            distance = haversine(waypoints[i], waypoints[i+1])
        end
        total_distance += distance
    end

    return total_distance
end

function distance_matrix(cluster_centroids::DataFrame)
    # Number of centroids
    num_centroids = nrow(cluster_centroids)

    # Initialize the distance matrix
    dist_matrix = Matrix{Float64}(undef, num_centroids, num_centroids)

    # Get the coordinates of the centroids
    centroid_coords = [(row.lat, row.lon) for row in eachrow(cluster_centroids)]

    # Compute distances between centroids
    for i in 1:num_centroids
        for j in 1:num_centroids
            dist_matrix[i, j] = haversine(centroid_coords[i], centroid_coords[j])
        end
    end

    return dist_matrix
end

"""
    nearest_neighbour(dist_matrix::Matrix{Float64})

Apply the nearest neighbor algorithm starting from the depot (1st row/col) and returning to the depot.
"""
function nearest_neighbour(cluster_centroids::DataFrame)

    dist_matrix = distance_matrix(cluster_centroids)

    num_clusters = size(dist_matrix, 1) - 1  # excludes the depot
    visited = falses(num_clusters + 1)
    tour = Int[]
    total_distance = 0.0

    # Start at the depot
    current_location = 1
    push!(tour, current_location)
    visited[current_location] = true

    while length(tour) <= num_clusters
        # Find the nearest unvisited neighbor
        distances = dist_matrix[current_location, :]
        distances[visited] .= Inf  # visited distances = inf
        nearest_idx = argmin(distances)

        # Update tour and total distance
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

    # Generate mothership waypoints
    # cluster_centroids = vcat(DataFrame(cluster_id = [0], lat = [depot[1]], lon = [depot[2]]), sort!(cluster_centroids, :cluster_id))
    ordered_centroids = cluster_centroids[[findfirst(==(id), cluster_centroids.cluster_id) for id in cluster_sequence], :]
    # mothership_waypoints = [(row.lat, row.lon) for row in eachrow(ordered_centroids)]

    return ordered_centroids, total_distance, dist_matrix
end

"""
    get_waypoints(centroid_sequence::DataFrame)::Vector{Tuple{Float64, Float64}}

Calculate mothership waypoints between sequential clusters.
Based on calc: midpoint of current and next cluster.

# Arguments
- `centroid_sequence` :  centroid lat long coordinates in sequence; including depot as the first and last cluster.

# Returns
Route as vector of lat long tuples.
"""
function get_waypoints(centroid_sequence::DataFrame)::Vector{Tuple{Float64, Float64}}
    cluster_seq_ids = centroid_sequence.cluster_id
    n_cluster_seqs = length(cluster_seq_ids)

    waypoints = Vector{Tuple{Float64, Float64}}(undef, n_cluster_seqs)

    # Add the depot as the first waypoint
    waypoints[1] = (centroid_sequence.lat[1], centroid_sequence.lon[1])

    for i in 1:(n_cluster_seqs - 1)
        current_centroid = first(centroid_sequence[centroid_sequence.cluster_id .== cluster_seq_ids[i], :])
        next_centroid = first(centroid_sequence[centroid_sequence.cluster_id .== cluster_seq_ids[i+1], :])

        # Calculate the midpoint between the current and next cluster centroids
        centroid = (
            (current_centroid.lat + next_centroid.lat) / 2,
            (current_centroid.lon + next_centroid.lon) / 2
        )
        waypoints[i+1] = centroid
    end

    # Add the depot as the last waypoint
    waypoints[n_cluster_seqs] = (centroid_sequence.lat[n_cluster_seqs], centroid_sequence.lon[n_cluster_seqs])

    return waypoints
end
