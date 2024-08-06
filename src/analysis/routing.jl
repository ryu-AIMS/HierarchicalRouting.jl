
function create_exclusion_zones(target_bathy::Raster, ms_depth)
    # Create exclusion zones based on the bathymetry data
    exclusion_zones = target_bathy .<= ms_depth
    return exclusion_zones
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

"""
    get_feasible_matrix(waypoints::Vector{Tuple{Float64, Float64}}, exclusions::DataFrame)

Create a distance matrix between waypoints accounting for environmental constraints.
"""
function get_feasible_matrix(waypoints::Vector{Tuple{Float64, Float64}}, exclusions::DataFrames.DataFrame)::Matrix{Float64}
    n_waypoints = length(waypoints)-1
    feasible_matrix = zeros(Float64, n_waypoints, n_waypoints)

    for j in 1:n_waypoints
        for i in 1:n_waypoints
            if i != j
                println(i,j)
                feasible_matrix[i, j] = min_feasible_dist(waypoints[i], waypoints[j], exclusions)   # haversine(waypoints[i], waypoints[j])
            end
        end
    end

    return feasible_matrix
end

function min_feasible_dist(start_pt::Tuple{Float64, Float64}, end_pt::Tuple{Float64, Float64}, env_constraint::DataFrames.DataFrame)::Float64
    line = LineString([Point2f(start_pt), Point2f(end_pt)])

    # If the line intersects with any polygon in env_constraint, find shortest path around polygon
    for row in eachrow(env_constraint)
        polygon = row.geometry
        if GO.intersects(line, polygon)
            int1, int2 = GO.intersection_points(line, polygon)       # If polygon is not 'closed', this will return only return one point

            vertices = extract_vertices(polygon)

            path_trav_anti, path_trav_clock = paths_around_poly(line, vertices, polygon)
            println(line)
            min(dist_traverse_path(line, path_trav_anti), dist_traverse_path(line, path_trav_clock))

            return min(dist_clock, dist_anti)  # Use Inf or a large value to indicate infeasibility
        end
    end
    return haversine(start_pt, end_pt)
end

function extract_vertices(polygon::Polygon)
    vertices = Any[]#Vector{Any}(undef, length(polygon.exterior))
    for line in polygon.exterior
        push!(vertices, Point2f(line[1]))#vertices[i] = Point2f(line[1])
    end
    return vertices
end

"""
    path_around_poly(side_idx::Int, polygon::Polygon)

Find the two paths around a polygon.

"""
function paths_around_poly(line, vertices, polygon::Polygon)
    side_idx = Int[] #Tuple{LineString, LineString}[]
    # find which side line intersects with polygon
    for s in 1:length(polygon.exterior)    # for side in polygon.exterior
        if GO.intersects(line, polygon.exterior[s])
            push!(side_idx, s)
        end
    end

    # Find vertex at end of side in side_idx
    side_enter = polygon.exterior[minimum(side_idx)]
    side_exit = polygon.exterior[maximum(side_idx)]

    vert_start_anti = side_enter[2]
    vert_start_clock = side_enter[1]

    vert_end_anti = side_exit[1]
    vert_end_clock = side_exit[2]

    path_anti = trav_vert_path(vert_start_anti, vert_end_anti, vertices, true)
    path_clock = trav_vert_path(vert_start_clock, vert_end_clock, vertices, false)

    return path_anti, path_clock
end

function trav_vert_path(vert_start::Point{2, Float32}, vert_end::Point{2, Float32}, vertices::Vector{Point{2, Float32}}, direction_flag::Bool)
        # add all vertices between start_vertex_idx and end_vertex_idx to side_int
        vert_start_idx = findfirst(x -> x == vert_start, vertices)
        vert_end_idx = findfirst(x -> x == vert_end, vertices)

        if direction_flag
            path = [vertices[i] for i in vert_start_idx:vert_end_idx]
        else
            path = [vertices[i] for i in vert_start_idx:-1:vert_end_idx]
        end
    return path
end

function dist_traverse_path(line, trav_path)
    dist = haversine(line[1][1], trav_path[1])

    for i in 1:length(trav_path)-1
        dist += haversine(trav_path[i], trav_path[i+1])
    end
    dist += haversine(trav_path[end], line[1][2])

    return dist
end
