
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
    # Create exclusion zones based on the bathymetry data
    exclusion_zones = env_constraint .<= threshold
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

    # Start at depot
    current_location = 1
    push!(tour, current_location)
    visited[current_location] = true

    while length(tour) <= num_clusters
        # Find the nearest unvisited neighbor
        distances = dist_matrix[current_location, :]
        distances[visited] .= Inf
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

# Arguments
- `waypoints` : Vector of lat long tuples.
- `exclusions` : DataFrame containing exclusion zones representing given vehicle's cumulative environmental constraints.

# Returns
Feasible distance matrix between waypoints.
"""
function get_feasible_matrix(waypoints::Vector{Tuple{Float64, Float64}}, exclusions::DataFrames.DataFrame)::Matrix{Float64}
    n_waypoints = length(waypoints)-1
    feasible_matrix = zeros(Float64, n_waypoints, n_waypoints)

    for j in 1:n_waypoints
        for i in 1:n_waypoints
            if i != j
                feasible_matrix[i, j] = min_feasible_dist(waypoints[i], waypoints[j], exclusions)
            end
        end
    end

    return feasible_matrix
end

"""
    min_feasible_dist(start_pt::Tuple{Float64, Float64}, end_pt::Tuple{Float64, Float64}, env_constraint::DataFrames.DataFrame)::Float64

Calculate the minimum distance between two points, avoiding exclusion zones, accounting for environmental constraints.

# Arguments
- `start_pt` : Tuple of lat long coordinates.
- `end_pt` : Tuple of lat long coordinates.
- `env_constraint` : DataFrame containing exclusion zones representing given vehicle's cumulative environmental constraints.

# Returns
Minimum feasible distance between two points.
"""
function min_feasible_dist(start_pt::Tuple{Float64, Float64}, end_pt::Tuple{Float64, Float64}, env_constraint::DataFrames.DataFrame)::Float64
    dist = 0.0
    line = LineString([Point2f(start_pt), Point2f(end_pt)])

    for row in eachrow(env_constraint)
        # TODO navigate around multiple polygons
        polygon = row.geometry
        # If the line intersects with any polygon in env_constraint, find shortest path around polygon
        if GO.intersects(line, polygon)
            # TODO If polygon is not 'closed', this will return only return one point

            # TODO: if line intersects with multiple polygons, find shortest cumulative path around all!

            dist += min_dist_around_poly(line, polygon)
        else
            dist += haversine(start_pt, end_pt)
        end
    end
    return dist
end

"""
    extract_vertices(polygon::Polygon)

Extract vector of polygon's vertices.

# Arguments
- `polygon` : Polygon to extract vertices from.

# Returns
Vector of all vertices of polygon.
"""
function extract_vertices(polygon::Polygon)
    vertices = Any[]#Vector{Any}(undef, length(polygon.exterior))
    for line in polygon.exterior
        push!(vertices, Point2f(line[1]))#vertices[i] = Point2f(line[1])
    end
    return vertices
end

"""
    paths_around_poly(line, vertices, polygon::Polygon)

Find the two paths around a polygon.

# Arguments
- `line` : LineString intersecting with polygon.
- `vertices` : Vector of polygon vertices.
- `polygon` : Polygon to find paths around.

# Returns
Two paths (anti- and clockwise) around polygon.
"""
function min_dist_around_poly(line, polygon::Polygon)
    side_ints_idx = Vector{Int}(undef,2)    #Tuple{Int, Int}[]
    int_pts = GO.intersection_points(line, polygon)

    # Find start and end indices of sides intersecting with line
    start_side_idx = argmin([haversine(line[1][1], int_pt) for int_pt in int_pts])
    side_ints_idx[1] = [i for i in 1:length(polygon.exterior) if GO.distance(int_pts[start_side_idx], polygon.exterior[i]) < 0.01][1]
    side_ints_idx[2] = [i for i in 1:length(polygon.exterior) if GO.distance(int_pts[3-start_side_idx], polygon.exterior[i]) < 0.01][1]

    vertices = extract_vertices(polygon)
    vert_start_idx = findfirst(x -> x == polygon.exterior[side_ints_idx[1]][2], vertices)
    vert_end_idx = findlast(x -> x == polygon.exterior[side_ints_idx[2]][1], vertices)

    path_a = trav_vert_path(vert_start_idx, vert_end_idx, vertices)
    path_b = append!([i for i in path_a[end]:length(vertices) if i ∉ path_a], [i for i in 1:path_a[1] if i ∉ path_a])

    # Reverse if first point in path is not closest point to start
    if argmin([haversine(line[1][1], vertices[i]) for i in path_a]) != 1
        path_a = reverse(path_a)
    end
    if argmin([haversine(line[1][1], vertices[i]) for i in path_b]) != 1
        path_b = reverse(path_b)
    end

    return min(dist_traverse_path(line, [vertices[i] for i in path_anti]), dist_traverse_path(line, [vertices[i] for i in path_clock]))
end

"""
    trav_vert_path(vert_start::Point{2, Float32}, vert_end::Point{2, Float32}, vertices::Vector{Point{2, Float32}}, direction_flag::Bool)

Traverse a path between two vertices.

# Arguments
- `vert_start` : Starting vertex.
- `vert_end` : Ending vertex.
- `vertices` : Vector of all vertices.
- `direction_flag` : Direction of traversal. 1 = anti-clockwise, 0 = clockwise.

# Returns
Sequence of vertices representing path between two vertices.
"""
function trav_vert_path(vert_start_idx::Point{2, Float32}, vert_end_idx::Point{2, Float32}, vertices::Vector{Point{2, Float32}})
    if vert_start_idx == vert_end_idx
        return [vert_start]
    else
        return [vertices[i] for i in vert_start_idx:vert_end_idx]
    end
end

"""
    dist_traverse_path(line, trav_path)

Calculate the distance of traversing a path.

# Arguments
- `line` : LineString intersecting with polygon.
- `trav_path` : Sequence of vertices representing path between two vertices.

# Returns
Haversine distance of traversing the path.
"""
function dist_traverse_path(line, trav_path)
    dist = haversine(line[1][1], trav_path[1])

    for i in 1:length(trav_path)-1
        dist += haversine(trav_path[i], trav_path[i+1])
    end
    dist += haversine(trav_path[end], line[1][2])

    return dist
end
