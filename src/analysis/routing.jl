
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
    get_feasible_matrix(waypoints::Vector{Point{2, Float32}}, exclusions::DataFrame)

Create a distance matrix between waypoints accounting for environmental constraints.

# Arguments
- `waypoints` : Vector of lat long tuples.
- `exclusions` : DataFrame containing exclusion zones representing given vehicle's cumulative environmental constraints.

# Returns
Feasible distance matrix between waypoints.
"""
function get_feasible_matrix(waypoints::Vector{Point{2, Float32}}, exclusions::DataFrames.DataFrame)::Matrix{Float64}
    n_waypoints = length(waypoints)-1
    feasible_matrix = zeros(Float64, n_waypoints, n_waypoints)

    for j in 1:n_waypoints
        for i in 1:n_waypoints
            if i != j
                println("($i, $j)")
                feasible_matrix[i, j] = min_feasible_dist(waypoints[i], waypoints[j], exclusions)
            end
        end
    end

    return feasible_matrix
end

"""
    min_feasible_dist(start_pt::Point{2, Float32}, end_pt::Point{2, Float32}, env_constraint::DataFrames.DataFrame)::Float64

Calculate the minimum distance between two points, avoiding exclusion zones, representing environmental constraints.

# Arguments
- `start_pt` : Lat long point coordinates.
- `end_pt` : Lat long point coordinates.
- `env_constraint` : DataFrame containing exclusion zones representing given vehicle's cumulative environmental constraints.

# Returns
Minimum feasible distance between two points.
"""
function min_feasible_dist(start_pt::Point{2, Float32}, end_pt::Point{2, Float32}, env_constraint::DataFrames.DataFrame)::Float64
    dist = 0.0
    line = LineString([Point2f(start_pt), Point2f(end_pt)])
    all_intersections = get_all_intersections(line, env_constraint)

    if (length(all_intersections) == 0)
        # No intersections
        return haversine(start_pt, end_pt)
    elseif (length(all_intersections) == 1)
        # Single intersection pair
        int_polygon = [j for i in all_intersections for j in eachrow(env_constraint).geometry if GO.distance(i[1], j) < 0.01][1]
        return min_dist_around_poly(line, int_polygon)
    else
        # Arrange intersections in order of distance from start point
        sort!(all_intersections, by = x -> haversine(start_pt, x[1]))

        # get vector of polygons in env_constraint ordered by intersection with all_intersections
        current_loc = start_pt
        int_polygons = [j for i in all_intersections for j in eachrow(env_constraint).geometry if GO.distance(i[1], j) < 0.01]
        next_loc = all_intersections[2][1]

        for poly_id in 1:length(int_polygons) - 1
            # TODO If polygon is not 'closed', this will return only return one point
            vertices = extract_vertices(int_polygons[poly_id])

            poly_dist, next_loc_idx = min_dist_around_polys(current_loc, next_loc, vertices, int_polygons[poly_id])
            dist += poly_dist

            current_loc = vertices[next_loc_idx] # all_intersections[poly_id + 1][1] # vertices[temp_loc]
            if poly_id != length(int_polygons) - 1
                # TODO: Correct this next_loc to account for mismatched polygon and intersection order
                next_loc = all_intersections[poly_id+2][1]
            end
        end

        vertices = extract_vertices(int_polygons[end])
        poly_dist, next_loc_idx = min_dist_around_polys(current_loc, end_pt, vertices, int_polygons[end])
        dist += poly_dist

        dist += haversine(vertices[next_loc_idx], end_pt)

        return dist
    end
end

function get_all_intersections(line::LineString{2, Float32}, env_constraint::DataFrame)
    all_intersections = []

    for row in eachrow(env_constraint)
        polygon = row.geometry
        intersections = GO.intersection_points(line, polygon)

        if !isempty(intersections)
            push!(all_intersections, [Point{2, Float32}(intersection[1],intersection[2]) for intersection in intersections])
        end
    end

    return all_intersections
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
    # TODO: Pre-allocate vector size
    vertices = Any[]#Vector{Any}(undef, length(polygon.exterior))
    for line in polygon.exterior
        push!(vertices, Point2f(line[1]))#vertices[i] = Point2f(line[1])
    end
    return vertices
end


"""
    min_dist_around_polys(current_loc::Point{2, Float32}, next_loc::Point{2, Float32}, vertices::Vector{Any}, polygon::Polygon)

Find shortest detour path from the start and end points of a line around an exclusion polygon.

# Arguments
- `current_loc` : Current location.
- `next_loc` : Next location.
- `vertices` : Vector of all vertices.
- `polygon` : Polygon to find paths around.

# Returns
- Min dist around polygon.
- Last index of path.
"""
function min_dist_around_polys(current_loc::Point{2, Float32}, next_loc::Point{2, Float32}, vertices::Vector{Any}, polygon::Polygon)
    side_idx = Vector{Int}(undef,2)
    # find which side line intersects with polygon
    int_pts = GO.intersection_points(LineString([Point2f(current_loc), Point2f(next_loc)]), polygon)

    # find which int_pt is closest to the start of the line
    start_side_idx = argmin([haversine(current_loc, int_pt) for int_pt in int_pts])

    # min dist threshold to replace intersection for rounding errors
    side_idx[1] = [i for i in 1:length(polygon.exterior) if GO.distance(int_pts[start_side_idx], polygon.exterior[i]) < 0.01][1]
    side_idx[2] = [i for i in 1:length(polygon.exterior) if GO.distance(int_pts[3-start_side_idx], polygon.exterior[i]) < 0.01][1]

    vert_start_idx = findfirst(x -> x == polygon.exterior[side_idx[1]][2], vertices)
    vert_end_idx = findlast(x -> x == polygon.exterior[side_idx[2]][1], vertices)

    path_a = trav_vert_path(vert_start_idx, vert_end_idx, vertices)
    path_b = append!([i for i in path_a[end]:length(vertices) if i ∉ path_a], [i for i in 1:path_a[1] if i ∉ path_a])

    # Check if first point in path is the closest point to line start
    if argmin([haversine(current_loc, vertices[i]) for i in path_a]) != 1
        path_a = reverse(path_a)
    end
    if argmin([haversine(current_loc, vertices[i]) for i in path_b]) != 1
        path_b = reverse(path_b)
    end

    paths = [path_a, path_b]
    min_dist = minimum([dist_traverse_path(current_loc, [vertices[i] for i in path]) for path in paths])
    last_idx = [stop for path in paths for stop in path if dist_traverse_path(current_loc, [vertices[i] for i in path]) == min_dist][end]

    return min_dist, last_idx
end

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

    return min(dist_traverse_path(line, [vertices[i] for i in path_a]), dist_traverse_path(line, [vertices[i] for i in path_b]))
end

"""
    trav_vert_path(vert_start_idx::Int64, vert_end_idx::Int64, vertices::Vector{Any})

Traverse a path of vartices between two vertices.

# Arguments
- `vert_start_idx` : Starting vertex index.
- `vert_end_idx` : Ending vertex index.
- `vertices` : Vector of all vertices.

# Returns
Vector sequence of vertices representing path between two vertices.
"""
function trav_vert_path(vert_start_idx::Int64, vert_end_idx::Int64, vertices::Vector{Any})
    if vert_start_idx == vert_end_idx
        return [vert_start_idx]
    elseif vert_start_idx < vert_end_idx
        return [i for i in vert_start_idx:vert_end_idx]
    else
        return vcat(vert_start_idx:length(vertices), 1:vert_end_idx)
    end
end

"""
    dist_traverse_path(current_loc::Point{2, Float32}, trav_path::Vector{Vector{Point{2, Float32}}})

Calculate the distance of traversing a path from current location to the final vertex on the path.

# Arguments
- `current_loc` : Current location.
- `trav_path` : Sequence of vertices representing path between two vertices.

# Returns
Haversine distance of traversing the path.
"""
function dist_traverse_path(current_loc::Point{2, Float32}, trav_path::Vector{Point{2, Float32}})
    # TODO: Check first point in path is the closest to line start
    dist = haversine(current_loc, trav_path[1])

    for i in 1:length(trav_path)-1
        dist += haversine(trav_path[i], trav_path[i+1])
    end

    return dist
end

"""
    dist_traverse_path(line::LineString{2, Float32}, trav_path::Vector{Point{2, Float32}})

Calculate the full distance of traversing a path from start of line through vertices to end of line.

# Arguments
- `line` : LineString intersecting with polygon.
- `trav_path` : Sequence of vertices representing path between two vertices.

# Returns
Haversine distance of traversing the path.
"""
function dist_traverse_path(line::LineString{2, Float32}, trav_path::Vector{Point{2, Float32}})
    dist = haversine(line[1][1], trav_path[1])

    for i in 1:length(trav_path)-1
        dist += haversine(trav_path[i], trav_path[i+1])
    end
    dist += haversine(trav_path[end], line[1][2])

    return dist
end
