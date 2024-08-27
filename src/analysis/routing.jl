
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
Route as vector of lat long points. Depot is first and last waypoints.
"""
function get_waypoints(centroid_sequence::DataFrame)::Vector{Point{2, Float32}}
    cluster_seq_ids = centroid_sequence.cluster_id
    n_cluster_seqs = length(cluster_seq_ids)

    waypoints = Vector{Tuple{Float64, Float64}}(undef, n_cluster_seqs)

    waypoints[1] = (centroid_sequence.lat[1], centroid_sequence.lon[1])

    for i in 1:(n_cluster_seqs - 1)
        current_centroid = first(centroid_sequence[centroid_sequence.cluster_id .== cluster_seq_ids[i], :])
        next_centroid = first(centroid_sequence[centroid_sequence.cluster_id .== cluster_seq_ids[i+1], :])

        centroid = (
            (current_centroid.lat + next_centroid.lat) / 2,
            (current_centroid.lon + next_centroid.lon) / 2
        )
        waypoints[i+1] = centroid
    end

    waypoints[n_cluster_seqs] = (centroid_sequence.lat[n_cluster_seqs], centroid_sequence.lon[n_cluster_seqs])

    return [GeometryBasics.Point{2, Float32}(wp...) for wp in waypoints]
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
function get_feasible_matrix(waypoints::Vector{Point{2, Float32}}, exclusions::DataFrame)::Matrix{Float64}
    n_waypoints = length(waypoints) - 1
    return [i != j ? shortest_feasible_path((waypoints[i], waypoints[j]), exclusions)[1] : 0.0 for j in 1:n_waypoints, i in 1:n_waypoints]
end

function shortest_feasible_path(line_pts::Tuple{Point{2, Float32}, Point{2, Float32}}, exclusions::DataFrame)
    pts = extract_unique_vertices(exclusions)
    insert!(pts, 1, line_pts[1])
    push!(pts, line_pts[2])

    g = build_graph(pts, exclusions)

    path = a_star(g, 1, length(pts), weights(g))
    dist = sum([g.weights[path[i].src, path[i].dst] for i in 1:length(path)])

    return dist, path
end

function extract_unique_vertices(exclusions::DataFrame)::Vector{Point{2, Float32}}
    unique_vertices = Set{Point{2, Float32}}()

    for row in eachrow(exclusions)
        polygon = row.geometry

        exterior_ring = AG.getgeom(polygon, 0)
        n_pts = AG.ngeom(exterior_ring)

        for i in 0:n_pts - 1
            x, y, _ = AG.getpoint(exterior_ring, i)
            point = Point{2, Float32}(x, y)

            if !(point in unique_vertices)
                push!(unique_vertices, point)
            end
        end
    end

    return collect(unique_vertices)  # Convert the Set back to a Vector
end

function build_graph(pts::Vector{Point{2, Float32}}, exclusions::DataFrame)::SimpleWeightedGraph{Int64, Float64}
    g = SimpleWeightedGraph(length(pts))

    for j in 1:length(pts)
        for i in 1:j-1
            line = LineString([pts[i], pts[j]])

            # TODO Eliminate pairs of points that do not intersect with any polygon using bounding box checks
            # if !can_intersect_bounding_box(pts[i], pts[j], exclusions)

            if !intersects_polygon(line, exclusions)
                dist = haversine(pts[i], pts[j])
                add_edge!(g, i, j, dist)
                # println(i, '\t', j, '\t', dist)
            end
        end
    end

    return g
end

function intersects_polygon(line::LineString{2, Float32}, exclusions::DataFrame)::Bool
    for row in eachrow(exclusions)
        polygon = row.geometry
        intersections = [Point{2, Float32}(float(x), float(y)) for (x, y) in GO.intersection_points(line, polygon)]

        is_tangent, only_vertex_intersection = check_tangency_and_vertex_intersection(line, polygon, intersections)

        if !isempty(intersections) && !is_tangent && !only_vertex_intersection
            return true
        end
    end
    return false
end

function check_tangency_and_vertex_intersection(line::LineString{2, Float32}, polygon::AG.IGeometry, intersections::Vector{Point{2, Float32}})::Tuple{Bool, Bool}
    exterior_ring = AG.getgeom(polygon, 0)  # 0 indicates the exterior ring
    n_points = AG.ngeom(exterior_ring)

    verts = Set{Point{2, Float32}}()

    is_tangent = false
    only_vertex_intersection = true

    for i in 0:n_points - 2
        x1, y1, _ = AG.getpoint(exterior_ring, i)
        x2, y2, _ = AG.getpoint(exterior_ring, i+1)
        edge = LineString(Point{2, Float32}[(x1, y1), (x2, y2)])

        if length(GO.intersection_points(line, edge)) > 1
            is_tangent = true
        end

        push!(verts, Point{2, Float32}(x1, y1))
    end

    x, y, _ = AG.getpoint(exterior_ring, n_points - 1)
    push!(verts, Point{2, Float32}(x, y))

    for intersection in intersections
        if !(intersection in verts)
            only_vertex_intersection = false
        end
    end

    return is_tangent, only_vertex_intersection
end
