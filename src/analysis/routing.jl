
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

    return [i != j ? shortest_feasible_path(LineString([waypoints[i], waypoints[j]]), exclusions)[1] : 0.0 for j in 1:n_waypoints, i in 1:n_waypoints]
end

function shortest_feasible_path(line::LineString{2, Float32, Point{2, Float32}}, exclusions::DataFrame)
    points = Vector{Point{2, Float32}}(undef, sum([length(row.geometry.exterior.points) for row in eachrow(exclusions)]) + 2)

    points[2:end-1] = extract_unique_vertices(exclusions)
    points[1] = line.points[1][1]
    points[end] = line.points[1][2]

    g = build_graph(points, exclusions)

    # TODO: Convert to haversine distances
    path = a_star(g, 1, length(points), weights(g))
    dist = sum([g.weights[path[i].src, path[i].dst] for i in 1:length(path)])
    return dist, path
end

function extract_unique_vertices(exclusions::DataFrame)::Vector{Point{2, Float32}}
    vertices = Set{Point{2, Float32}}()

    for row in eachrow(exclusions)

        for line in row.geometry.exterior
            for point in line
                if !(point in vertices)
                    push!(vertices, point)
                end
            end
        end
    end

    return collect(vertices)  # Convert the Set back to a Vector
end

function build_graph(vertices::Vector{Point{2, Float32}}, exclusions::DataFrame)::SimpleWeightedGraph{Int64, Float64}
    g = SimpleWeightedGraph(length(vertices))

    for j in 1:length(vertices)
        for i in 1:j-1
            if !intersects_polygon(LineString([vertices[i], vertices[j]]), exclusions)
                add_edge!(g, i, j, haversine(vertices[i], vertices[j]))
            end
        end
    end

    return g
end

function intersects_polygon(line::LineString{2, Float32}, exclusions::DataFrame)::Bool

    for row in eachrow(exclusions)
        polygon = row.geometry
        ints = GO.intersection_points(line, polygon)

        if !isempty(ints) && !is_tangent(line, polygon) && !only_vertex_int([Point{2, Float32}(Float32(x), Float32(y)) for (x, y) in ints], polygon)
            return true
        end
    end
    return false
end

function is_tangent(line::LineString, polygon::Polygon)
    for edge in polygon.exterior
        if length(GO.intersection_points(line, edge)) > 1
            return true
        end
   end
    return false
end

function only_vertex_int(ints::Vector{Point{2, Float32}}, polygon::Polygon)
    for int in ints
        if !any([int in pt for pt in [line for line in polygon.exterior]])
            return false
        end
    end
    return true
end
