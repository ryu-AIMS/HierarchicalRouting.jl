
struct Route
    nodes::Vector{Point{2, Float64}}
    dist_matrix::Matrix{Float64}
    line_strings::Vector{LineString{2, Float64}}
end

@kwdef struct MothershipSolution
    cluster_sequence::DataFrame
    route::Route
end

struct TenderSolution
    id::Int
    start::Point{2, Float64}
    finish::Point{2, Float64}
    sorties::Vector{Route}
    dist_matrix::Matrix{Float64}
end

struct MSTSolution
    clusters::Vector{Cluster}
    mothership::MothershipSolution
    tenders::Vector{TenderSolution}
end

"""
    create_exclusion_zones(env_constraint::Raster, threshold::Float64)

Create exclusion zones based on environmental raster data and vessel threshold.

# Arguments
- `env_constraint` : Environmental constraint raster.
- `threshold` : Threshold for given vessel's environmental constraint.

# Returns
Exclusion zones for environmental constraint and vessel threshold provided.
"""
function create_exclusion_zones(env_constraint::Raster, threshold::Float64)
    # TODO: kwargs to specify max/min threshold values
    return (env_constraint .!= env_constraint.missingval) .&& (env_constraint .>= threshold)
end

"""
    adjust_waypoint(
        waypoint::Point{2, Float64},
        exclusions::DataFrame,
    )::Point{2, Float64}

Adjust waypoint if inside exclusion zone to closest boundary point outside exclusion zone.

# Arguments
- `waypoint` : Point{2, Float64} waypoint.
- `exclusions` : DataFrame containing exclusion zones.

# Returns
Adjusted waypoint if inside exclusion zone, else original waypoint.
"""
function adjust_waypoint(
    waypoint::Point{2, Float64},
    exclusions::DataFrame,
)::Point{2, Float64}

    waypoint_geom = AG.createpoint(waypoint[1], waypoint[2])

    containing_polygons = [
        AG.convexhull(polygon) for polygon in exclusions.geometry
            if AG.contains(AG.convexhull(polygon), waypoint_geom)
    ]

    if isempty(containing_polygons)
        return waypoint
    end

    # TODO Consider cases for multiple exclusion zones
    union_poly = reduce(AG.union, containing_polygons)

    exterior_ring = AG.getgeom(union_poly, 0)
    n_points = AG.ngeom(exterior_ring)

    boundary_points = [
        Point(AG.getpoint(exterior_ring, i)[1:2]...) for i in 0:n_points - 1
    ]

    valid_boundary_points = [
        pt for pt in boundary_points
            if !point_in_exclusion(pt, union_poly)
    ]

    closest_point = argmin(
        p -> sqrt((p[1] - waypoint[1])^2 + (p[2] - waypoint[2])^2),
        valid_boundary_points
    )

    # Recursively call adjust_waypoint() until waypoint is outside exclusion zone
    if point_in_exclusion(closest_point, exclusions)
        return adjust_waypoint(closest_point, exclusions)
    end

    return closest_point
end

"""
    get_waypoints(sequence::DataFrame, exclusions::DataFrame)::DataFrame

Calculate mothership waypoints between sequential clusters.
For each cluster, waypoint 1/3 dist before and after cluster centroid,
unless within exclusion zone, then adjust to closest boundary point.

# Arguments
- `sequence` : id, amd centroid lat, long coordinates in sequence; including depot as the first and last rows wwith id=0.
- `exclusions` : DataFrame containing exclusion zones.

# Returns
- `waypoint_df` : DataFrame for each waypoint on route. Depot included as first and last rows.
                    Cols: waypoint::Point{2, Float64}, connecting_clusters::NTuple{2, Int64} reference to previous and next id.
"""
function get_waypoints(sequence::DataFrame, exclusions::DataFrame)::DataFrame
    # TODO: Implement convex hull exclusion zones
    # TODO: Use graph to determine feasible paths and waypoints along path
    n_cluster_seqs = nrow(sequence)

    waypoints = Vector{Point{2, Float64}}(undef, 2*(n_cluster_seqs-2)+2)
    connecting_clusters = Vector{NTuple{2, Int64}}(undef, 2*(n_cluster_seqs-2)+2)

    waypoints[1] = Point{2, Float64}(sequence.lon[1], sequence.lat[1])
    connecting_clusters[1] = (sequence.id[1], sequence.id[1])

    for i in 2:(n_cluster_seqs - 1)
        prev_lon, prev_lat, prev_clust = sequence.lon[i-1], sequence.lat[i-1], sequence.id[i-1]
        current_lon, current_lat, current_clust = sequence.lon[i], sequence.lat[i], sequence.id[i]
        next_lon, next_lat, next_clust = sequence.lon[i+1], sequence.lat[i+1], sequence.id[i+1]

        # Compute waypoints before and after the current cluster centroid
        prev_waypoint = Point{2, Float64}(2/3 * current_lon + 1/3 * prev_lon, 2/3 * current_lat + 1/3 * prev_lat)
        next_waypoint = Point{2, Float64}(2/3 * current_lon + 1/3 * next_lon, 2/3 * current_lat + 1/3 * next_lat)

        # Adjust waypoints if they are inside exclusion polygons
        prev_waypoint = point_in_exclusion(prev_waypoint, exclusions) ?
            adjust_waypoint(prev_waypoint, exclusions) : prev_waypoint
        next_waypoint = point_in_exclusion(next_waypoint, exclusions) ?
            adjust_waypoint(next_waypoint, exclusions) : next_waypoint

        waypoints[2*i-2] = prev_waypoint
        connecting_clusters[2*i-2] = (prev_clust, current_clust)

        waypoints[2*i-1] = next_waypoint
        connecting_clusters[2*i-1] = (current_clust, next_clust)
    end

    waypoints[2*(n_cluster_seqs-2)+2] = (sequence.lon[end], sequence.lat[end])
    connecting_clusters[2*(n_cluster_seqs-2)+2] = (sequence.id[end], sequence.id[end])

    waypoint_df = DataFrame(
        waypoint = waypoints,
        connecting_clusters = connecting_clusters
    )
    return waypoint_df
end

"""
    nearest_neighbour(
        nodes::DataFrame,
        exclusions_mothership::DataFrame,
        exclusions_tender::DataFrame,
    )

Apply the nearest neighbor algorithm starting from the depot (1st row/col) and returning to the depot.

# Arguments
- `nodes` : DataFrame containing id, lat, lon. Depot is cluster 0 in row 1.
- `exclusions_mothership` : DataFrame containing exclusion zones for mothership.
- `exclusions_tender` : DataFrame containing exclusion zones for tenders.

# Returns
- `solution` : MothershipSolution object containing:
    - `cluster_sequence` : Centroid sequence DataFrame (containing id, lat, lon).
    - `route` : Route object containing waypoints, distance matrix, and line strings.
        - `nodes` : Vector of Point{2, Float64} waypoints.
        - `dist_matrix` : Distance matrix between centroids.
        - `line_strings` : Vector of LineString objects for each path.
"""
function nearest_neighbour(
    cluster_centroids::DataFrame,
    exclusions_mothership::DataFrame,
    exclusions_tender::DataFrame,
)
    # TODO: Use vector rather than DataFrame for cluster_centroids
    # adjust_waypoints to ensure not within exclusion zones - allows for feasible path calc
    feasible_centroids = HierarchicalRouting.adjust_waypoint.(
        Point{2,Float64}.(cluster_centroids.lon, cluster_centroids.lat),
        [exclusions_mothership]
    )

    cluster_centroids[!, :lon] = [pt[1] for pt in feasible_centroids]
    cluster_centroids[!, :lat] = [pt[2] for pt in feasible_centroids]

    # Create distance matrix between feasible nodes - cluster centroids
    dist_matrix = get_feasible_matrix(
        feasible_centroids,
        exclusions_mothership
    )[1]

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

    # Return to the depot and adjust cluster_sequence to zero-based indexing
    push!(tour, 1)
    total_distance += dist_matrix[current_location, 1]
    cluster_sequence = tour .- 1

    ordered_centroids = cluster_centroids[[findfirst(==(id), cluster_centroids.id) for id in cluster_sequence], :]

    # combine exclusions for mothership and tenders
    exclusions_all = vcat(exclusions_mothership, exclusions_tender)
    waypoints = get_waypoints(ordered_centroids, exclusions_all)

    # Calc feasible path between waypoints.
    waypoint_dist_vector, waypoint_path_vector = get_feasible_vector(waypoints.waypoint, exclusions_mothership)
    path = vcat(waypoint_path_vector...)

    return MothershipSolution(
        cluster_sequence=ordered_centroids,
        route=Route(waypoints.waypoint, dist_matrix, path)
    )
end

"""
    two_opt(
        ms_soln_current::MothershipSolution,
        exclusions_mothership::DataFrame,
        exclusions_tender::DataFrame,
    )

Apply the 2-opt heuristic to improve the current MothershipSolution route (by uncrossing crossed links) between waypoints.

# Arguments
- `ms_soln_current` : Current MothershipSolution - from nearest_neighbour.
- `exclusions_mothership` : DataFrame containing exclusion zones for mothership.
- `exclusions_tender` : DataFrame containing exclusion zones for tenders.

# Returns
- `solution` : MothershipSolution object containing:
    - `cluster_sequence` : DataFrame of cluster centroids in sequence (containing id, lon, lat).
    - `route` : Route object containing waypoints, distance matrix, and line strings.
        - `nodes` : Vector of Point{2, Float64} waypoints.
        - `dist_matrix` : Distance matrix between centroids.
        - `line_strings` : Vector of LineString objects for each path.
"""
function two_opt(
    ms_soln_current::MothershipSolution,
    exclusions_mothership::DataFrame,
    exclusions_tender::DataFrame,
)

    cluster_centroids = ms_soln_current.cluster_sequence
    dist_matrix = ms_soln_current.route.dist_matrix

    # If depot is last row, remove
    if cluster_centroids.id[1] == cluster_centroids.id[end]
        cluster_centroids = cluster_centroids[1:end-1, :]
    end

    # Initialize route as ordered waypoints
    best_route = [row.id+1 for row in eachrow(cluster_centroids)] # cluster_centroids = cluster_centroids[1:end-1, :]
    best_distance = return_route_distance(best_route, dist_matrix)
    improved = true

    while improved
        improved = false
        for j in 1:(length(best_route) - 1)
            for i in (j + 1):length(best_route)
                new_route = two_opt_swap(best_route, j, i)
                new_distance = return_route_distance(new_route, dist_matrix)

                if new_distance < best_distance
                    best_route = new_route
                    best_distance = new_distance
                    improved = true
                end
            end
        end
    end

    # Re-orient route to start from and end at the depot, and adjust to zero-based indexing
    best_route = orient_route(best_route)
    push!(best_route, best_route[1])
    best_route .-= 1

    ordered_nodes = cluster_centroids[[findfirst(==(id), cluster_centroids.id) for id in best_route], :]
    exclusions_all = vcat(exclusions_mothership, exclusions_tender)
    waypoints = get_waypoints(ordered_nodes, exclusions_all)

    waypoint_dist_vector, waypoint_path_vector = get_feasible_vector(
        waypoints.waypoint,
        exclusions_mothership
    )

    path = vcat(waypoint_path_vector...)

    return MothershipSolution(
        cluster_sequence=ordered_nodes,
        route=Route(waypoints.waypoint, dist_matrix, path)
    )
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

"""
    tender_sequential_nearest_neighbour(
        cluster::Cluster,
        waypoints::NTuple{2, Point{2, Float64}},
        n_tenders::Int,
        t_cap::Int,
        exclusions::DataFrame
    )

Assign nodes to tenders sequentially (stop-by-stop) based on nearest neighbor.

# Arguments
- `cluster` : Cluster object containing nodes.
- `waypoints` : Tuple of start and end waypoints.
- `n_tenders` : Number of tenders.
- `t_cap` : Tender capacity.
- `exclusions` : DataFrame containing exclusion zones.

# Returns
- `solution` : TenderRoutingSolution object containing:
    - `cluster_id` : Cluster ID.
    - `start` : Start waypoint.
    - `finish` : End waypoint.
    - `sorties` : Vector of Sortie objects containing nodes and sortie distance.
        - `nodes` : Vector of Point{2, Float64} waypoints.
        - `dist_matrix` : Distance matrix between nodes.
        - `line_strings` : Vector of LineString objects for each path.
"""
function tender_sequential_nearest_neighbour(
    cluster::Cluster,
    waypoints::NTuple{2, Point{2, Float64}},
    n_tenders::Int,
    t_cap::Int,
    exclusions::DataFrame
)
    nodes = [[waypoints[1]]; cluster.nodes]
    full_nodes = vcat(nodes, [waypoints[2]])

    # Compute the full feasible matrix and associated paths once.
    dist_matrix, path_matrix = get_feasible_matrix(full_nodes, exclusions)

    tender_tours = [Int[] for _ in 1:n_tenders]
    visited = falses(length(nodes))
    visited[1] = true
    visited[2:end] .= isinf.(dist_matrix[1, 2:length(nodes)])

    # sequentially assign closest un-visited nodes stop-by-stop, tender-by-tender
    for _ in 1:t_cap
        for t in 1:n_tenders
            if all(visited)
                break
            end

            current_node = isempty(tender_tours[t]) ? 1 : last(tender_tours[t])

            distances = dist_matrix[current_node, 1:length(nodes)]
            distances[visited] .= Inf
            nearest_idx = argmin(distances)

            push!(tender_tours[t], nearest_idx)
            visited[nearest_idx] = true
        end
    end

    # delete excess elements and remove empty tours
    tender_tours = filter(!isempty, tender_tours)

    routes = Route[]
    for tour in tender_tours
        # Route indices incl waypoints: (1), tour nodes, end (last index of full_nodes)
        route_indices = vcat(1, tour, length(full_nodes))

        # Submatrix for this route
        route_matrix = dist_matrix[route_indices, route_indices]

        route_paths = [
            (i < j ?
            path_matrix[i, j] :
            path_matrix[j, i])
            for (i, j) in zip(route_indices[1:end-1], route_indices[2:end])
        ]

        push!(routes, Route(
            cluster.nodes[tour.-1],
            route_matrix,
            vcat(route_paths...)
        ))
    end

    return TenderSolution(
        cluster.id,
        waypoints[1],
        waypoints[2],
        routes,
        dist_matrix
    )
end
