
@kwdef struct MothershipSolution
    cluster_sequence::DataFrame
    route::DataFrame
    cost::Float64
    line_strings::Vector{LineString{2, Float64}} = []
end

struct Sortie
    nodes::Vector{Point{2, Float64}}
    cost::Float64
end

struct TenderSolution
    id::Int
    sorties::Vector{Sortie}
    # cost::Float64   # TODO: delete? redundant with sum(...) and critical path metric...
    start::Point{2, Float64}
    finish::Point{2, Float64}
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
    exclusion_zones = (env_constraint .<= threshold) # .| (env_constraint .== env_constraint.missingval)
    return exclusion_zones
end

"""
    nearest_neighbour(nodes::DataFrame, exclusions::DataFrame)

Apply the nearest neighbor algorithm starting from the depot (1st row/col) and returning to the depot.

# Arguments
- `nodes` : DataFrame containing id, lat, lon. Depot is cluster 0 in row 1.
- `exclusions` : DataFrame containing exclusion zones.
# Returns
- `solution` : MSRoutingSolution object containing:
    - `cluster_sequence` : Centroid sequence DataFrame (containing id, lat, lon).
    - `route` : DataFrame of waypoints on the route.
    - `distance` : Total distance of the route.
- `dist_matrix` : Distance matrix between centroids.
- `feasible_path` : A vector of tuples containing the graph, point to index mapping, and edges for each pair of nodes.
"""
function nearest_neighbour(nodes::DataFrame, exclusions::DataFrame)
    dist_matrix, feasible_path = get_feasible_matrix([Point{2, Float64}(row.lon, row.lat) for row in eachrow(nodes)], exclusions)

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

    ordered_nodes = nodes[[findfirst(==(id), nodes.id) for id in cluster_sequence], :]
    waypoints = get_waypoints(ordered_nodes)

    _, waypoint_feasible_path = get_feasible_matrix(waypoints.waypoint, exclusions)
    paths = get_linestrings(waypoint_feasible_path, waypoints.waypoint)

    return MothershipSolution(cluster_sequence=ordered_nodes, route=waypoints, cost=total_distance, line_strings=paths), dist_matrix, feasible_path
end

"""
    get_waypoints(sequence::DataFrame)::DataFrame

Calculate mothership waypoints between sequential clusters.
For each cluster, waypoint 1/3 dist before and after cluster centroid.

# Arguments
- `sequence` : id, amd centroid lat, long coordinates in sequence; including depot as the first and last rows wwith id=0.

# Returns
- `waypoint_df` : DataFrame for each waypoint on route.
                    Cols: waypoint::Point{2, Float64}, connecting_clusters::NTuple{2, Int64} reference to previous and next id.
                    Depot is included as first and last rows.
"""
function get_waypoints(sequence::DataFrame)::DataFrame
    # TODO: Implement convex hull exclusion zones
    # TODO: Use graph to determine feasible paths and waypoints along path
    n_cluster_seqs = nrow(sequence)

    waypoints = Vector{Point{2, Float64}}(undef, 2*(n_cluster_seqs-2)+2)
    connecting_clusters = Vector{NTuple{2, Int64}}(undef, 2*(n_cluster_seqs-2)+2)

    waypoints[1] = (sequence.lon[1], sequence.lat[1])
    connecting_clusters[1] = (sequence.id[1], sequence.id[1])

    for i in 2:(n_cluster_seqs - 1)
        prev_lat, prev_lon, prev_clust = sequence.lat[i-1], sequence.lon[i-1], sequence.id[i-1]
        current_lat, current_lon, current_clust = sequence.lat[i], sequence.lon[i], sequence.id[i] #first(sequence[sequence.id .== cluster_seq_ids[i], :])
        next_lat, next_lon, next_clust = sequence.lat[i+1], sequence.lon[i+1], sequence.id[i+1] #first(sequence[sequence.id .== cluster_seq_ids[i+1], :])

        prev_waypoint_lat = (2/3 * (current_lat) + 1/3 * (prev_lat))
        prev_waypoint_lon = (2/3 * (current_lon) + 1/3 * (prev_lon))

        next_waypoint_lat = (2/3 * (current_lat) + 1/3 * (next_lat))
        next_waypoint_lon = (2/3 * (current_lon) + 1/3 * (next_lon))

        waypoints[2*i-2] = (prev_waypoint_lon, prev_waypoint_lat)
        connecting_clusters[2*i-2] = (prev_clust, current_clust)

        waypoints[2*i-1] = (next_waypoint_lon, next_waypoint_lat)
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
    two_opt(ms_soln_current::MothershipSolution, dist_matrix::Matrix{Float64}, feasible_path, exclusions::DataFrame)

Apply the 2-opt heuristic to improve the current MothershipSolution route (by uncrossing crossed links) between waypoints.

# Arguments
- `ms_soln_current` : Current MothershipSolution - from nearest_neighbour.
- `dist_matrix` : Distance matrix between waypoints. Depot is the first, but not last point.
- `feasible_path` : A vector of tuples containing the graph, point to index mapping, and edges for each pair of nodes.
- `exclusions` : DataFrame containing exclusion zones.

# Returns
- `solution` : MothershipSolution object containing:
    - `cluster_sequence` : DataFrame of cluster centroids in sequence (containing id, lon, lat).
    - `route` : DataFrame of waypoints characterising route.
    - `cost` : Total distance of the route.
    - `line_strings` : Vector of LineString objects for each path.
"""
function two_opt(ms_soln_current::MothershipSolution, dist_matrix::Matrix{Float64}, feasible_path, exclusions::DataFrame)

    nodes = ms_soln_current.cluster_sequence

    # If depot is last row, remove
    if nodes.id[1] == nodes.id[end]
        nodes = nodes[1:end-1, :]
    end

    # Initialize route as ordered waypoints
    best_route = [row.id+1 for row in eachrow(nodes)]
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

    # Re-orient route to start from the depot (1), and add the depot as final point
    best_route = orient_route(best_route)
    push!(best_route, best_route[1])

    # Adjust sequence to zero-based indexing where depot = 0
    best_route .-= 1

    ordered_nodes = nodes[[findfirst(==(id), nodes.id) for id in best_route], :]

    waypoints = get_waypoints(ordered_nodes)
    _, waypoint_feasible_path = get_feasible_matrix(waypoints.waypoint, exclusions)

    paths = get_linestrings(waypoint_feasible_path, waypoints.waypoint)

    return MothershipSolution(cluster_sequence=ordered_nodes, route=waypoints, cost=best_distance, line_strings=paths)
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
    tender_sequential_nearest_neighbour(cluster::Cluster, waypoints::NTuple{2, Point{2, Float64}}, n_tenders::Int, t_cap::Int, exclusions::DataFrame)

Assign nodes to tenders sequentially based on nearest neighbor.

# Arguments
- `cluster` : Cluster object containing nodes.
- `waypoints` : Tuple of start and end waypoints.
- `n_tenders` : Number of tenders.
- `t_cap` : Tender capacity.
- `exclusions` : DataFrame containing exclusion zones.

# Returns
- `solution` : TenderRoutingSolution object containing:
    - `cluster_id` : Cluster ID.
    - `sorties` : Vector of Sortie objects containing nodes and sortie distance.
    - `cost` : Total distance of the route.
    - `start` : Start waypoint.
    - `finish` : End waypoint.
"""
function tender_sequential_nearest_neighbour(cluster::Cluster, waypoints::NTuple{2, Point{2, Float64}}, n_tenders::Int, t_cap::Int, exclusions::DataFrame)
    nodes = [waypoints[1]]
    append!(nodes, cluster.nodes)

    dist_matrix = get_feasible_matrix([Point{2, Float64}(node[2], node[1]) for node in nodes], exclusions)[1]

    tender_tours = [fill(0, t_cap) for _ in 1:n_tenders]

    visited = falses(length(nodes))
    visited[1] = true

    # for each tender in number of tenders, sequentially assign closest nodes tender-by-tender stop-by-stop
    for i in 1:t_cap
        for j in 1:n_tenders
            # TODO: check conditions
            if count(visited) == length(visited)
                break
            end

            current_node = i == 1 ? 1 : tender_tours[j][i-1]

            distances = dist_matrix[current_node, :]
            vis_idxs = findall(v -> v != 0, visited)
            distances[vis_idxs] .= Inf
            nearest_idx = argmin(distances)

            tender_tours[j][i] = nearest_idx
            visited[nearest_idx] = true
        end
    end

    # delete excess elements and remove empty tours
    tender_tours = [count(tour .== 0) > 0 ? tour[1:findfirst(==(0), tour)-1] : tour for tour in tender_tours if any(!=(0), tour)]

    # TODO: Handling for empty tender tours
    sortie_dist = sortie_cost(tender_tours, dist_matrix)

    total_distance = sum(sortie_dist)

    return TenderSolution(cluster.id, [Sortie([nodes[stop] for stop in tender_tours[t]], sortie_dist[t]) for t in 1:length(tender_tours)], waypoints[1], waypoints[2]), dist_matrix
end
