
using Optim

const VALID_SECTIONS::Vector{Symbol} = [:from, :to]

struct Route
    nodes::Vector{Point{2,Float64}}
    dist_matrix::AbstractArray{Float64}
    line_strings::Vector{LineString{2,Float64}}
end

@kwdef struct MothershipSolution
    cluster_sequence::DataFrame
    route::Route
end

struct TenderSolution
    id::Int64
    start::Point{2,Float64}
    finish::Point{2,Float64}
    sorties::Vector{Route}
end
function TenderSolution(t::TenderSolution, sorties::Vector{Route})
    return TenderSolution(
        t.id,
        t.start,
        t.finish,
        sorties,
    )
end

function generate_letter_id(t::TenderSolution)
    return generate_letter_id(t.id)
end

struct MSTSolution
    cluster_sets::Vector{Vector{Cluster}}
    mothership_routes::Vector{MothershipSolution}
    tenders::Vector{Vector{TenderSolution}}
end

"""
    rebuild_sortie(
        route::Route,
        start_new::Point{2,Float64},
        finish_new::Point{2,Float64},
        exclusions::POLY_VEC,
        sortie_start_has_moved,
        sortie_end_has_moved
    )::Route

Rebuild a sortie with updated start and finish points, adjusting the start and/or finish
segments of the line strings, ensuring feasibility by avoiding exclusion zones.

# Arguments
- `route`: The **existing** route to be updated.
- `start_new`: The new start point for the sortie.
- `finish_new`: The new finish point for the sortie.
- `exclusions`: Exclusion zone polygons for the sortie.
- `sortie_start_has_moved`: Boolean indicating if the start point has moved.
- `sortie_end_has_moved`: Boolean indicating if the finish point has moved.

# Returns
- A new Route object with updated nodes, distance matrix, and line strings.
"""
function rebuild_sortie(
    route::Route,
    start_new::Point{2,Float64},
    finish_new::Point{2,Float64},
    exclusions::POLY_VEC,
    sortie_start_has_moved,
    sortie_end_has_moved,
)::Route
    segment_to_keep = route.line_strings
    updated_dist_matrix::Vector{Float64} = get_distance_vector(route.dist_matrix)

    # Keep the segments of the line strings that contain matching start/end points
    if sortie_start_has_moved
        segment_dist, segment_to_keep = update_segment(
            segment_to_keep,
            :from,
            start_new,
            route.nodes[1],
            exclusions
        )
        updated_dist_matrix[1] = segment_dist
    end

    if sortie_end_has_moved
        segment_dist, segment_to_keep = update_segment(
            segment_to_keep,
            :to,
            route.nodes[end],
            finish_new,
            exclusions
        )
        updated_dist_matrix[end] = segment_dist
    end

    return Route(route.nodes, updated_dist_matrix, segment_to_keep)
end

"""
    update_segment(
        segment::Vector{LineString{2,Float64}},
        section::Symbol,
        point_from::Point{2,Float64},
        point_to::Point{2,Float64},
        exclusions::POLY_VEC
    )::Tuple{Float64,Vector{LineString{2,Float64}}}

Update the segment of line strings based on the specified section and points.

# Arguments
- `segment`: Vector of line strings to be updated.
- `section`: Symbol indicating whether to update from the start (`:from`) or to the end
    (`:to`).
- `point_from`: The start point to update linestring from.
- `point_to`: The end point to update linestring to.
- `exclusions`: Exclusion zone polygons for tenders to avoid.

# Returns
- A tuple containing:
  - the total segment distance, and
  - a vector of updated line strings.
"""
function update_segment(
    segment::Vector{LineString{2,Float64}},
    section::Symbol,
    point_from::Point{2,Float64},
    point_to::Point{2,Float64},
    exclusions::POLY_VEC
)::Tuple{Float64,Vector{LineString{2,Float64}}}

    if section ∉ VALID_SECTIONS
        throw(ArgumentError("Invalid section: $section. Use one of $VALID_SECTIONS"))
    end

    if section == :from
        target_point = point_to
    elseif section == :to
        target_point = point_from
    end

    remaining_segment::Vector{LineString{2,Float64}} = linestring_segment_to_keep(
        section,
        target_point,
        segment
    )

    leg_to_update = [point_from, point_to]
    dist::Vector{Float64}, leg::Vector{Vector{LineString{2,Float64}}} = get_feasible_vector(
        leg_to_update, exclusions
    )
    if section == :from
        return (sum(dist), vcat(leg..., remaining_segment...))
    elseif section == :to
        return (sum(dist), vcat(remaining_segment..., leg...))
    end
end

"""
    generate_tender_sorties(
        soln::MSTSolution,
        tmp_wpts::Vector{Point{2,Float64}},
        exclusions::POLY_VEC
    )::Vector{TenderSolution}

Generate tender sorties based on the updated waypoints, including all updated attributes.

# Arguments
- `soln`: The existing MSTSolution containing tenders.
- `tmp_wpts`: Vector of temporary waypoints to update the tender sorties.
- `exclusions`: Exclusion zone polygons for tenders.

# Returns
Updated TenderSolution objects with new waypoints.
"""
function generate_tender_sorties(
    soln::MSTSolution,
    tmp_wpts::Vector{Point{2,Float64}},
    exclusions::POLY_VEC
)::Vector{TenderSolution}
    # Update the tender solutions with the new waypoints
    tenders_old = soln.tenders[end]
    n = length(tenders_old)
    length(tmp_wpts) != (2n + 2) && throw(DimensionMismatch(
        "Expected 2 points per tender plus depot start/end"
    ))

    js = 1:n
    starts_new = tmp_wpts[2 .* js]
    finishes_new = tmp_wpts[2 .* js.+1]

    tenders_new = Vector{TenderSolution}(undef, n)

    for j in js
        tender_old = tenders_old[j]
        start_new, finish_new = starts_new[j], finishes_new[j]

        sortie_start_has_moved = tender_old.start != start_new
        sortie_end_has_moved = tender_old.finish != finish_new

        # If neither start nor finish has moved, keep the old tender solution and continue
        if !sortie_start_has_moved && !sortie_end_has_moved
            tenders_new[j] = tender_old
            continue
        end

        # Rebuild the sorties with new start/finish points
        sorties_new = rebuild_sortie.(
            tender_old.sorties,
            Ref(start_new),
            Ref(finish_new),
            Ref(exclusions),
            Ref(sortie_start_has_moved),
            Ref(sortie_end_has_moved)
        )

        tenders_new[j] = TenderSolution(
            tender_old.id,
            start_new,
            finish_new,
            sorties_new,
        )
    end

    return tenders_new
end

"""
    nearest_neighbour(
        cluster_centroids::DataFrame,
        exclusions_mothership::POLY_VEC,
        exclusions_tender::POLY_VEC,
    )::MothershipSolution
    nearest_neighbour(
        cluster_centroids::DataFrame,
        exclusions_mothership::POLY_VEC,
        exclusions_tender::POLY_VEC,
        current_point::Point{2,Float64},
        ex_ms_route::MothershipSolution,
        cluster_seq_idx::Int64
    )::MothershipSolution

Apply the nearest neighbor algorithm:
- starting from the depot (1st row/col) and returning to the depot, or
- starting from the current point and returning to the depot.

# Arguments
- `cluster_centroids`: DataFrame containing id, lat, lon. Depot has `id=0` in row 1.
- `exclusions_mothership`: Exclusion zone polygons for mothership.
- `exclusions_tender`: Exclusion zone polygons for tenders.
- `current_point`: Point{2, Float64} representing the current location of the mothership.
- `ex_ms_route`: MothershipSolution object containing the existing route.
- `cluster_seq_idx`: Index denoting the mothership position by cluster sequence index.

# Returns
MothershipSolution object containing:
- `cluster_sequence`: Centroid sequence DataFrame (containing id, lat, lon).
- `route`: Route object containing waypoints, distance matrix, and line strings.
    - `nodes`: Vector of Point{2, Float64} waypoints.
    - `dist_matrix`: Distance matrix between centroids.
    - `line_strings`: Vector of LineString objects for each path.
"""
function nearest_neighbour(
    cluster_centroids::DataFrame,
    exclusions_mothership::POLY_VEC,
    exclusions_tender::POLY_VEC,
)::MothershipSolution
    # TODO: Use vector rather than DataFrame for cluster_centroids
    # adjust_waypoints to ensure not within exclusion zones - allows for feasible path calc
    feasible_centroids::Vector{Point{2,Float64}} = adjust_waypoint.(
        Point{2,Float64}.(cluster_centroids.lon, cluster_centroids.lat),
        Ref(exclusions_mothership)
    )

    cluster_centroids[!, :lon] = [pt[1] for pt in feasible_centroids]
    cluster_centroids[!, :lat] = [pt[2] for pt in feasible_centroids]

    # Create distance matrix between feasible nodes - cluster centroids
    dist_matrix = get_feasible_matrix(
        feasible_centroids,
        exclusions_mothership
    )[1]

    tour_length = size(dist_matrix, 1)
    visited = falses(tour_length)
    tour = Vector{Int64}(undef, tour_length + 1)
    total_distance = 0.0

    idx = 1
    current_location = 1
    tour[1] = current_location
    visited[current_location] = true

    while idx < tour_length
        idx += 1
        distances = dist_matrix[current_location, :]
        distances[visited] .= Inf
        nearest_idx = argmin(distances)

        tour[idx] = nearest_idx
        total_distance += dist_matrix[current_location, nearest_idx]
        visited[nearest_idx] = true
        current_location = nearest_idx
    end

    # Return to the depot and adjust cluster_sequence to zero-based indexing
    tour[end] = 1
    total_distance += dist_matrix[current_location, 1]
    tour .-= 1

    ordered_centroids = cluster_centroids[[findfirst(==(id), cluster_centroids.id) for id in tour], :]

    # combine exclusions for mothership and tenders
    exclusions_all = vcat(exclusions_mothership, exclusions_tender)
    waypoints = get_waypoints(ordered_centroids, exclusions_all)

    # Calc feasible path between waypoints.
    waypoint_dist_vector, waypoint_path_vector = get_feasible_vector(
        waypoints.waypoint, exclusions_mothership
    )
    path = vcat(waypoint_path_vector...)

    return MothershipSolution(
        cluster_sequence=ordered_centroids,
        route=Route(waypoints.waypoint, waypoint_dist_vector, path)
    )
end
function nearest_neighbour(
    cluster_centroids::DataFrame,
    exclusions_mothership::POLY_VEC,
    exclusions_tender::POLY_VEC,
    current_point::Point{2,Float64},
    ex_ms_route::MothershipSolution,
    cluster_seq_idx::Int64
)::MothershipSolution
    # TODO: Use vector rather than DataFrame for cluster_centroids
    # adjust_waypoints to ensure not within exclusion zones - allows for feasible path calc
    feasible_centroids::Vector{Point{2,Float64}} = adjust_waypoint.(
        Point{2,Float64}.(cluster_centroids.lon, cluster_centroids.lat),
        Ref(exclusions_mothership)
    )
    cluster_centroids[!, :lon] = [pt[1] for pt in feasible_centroids]
    cluster_centroids[!, :lat] = [pt[2] for pt in feasible_centroids]

    # Create distance matrix between start, end, and feasible cluster centroids
    dist_matrix = get_feasible_matrix(
        vcat([current_point], feasible_centroids),
        exclusions_mothership
    )[1]
    centroid_matrix = dist_matrix[3:end, 3:end]
    current_dist_vector = dist_matrix[1, 3:end]
    end_dist_vector = dist_matrix[2, 3:end]

    tour_length = size(centroid_matrix, 1)
    visited = falses(tour_length)
    tour = Vector{Int64}(undef, tour_length)

    current_location = argmin(current_dist_vector)
    total_distance = current_dist_vector[current_location]
    tour[1] = current_location
    visited[current_location] = true
    idx = 1

    while !all(visited)
        idx += 1
        distances = centroid_matrix[current_location, :]
        distances[visited] .= Inf
        nearest_idx = argmin(distances)

        tour[idx] = nearest_idx
        total_distance += distances[nearest_idx]
        visited[nearest_idx] = true
        current_location = nearest_idx
    end

    # Return to the depot and adjust cluster_sequence to zero-based indexing
    push!(tour, 0)
    total_distance += end_dist_vector[current_location]

    sequenced_remaining_centroids = cluster_centroids[tour.+1, :]

    # combine exclusions for mothership and tenders
    exclusions_all = vcat(exclusions_mothership, exclusions_tender)
    waypoints = get_waypoints(
        current_point,
        vcat(DataFrame(ex_ms_route.cluster_sequence[cluster_seq_idx, :]), sequenced_remaining_centroids),
        exclusions_all
    )

    cluster_sequence = vcat(
        ex_ms_route.cluster_sequence[1:cluster_seq_idx, :],
        sequenced_remaining_centroids
    )
    ex_path = ex_ms_route.route.line_strings[
        1:findfirst(
        x -> x == ex_ms_route.route.nodes[2*cluster_seq_idx-1],
        getindex.(getproperty.(ex_ms_route.route.line_strings, :points), 2)
    )]

    # Calc feasible path between waypoints.
    waypoint_dist_vector, waypoint_path_vector = get_feasible_vector(
        waypoints.waypoint, exclusions_mothership
    )
    new_path = vcat(waypoint_path_vector...)
    full_path = vcat(ex_path, new_path...)

    return MothershipSolution(
        cluster_sequence=cluster_sequence,
        route=Route(
            vcat(
                ex_ms_route.route.nodes[1:2*cluster_seq_idx-2],
                waypoints.waypoint
            ),
            waypoint_dist_vector,
            full_path
        )
    )
end

"""
Core two-opt optimization function that can be used by all solution types.
"""
function optimize_route_two_opt(
    route::Vector{Int},
    dist_matrix::AbstractArray{Float64},
    distance_fn::Function;
    start_idx::Int=1,
    end_idx::Int=length(route)
)::Vector{Int}
    best_route = copy(route)
    best_distance = distance_fn(best_route, dist_matrix)
    improved = true

    while improved
        improved = false
        for i in start_idx+1:(end_idx-1)
            for j in (i+1):end_idx
                new_route = two_opt_swap(best_route, i, j)
                new_distance = distance_fn(new_route, dist_matrix)

                if new_distance < best_distance
                    best_route = new_route
                    best_distance = new_distance
                    improved = true
                    @debug "Improved by swapping $(i) and $(j)"
                end
            end
        end
    end

    return best_route
end

"""
Process optimized route for mothership solutions.
"""
function process_mothership_route(
    optimized_route::Vector{Int},
    cluster_centroids::DataFrame,
    exclusions_mothership::POLY_VEC,
    exclusions_tender::POLY_VEC,
    existing_waypoints::Union{Vector,Nothing}=nothing,
    cluster_seq_idx::Union{Int,Nothing}=nothing
)::MothershipSolution
    # Re-orient route to start from and end at the depot, and adjust to zero-based indexing
    processed_route = orient_route(optimized_route)
    push!(processed_route, processed_route[1])
    processed_route .-= 1

    # Handle partial route optimization
    if cluster_seq_idx !== nothing
        processed_route = orient_route(
            processed_route,
            cluster_centroids.id[1:cluster_seq_idx]
        )
    end

    ordered_nodes = cluster_centroids[
        [findfirst(==(id), cluster_centroids.id) for id in processed_route], :
    ]

    exclusions_all = vcat(exclusions_mothership, exclusions_tender)
    waypoints = get_waypoints(ordered_nodes, exclusions_all)

    # Handle waypoint merging for partial optimization
    final_waypoints = if existing_waypoints !== nothing && cluster_seq_idx !== nothing
        vcat(
            existing_waypoints[1:2*cluster_seq_idx-1],
            waypoints.waypoint[2*cluster_seq_idx:end]
        )
    else
        waypoints.waypoint
    end

    waypoint_dist_vec, waypoint_path_vector = get_feasible_vector(
        final_waypoints,
        exclusions_mothership
    )

    path = vcat(waypoint_path_vector...)

    return MothershipSolution(
        cluster_sequence=ordered_nodes,
        route=Route(final_waypoints, waypoint_dist_vec, path)
    )
end

"""
    two_opt(
        ms_soln_current::MothershipSolution,
        exclusions_mothership::Vector{IGeometry{wkbPolygon}},
        exclusions_tender::Vector{IGeometry{wkbPolygon}},
    )::MothershipSolution

Apply the 2-opt heuristic to improve current routes by uncrossing crossed links between
waypoints for the whole route.
"""
function two_opt(
    ms_soln_current::MothershipSolution,
    exclusions_mothership::POLY_VEC,
    exclusions_tender::POLY_VEC,
)::MothershipSolution
    cluster_centroids = ms_soln_current.cluster_sequence

    ordered_centroids = sort(cluster_centroids[1:end-1, :], :id)  # drop last row
    centroid_nodes = Point{2,Float64}.(ordered_centroids.lon, ordered_centroids.lat)
    push!(centroid_nodes, first(centroid_nodes))  # return to depot

    dist_matrix::Matrix{Float64} = get_feasible_matrix(centroid_nodes, exclusions_mothership)[1]

    # If depot is last row, remove
    if cluster_centroids.id[1] == cluster_centroids.id[end]
        cluster_centroids = cluster_centroids[1:end-1, :]
    end

    # Initialize route as ordered waypoints
    initial_route = cluster_centroids.id .+ 1

    # Optimize the entire route
    optimized_route = optimize_route_two_opt(
        initial_route,
        dist_matrix,
        route_distance
    )

    return process_mothership_route(
        optimized_route,
        cluster_centroids,
        exclusions_mothership,
        exclusions_tender
    )
end

"""
    two_opt(
        ms_soln_current::MothershipSolution,
        exclusions_mothership::Vector{IGeometry{wkbPolygon}},
        exclusions_tender::Vector{IGeometry{wkbPolygon}},
        cluster_seq_idx::Int64
    )::MothershipSolution

Apply the 2-opt heuristic to improve current routes by uncrossing crossed links between
waypoints from the current cluster position to the depot.
"""
function two_opt(
    ms_soln_current::MothershipSolution,
    exclusions_mothership::POLY_VEC,
    exclusions_tender::POLY_VEC,
    cluster_seq_idx::Int64
)::MothershipSolution
    cluster_centroids = ms_soln_current.cluster_sequence

    ordered_centroids = sort(cluster_centroids[1:end-1, :], :id)  # drop last row
    centroid_nodes = Point{2,Float64}.(ordered_centroids.lon, ordered_centroids.lat)
    push!(centroid_nodes, first(centroid_nodes))  # return to depot

    dist_matrix::Matrix{Float64} = get_feasible_matrix(centroid_nodes, exclusions_mothership)[1]

    # If depot is last row, remove
    if cluster_centroids.id[1] == cluster_centroids.id[end]
        cluster_centroids = cluster_centroids[1:end-1, :]
    end

    # Initialize route as ordered waypoints
    initial_route = cluster_centroids.id .+ 1

    # Optimize only from cluster_seq_idx to the end
    optimized_route = optimize_route_two_opt(
        initial_route,
        dist_matrix,
        route_distance;
        start_idx=cluster_seq_idx,
        end_idx=length(initial_route)
    )

    return process_mothership_route(
        optimized_route,
        cluster_centroids,
        exclusions_mothership,
        exclusions_tender,
        ms_soln_current.route.nodes,
        cluster_seq_idx
    )
end

"""
    two_opt(
        tender_soln_current::TenderSolution,
        exclusions_tender::DataFrame,
    )::TenderSolution

Apply the 2-opt heuristic to improve current tender sortie routes.
"""
function two_opt(
    tender_soln_current::TenderSolution,
    exclusions_tender::POLY_VEC,
)::TenderSolution
    sorties_current::Vector{Route} = tender_soln_current.sorties
    sorties_new::Vector{Route} = Vector{Route}(undef, length(sorties_current))

    for (idx, sortie) in enumerate(sorties_current)
        if length(sortie.nodes) ≤ 1
            # If there is only 1 node between start and finish, no change is possible
            # update to distance vector
            dist_vector::Vector{Float64} = get_distance_vector(sortie.dist_matrix)
            sorties_new[idx] = Route(
                sortie.nodes,
                dist_vector,
                vcat(sortie.line_strings...)
            )
            continue
        end

        nodes = vcat([tender_soln_current.start], sortie.nodes, [tender_soln_current.finish])
        n = length(nodes)
        initial_sequence = collect(1:n)

        # Optimize sequence, but don't modify start (index 1) and finish (index n)
        optimized_sequence = optimize_route_two_opt(
            initial_sequence,
            sortie.dist_matrix,
            tender_sortie_dist;
            start_idx=2,
            end_idx=n - 1
        )

        # Build new sortie ordering, ignoring start and finish
        nodes_new = nodes[optimized_sequence[2:end-1]]

        # Recompute feasible path through exclusions
        dist_vector, path_vector = get_feasible_vector(
            nodes[optimized_sequence],
            exclusions_tender
        )
        new_path = vcat(path_vector...)

        sorties_new[idx] = Route(nodes_new, dist_vector, new_path)
    end

    return TenderSolution(
        tender_soln_current.id,
        tender_soln_current.start,
        tender_soln_current.finish,
        sorties_new,
    )
end

"""
    two_opt_swap(route::Vector{Int64}, i::Int, j::Int)::Vector{Int64}

Swap two points in a route.

# Arguments
- `route`: Vector of cluster indices.
- `i`: Index of the first point to swap.
- `j`: Index of the second point to swap.

# Returns
Route with swapped points.
"""
function two_opt_swap(route::Vector{Int64}, i::Int, j::Int)::Vector{Int64}
    return vcat(route[1:(i-1)], reverse(route[i:j]), route[(j+1):end])
end

"""
    orient_route(route::Vector{Int64})::Vector{Int64}

Orient the route such that it starts with the depot (1).

# Arguments
- `route`: Vector of cluster indices.

# Returns
Route starting with the depot.
"""
function orient_route(route::Vector{Int64})::Vector{Int64}
    idx = findfirst(==(1), route)
    return vcat(route[idx:end], route[1:idx-1])
end
function orient_route(route::Vector{Int64}, start_segment::Vector{Int64})::Vector{Int64}
    if length(start_segment) < 2
        throw(ArgumentError("start_segment must have at least 2 elements"))
    end
    if start_segment[2] == route[2]
        return route
    elseif start_segment[2] == route[end-1]
        return reverse(route)
    else
        throw(ArgumentError("start_segment does not match route"))
    end
    return
end

"""
    tender_sequential_nearest_neighbour(
        cluster::Cluster,
        waypoints::NTuple{2,Point{2,Float64}},
        n_tenders::Int8,
        t_cap::Int16,
        exclusions::Vector{IGeometry{wkbPolygon}}
    )::TenderSolution

Assign nodes to tenders sequentially (stop-by-stop) based on nearest neighbor.

# Arguments
- `cluster`: Cluster object containing nodes.
- `waypoints`: Tuple of start and end waypoints.
- `n_tenders`: Number of tenders.
- `t_cap`: Tender capacity.
- `exclusions`: Exclusion zone polygons.

# Returns
- `solution`: TenderRoutingSolution object containing:
    - `cluster_id`: Cluster ID.
    - `start`: Start waypoint.
    - `finish`: End waypoint.
    - `sorties`: Vector of Sortie objects containing nodes and sortie distance.
        - `nodes`: Vector of Point{2, Float64} waypoints.
        - `dist_matrix`: Distance matrix between nodes.
        - `line_strings`: Vector of LineString objects for each path.
"""
function tender_sequential_nearest_neighbour(
    cluster::Cluster,
    waypoints::NTuple{2,Point{2,Float64}},
    n_tenders::Int8,
    t_cap::Int16,
    exclusions::POLY_VEC
)::TenderSolution
    nodes = [[waypoints[1]]; cluster.nodes]
    full_nodes = vcat(nodes, [waypoints[2]])

    # Compute the full feasible matrix and associated paths once.
    dist_matrix, path_matrix = get_feasible_matrix(full_nodes, exclusions)

    tender_tours = [Vector{Int}(undef, t_cap) for _ in 1:n_tenders]
    tour_lengths = zeros(Int, n_tenders)
    n_nodes = length(nodes)

    visited = falses(n_nodes)
    visited[1] = true
    visited[2:end] .= isinf.(dist_matrix[1, 2:n_nodes])

    # sequentially assign closest un-visited nodes stop-by-stop, tender-by-tender
    for _ in 1:t_cap
        for t in 1:n_tenders
            if all(visited)
                break
            end

            current_node = tour_lengths[t] == 0 ? 1 : tender_tours[t][tour_lengths[t]]

            distances = dist_matrix[current_node, 1:n_nodes]
            distances[visited] .= Inf
            nearest_idx = argmin(distances)

            tour_lengths[t] += 1
            tender_tours[t][tour_lengths[t]] = nearest_idx

            visited[nearest_idx] = true
        end
    end

    # delete excess elements and remove empty tours by slicing tender_tours to tour_lengths
    tender_tours = filter(
        !isempty,
        [tour[1:tour_lengths[t]] for (t, tour) in enumerate(tender_tours)]
    )

    routes = Vector{Route}(undef, length(tender_tours))
    for (t, tour) in enumerate(tender_tours)
        # Route indices incl waypoints: (1), tour nodes, end (last index of full_nodes)
        route_indices = Vector{Int}(undef, length(tour) + 2)
        route_indices[1] = 1
        route_indices[2:end-1] = tour
        route_indices[end] = length(full_nodes)

        # Submatrix for this route
        route_matrix = dist_matrix[route_indices, route_indices]

        route_paths = [
            (i < j ? path_matrix[i, j] : path_matrix[j, i])
            for (i, j) in zip(route_indices[1:end-1], route_indices[2:end])
        ]

        routes[t] = Route(
            cluster.nodes[tour.-1],
            route_matrix,
            vcat(route_paths...)
        )
    end

    initial_sortie = TenderSolution(
        cluster.id,
        waypoints[1],
        waypoints[2],
        routes,
    )
    improved_sortie = two_opt(
        initial_sortie,
        exclusions
    )

    return improved_sortie
end
