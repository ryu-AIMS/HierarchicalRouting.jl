
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
