
"""
    get_feasible_matrix(
        nodes::Vector{Point{2,Float64}},
        exclusions::Vector{IGeometry{wkbPolygon}}
    )::Tuple{Matrix{Float64} Matrix{Vector{LineString{2,Float64}}}}

Create matrices of distances and paths for feasible routes between waypoints accounting for
    (avoiding) environmental constraints.

# Arguments
- `nodes`: Vector of lat long tuples.
- `exclusions`: Exclusion zones representing vehicle's environmental constraints.

# Returns
- `dist_matrix`: A matrix of distances (in metres) between waypoints.
- `path_matrix`: A vector of paths between
    waypoints.
"""
function get_feasible_matrix(
    nodes::Vector{Point{2,Float64}},
    exclusions::POLY_VEC
)::Tuple{Matrix{Float64},Matrix{Vector{LineString{2,Float64}}}}
    n_points = length(nodes)
    dist_matrix = zeros(Float64, n_points, n_points)
    path_matrix = fill(Vector{LineString{2,Float64}}(), n_points, n_points)

    for point_j_idx in 1:n_points, point_i_idx in 1:point_j_idx-1
        point_nodes = getindex(nodes, [point_i_idx, point_j_idx])
        if point_nodes[1] != point_nodes[2]
            # TODO: Process elsewhere
            # Check if any of the points are within an exclusion zone
            if any(point_in_exclusion.(point_nodes, [exclusions]))
                dist_matrix[point_i_idx, point_j_idx] =
                    dist_matrix[point_j_idx, point_i_idx] =
                        Inf
            else
                dist_matrix[point_i_idx, point_j_idx],
                path_matrix[point_i_idx, point_j_idx] =
                    shortest_feasible_path(
                        point_nodes[1], point_nodes[2], exclusions
                    )
                dist_matrix[point_j_idx, point_i_idx] =
                    dist_matrix[point_i_idx, point_j_idx]
            end
        end
    end

    return dist_matrix, path_matrix
end

"""
    get_feasible_vector(
        nodes::Vector{Point{2,Float64}}, exclusions::Vector{IGeometry{wkbPolygon}}
    )::Tuple{Vector{Float64},Vector{Vector{LineString{2,Float64}}}}

Create vectors of feasible:
- distances, and
- paths
between sequential nodes, avoiding exclusions.

# Arguments
- `nodes`: Vector of lat long tuples.
- `exclusions`: Geometries of exclusion zones representing given vehicle's cumulative
environmental constraints.

# Returns
- `dist_vector`: A vector of distances (in metres) between waypoints.
- `path_vector`: A vector of paths between waypoints, represented as LineStrings.
"""
function get_feasible_vector(
    nodes::Vector{Point{2,Float64}}, exclusions::POLY_VEC
)::Tuple{Vector{Float64},Vector{Vector{LineString{2,Float64}}}}
    n_points = length(nodes) - 1
    dist_vector = zeros(Float64, n_points)
    path_vector = fill(Vector{LineString{2,Float64}}(), n_points)

    for point_i_idx in 1:n_points
        point_nodes = getindex(nodes, [point_i_idx, point_i_idx + 1])
        if nodes[point_i_idx] != nodes[point_i_idx+1]
            # TODO: Process elsewhere
            # Check if any of the points are within an exclusion zone
            if !any(point_in_exclusion.(point_nodes, Ref(exclusions)))
                dist_vector[point_i_idx], path_vector[point_i_idx] =
                    shortest_feasible_path(
                        point_nodes[1], point_nodes[2], exclusions
                    )
            else
                dist_vector[point_i_idx] = Inf
            end
        end
    end
    any(isinf.(dist_vector)) && throw(DomainError(dist_vector[findall(isinf, dist_vector)],
        "Distance vector contains Inf values, indicating a waypoint in an exclusion zone."))
    return dist_vector, path_vector
end

"""
    get_feasible_distances(
        current_location::Point{2,Float64},
        targets::Vector{Point{2,Float64}},
        exclusions::Vector{IGeometry{wkbPolygon}}
    )::Vector{Float64}

Create a vector of feasible distances between the current location and each target point.

# Arguments
- `current_location`: Current location of the vehicle.
- `targets`: Vector of target points.
- `exclusions`: Exclusion zones representing vehicle's environmental constraints.

# Returns
- A vector of feasible distances (m) between the current location and each target point.
"""
function get_feasible_distances(
    current_location::Point{2,Float64},
    targets::Vector{Point{2,Float64}},
    exclusions::POLY_VEC
)::Vector{Float64}
    n_points = length(targets)
    dist_vector = zeros(Float64, n_points)

    if point_in_exclusion(current_location, exclusions)
        #! Should never be in exclusion zone
        dist_vector .= Inf
        return dist_vector
    end

    for point_i_idx in 1:n_points
        point_nodes = [current_location, targets[point_i_idx]]

        # Check if any of the points are within an exclusion zone
        if point_in_exclusion(targets[point_i_idx], exclusions)
            #! Should never be in exclusion zone
            dist_vector[point_i_idx] = Inf
        else
            dist_vector[point_i_idx], _ = shortest_feasible_path(
                point_nodes[1], point_nodes[2], exclusions
            )
        end
    end

    return dist_vector
end

"""
    point_in_exclusion(point::Point{2,Float64}, exclusion::IGeometry{wkbPolygon})::Bool
    point_in_exclusion(point::Point{2,Float64}, exclusions::Vector{IGeometry{wkbPolygon}})::Bool

Check if a point is within an exclusion zone.

# Arguments
- `point`: Point to check.
- `exclusion::IGeometry`: One exclusion zone polygon/geometry from a DataFrame.
- `exclusions::Vector{IGeometry{wkbPolygon}}`: A vector of exclusion zone polygons/geometry.

# Returns
- `true` if point is within an exclusion zone, `false` otherwise.
"""
function point_in_exclusion(point::Point{2,Float64}, exclusion::IGeometry{wkbPolygon})::Bool
    return any(AG.contains(exclusion, AG.createpoint(point.data)))
end
function point_in_exclusion(point::Point{2,Float64}, exclusions::POLY_VEC)::Bool
    return any(AG.contains.(exclusions, Ref(AG.createpoint(point.data))))
end

"""
    containing_exclusion(point::Point{2,Float64}, exclusions::Vector{IGeometry{wkbPolygon}})::Int

Return the index of the exclusion zone that contains the point.

# Arguments
- `point`: Point to check.
- `exclusions`: A DataFrame containing exclusion zone polygons.

# Returns
- Index of the first exclusion zone that contains the point, or 0 if not found.
"""
function containing_exclusion(point::Point{2,Float64}, exclusions::POLY_VEC)::Int
    point_ag = AG.createpoint(point.data)
    exclusion_idx = findfirst(AG.contains.(exclusions, Ref(point_ag)))
    return isnothing(exclusion_idx) ? 0 : exclusion_idx
end

"""
    shortest_feasible_path(
        initial_point::Point{2,Float64},
        final_point::Point{2,Float64},
        exclusions::Vector{IGeometry{wkbPolygon}}
    )::Tuple{Float64,Vector{LineString{2,Float64}}}

Find the shortest feasible path between two points.
- Build network of all polygon vertices that recursively intersect with line to finish, from
start pt and any other intersecting polygons.
- Build graph from network of points.
- Use A* algorithm to find shortest path.

# Arguments
- `initial_point`: Starting point of path.
- `final_point`: Ending point of path.
- `exclusions`: A vector of exclusion zones.

# Returns
- The distance of the shortest feasible path (in metres).
- The shortest feasible path as a vector of LineStrings.
"""
function shortest_feasible_path(
    initial_point::Point{2,Float64},
    final_point::Point{2,Float64},
    exclusions::POLY_VEC;
    _tried_reverse::Bool=false
)::Tuple{Float64,Vector{LineString{2,Float64}}}
    final_exclusion_idx = point_in_convexhull(final_point, exclusions)
    initial_exclusion_idx = point_in_convexhull(initial_point, exclusions)

    points_from = Point{2,Float64}[]
    points_to = Point{2,Float64}[]
    exclusion_idxs = Int64[]

    if iszero(initial_exclusion_idx)
        # If initial point is not within an exclusion zone
        build_network!(
            points_from,
            points_to,
            exclusion_idxs,
            initial_point,
            final_point,
            exclusions,
            final_exclusion_idx
        )
    elseif is_visible(initial_point, final_point, exclusions)
        # If initial point is within an exclusion zone and visible to final point
        push!(points_from, initial_point)
        push!(points_to, final_point)
        push!(exclusion_idxs, initial_exclusion_idx)
    else
        # Collect all visible polygon vertices
        initial_polygon::IGeometry{wkbPolygon} = exclusions[initial_exclusion_idx]
        poly_vertices::Vector{Point{2,Float64}} = collect_polygon_vertices(
            initial_polygon
        )
        visible_vertices = unique(
            filter(
                v ->
                    v != initial_point &&
                        is_visible(initial_point, v, exclusions),
                poly_vertices
            )
        )
        if !isempty(visible_vertices)
            # Connect initial point to every visible polygon vertex
            n_vis_verts = length(visible_vertices)

            points_from = fill(initial_point, n_vis_verts)
            points_to = visible_vertices
            exclusion_idxs = fill(initial_exclusion_idx, n_vis_verts)

            #? Connect initial point to all visible vertices on final polygon?
        end

        # Reuseable cache for mask to avoid repeated allocations
        visibility_mask = falses(length(poly_vertices))

        # For each polygon vertex, add edges to all visible vertices
        for i in poly_vertices
            # Check visibility of all polygon vertices from the current vertex `i`
            visibility_mask .= is_visible.(Ref(i), poly_vertices, Ref(exclusions))
            vis_pts_from_i = poly_vertices[visibility_mask]
            n_vis_pts_from_i = length(vis_pts_from_i)

            if !isempty(vis_pts_from_i)
                append!(points_from, fill(i, n_vis_pts_from_i))
                append!(points_to, vis_pts_from_i)
                append!(
                    exclusion_idxs,
                    fill(initial_exclusion_idx, n_vis_pts_from_i)
                )
            end
            if is_visible(i, final_point, exclusions)
                # Connect initial point to final point if visible
                push!(points_from, i)
                push!(points_to, final_point)
                push!(exclusion_idxs, initial_exclusion_idx)
            end
        end

        # Connect widest_verts (from final point) to existing network
        #! Note that final/initial points are inverted
        widest_verts = find_widest_points(
            final_point,
            initial_point,
            [exclusions[initial_exclusion_idx]]
        )[1]

        # visited = Set(vcat(points_from, points_to))
        build_network!.(
            Ref(points_from),
            Ref(points_to),
            Ref(exclusion_idxs),
            widest_verts,
            Ref(final_point),
            Ref(exclusions),
            Ref(final_exclusion_idx);
            # visited=visited
        )
    end

    # Build graph from network of points
    graph, idx_to_point, initial_point_idx, final_point_idx = build_graph(
        points_from,
        points_to,
        exclusions,
        iszero(final_exclusion_idx) ? Point{2,Float64}[] : exclusions[final_exclusion_idx],
        initial_point,
        final_point
    )

    # Use A* algorithm to find shortest path
    path = a_star(graph, initial_point_idx, final_point_idx, graph.weights)
    if iszero(length(path)) && !_tried_reverse
        dist_rev, segs_rev = shortest_feasible_path(
            final_point,
            initial_point,
            exclusions;
            _tried_reverse=true
        )
        # Flip segments original forward direction (and order)
        segs_fwd = reverse([LineString([last(ls), first(ls)]) for ls in getfield.(segs_rev, :points)])
        return dist_rev, segs_fwd
    elseif iszero(length(path))
        throw(ErrorException(
            "No path ($initial_point -> $final_point) because network/graph incomplete"
        ))
    end
    dist = sum(graph.weights[p.src, p.dst] for p in path)

    linestring_path::Vector{LineString} = [
        LineString([idx_to_point[s.src], idx_to_point[s.dst]]) for s in path
    ]

    return dist, linestring_path
end

"""
    point_in_convexhull(
        point::Point{2,Float64},
        exclusions::Vector{IGeometry{wkbPolygon}}
    )::Int

Check if a point is within a convex hull of exclusion zones.

# Arguments
- `point`: Point to check.
- `exclusions`: Exclusion zones.

# Returns
- Index of exclusion zone if point is within a convex hull, 0 otherwise.
"""
function point_in_convexhull(point::Point{2,Float64}, exclusions::POLY_VEC)
    point_ag::AG.IGeometry{wkbPoint} = AG.createpoint(point.data)
    convex_exclusions_ag::POLY_VEC = AG.convexhull.(exclusions)

    point_in_exclusion_zone = AG.contains.(convex_exclusions_ag, [point_ag])

    # If point is within an exclusion zone, add all polygon vertices to graph
    # ? Consider cases where point is within the convex hull of multiple polygons
    # findall(final_point_in_exclusion_zone)
    exclusion_index = findfirst(point_in_exclusion_zone)

    return isnothing(exclusion_index) ? 0 : exclusion_index
end

"""
    build_network!(
        points_from::Vector{Point{2,Float64}},
        points_to::Vector{Point{2,Float64}},
        exclusion_idxs::Vector{Int},
        current_point::Point{2,Float64},
        final_point::Point{2,Float64},
        exclusions::Vector{IGeometry{wkbPolygon}},
        final_exclusion_idx::Int
    )::Nothing

Build a network of points to connect to each other.
Vectors `points_from`, `points_to`, and `exclusion_idxs` are modified in place.

# Arguments
- `points_from`: Vector of points to connect from.
- `points_to`: Vector of points to connect to.
- `exclusion_idxs`: Vector of exclusion zone indices.
- `current_point`: Current point to connect from.
- `final_point`: Final point to connect to.
- `exclusions`: DataFrame containing exclusion zones.
- `final_exclusion_idx`: Index of the final exclusion zone.
"""
function build_network!(
    points_from::Vector{Point{2,Float64}},
    points_to::Vector{Point{2,Float64}},
    exclusion_idxs::Vector{Int64},
    current_point::Point{2,Float64},
    final_point::Point{2,Float64},
    exclusions::POLY_VEC,
    final_exclusion_idx::Int64;
    visited::Set{Point{2,Float64}}=Set{Point{2,Float64}}(),
    existing_edges::Set{Tuple{Point{2,Float64},Point{2,Float64}}}=Set{Tuple{Point{2,Float64},Point{2,Float64}}}()
)::Nothing
    if current_point == final_point
        return nothing
    end

    candidates, next_exclusion_idxs = find_widest_points(
        current_point,
        final_point,
        exclusions
    )

    # Preallocate enough extra capacity assuming candidate point is added
    length_vectors_new = length(points_from) + length(candidates)
    sizehint!(points_from, length_vectors_new)
    sizehint!(points_to, length_vectors_new)
    sizehint!(exclusion_idxs, length_vectors_new)

    for (candidate, next_exclusion_idx) in zip(candidates, next_exclusion_idxs)
        # Skip invalid or previously visited candidates
        should_skip = isnothing(candidate) ||
                      candidate == current_point ||
                      candidate ∈ visited ||
                      (current_point, candidate) ∈ existing_edges
        if should_skip
            continue
        end

        # Record new point/edge
        push!(existing_edges, (current_point, candidate))
        push!(points_from, current_point)
        push!(points_to, candidate)
        push!(exclusion_idxs, next_exclusion_idx)

        # Continue building network if not at final exclusion zone
        keep_building = candidate ∉ visited &&
                        !(next_exclusion_idx == final_exclusion_idx &&
                          !isnothing(next_exclusion_idx))
        if keep_building
            build_network!(
                points_from,
                points_to,
                exclusion_idxs,
                candidate,
                final_point,
                exclusions,
                final_exclusion_idx;
                # visited=visited,
                existing_edges=existing_edges
            )
        end
    end
    # TODO: check if path from current_point to next/intermediate point is feasible (no intersecting polygons in between)
end

"""
    build_graph(
        points_from::Vector{Point{2,Float64}},
        points_to::Vector{Point{2,Float64}},
        exclusions::Vector{IGeometry{wkbPolygon}},
        final_polygon::IGeometry{wkbPolygon},
        initial_point::Point{2,Float64},
        final_point::Point{2,Float64}
    )::Tuple{SimpleWeightedGraph{Int64,Float64}, Vector{Point{2,Float64}},Int64,Int64}

Build a graph network from a set of points.
If polygon is provided, connect it to visible points in graph and visible vertices.

# Arguments
- `points_from`: Vector of points to connect from.
- `points_to`: Vector of points to connect to.
- `exclusions`: DataFrame containing exclusion zones.
- `final_polygon`: Polygon to connect to.
- `initial_point`: Initial point to connect from.
- `final_point`: Final point to connect to.

# Returns
- Graph object with edges connecting points.
- Vector of unique points.
- Index of initial point. Used in a_star to specify starting point.
- Index of final point. Used in a_star to specify target point.
"""
function build_graph(
    points_from::Vector{Point{2,Float64}},
    points_to::Vector{Point{2,Float64}},
    exclusions::POLY_VEC,
    final_polygon::IGeometry{wkbPolygon},
    initial_point::Point{2,Float64},
    final_point::Point{2,Float64}
)::Tuple{SimpleWeightedGraph{Int64,Float64},Vector{Point{2,Float64}},Int64,Int64}
    final_poly_verts::Vector{Point{2,Float64}} = collect_polygon_vertices(final_polygon)

    # Collect connected points from network
    #? Use `Set`?
    chain_points::Vector{Point{2,Float64}} = unique(
        vcat(points_from, points_to, [final_point])
    )

    unique_points::Vector{Point{2,Float64}} = unique(vcat(chain_points, final_poly_verts))
    n_points::Int = length(unique_points)

    # Create the graph with one vertex per unique point.
    graph = SimpleWeightedGraph(n_points)

    # Map points to indices
    pt_to_idx = Dict{Point{2,Float64},Int}(unique_points .=> 1:n_points)

    # Add candidate (chain) edges.
    add_edge!.(Ref(graph),
        map(p -> pt_to_idx[p], points_from),
        map(q -> pt_to_idx[q], points_to),
        haversine.(points_from, points_to)
    )

    # Only add polygon edges and extra visible connections if final_point is not in network
    if final_point ∉ points_to
        # Cache to reuse for masking visible vertices
        visibility_mask = falses(length(final_poly_verts))

        # Add edges connecting any visible polygon vertices to other polyogn vertices
        for i in final_poly_verts
            visibility_mask .= (final_poly_verts .!= i) .&
                               is_visible.(Ref(i), final_poly_verts, Ref(exclusions))

            visible_points = final_poly_verts[visibility_mask]

            if !isempty(visible_points)
                # Compute indices for all visible points and distances in one go.
                add_edge!.(
                    Ref(graph),
                    Ref(pt_to_idx[i]),
                    map(pt -> pt_to_idx[pt], visible_points),
                    haversine.(Ref(i), visible_points)
                )
            end
        end

        # Connect all chain points to visible polygon vertices
        for i in final_poly_verts
            visible_points = filter(pt -> is_visible(i, pt, exclusions), chain_points)
            if !isempty(visible_points)
                # Compute indices for all visible points and distances in one go.
                add_edge!.(
                    Ref(graph),
                    Ref(pt_to_idx[i]),
                    map(pt -> pt_to_idx[pt], visible_points),
                    haversine.(Ref(i), visible_points)
                )
            end
        end

        # Connect polygon vertices if not complete loop
        if final_poly_verts[1] != final_poly_verts[end]
            add_edge!(
                graph,
                pt_to_idx[final_poly_verts[end]],
                pt_to_idx[final_poly_verts[1]],
                haversine(final_poly_verts[end], final_poly_verts[1])
            )
        end
    end

    return graph, unique_points, pt_to_idx[initial_point], pt_to_idx[final_point]
end
function build_graph(
    points_from::Vector{Point{2,Float64}},
    points_to::Vector{Point{2,Float64}},
    exclusions::POLY_VEC,
    _::Vector{Point{2,Float64}},
    initial_point::Point{2,Float64},
    final_point::Point{2,Float64}
)::Tuple{SimpleWeightedGraph{Int64,Float64},Vector{Point{2,Float64}},Int64,Int64}
    # Collect connected points from network
    chain_points::Vector{Point{2,Float64}} = unique(
        vcat(points_from, points_to, [final_point])
    )
    n_points::Int = length(chain_points)

    # Create the graph with one vertex per unique point.
    graph::SimpleWeightedGraph{Int64,Float64} = SimpleWeightedGraph(n_points)
    visible_mask::Vector{Bool} = is_visible.(chain_points, Ref(final_point), Ref(exclusions))
    visible_idx::Vector{Int64} = findall(visible_mask)
    initial_pnt_idx::Vector{Int64} = findall(chain_points .== initial_point)
    final_pnt_idx::Vector{Int64} = findall(chain_points .== final_point)

    # Add candidate (chain) edges.
    add_edge!.(
        Ref(graph),
        vcat(map(p -> findall(chain_points .== p), points_from)...)::Vector{Int64},
        vcat(map(p -> findall(chain_points .== p), points_to)...)::Vector{Int64},
        haversine.(points_from, points_to)::Vector{Float64}
    )

    # Final_point is not in a convex hull, try to connect it to any visible point
    visible_pts_to_final = chain_points[visible_idx]
    if !isempty(visible_pts_to_final)
        add_edge!.(
            Ref(graph),
            visible_idx,
            final_pnt_idx,
            haversine.(visible_pts_to_final, Ref(final_point))
        )
    end

    return graph, chain_points, first(initial_pnt_idx), first(final_pnt_idx)
end

"""
    collect_polygon_vertices(polygon::IGeometry{wkbPolygon})::Vector{Point{2,Float64}}

Collect all vertices of a polygon.

# Arguments
- `polygon`: Polygon to collect vertices from.

# Returns
Vector of polygon vertices.
"""
function collect_polygon_vertices(polygon::IGeometry{wkbPolygon})::Vector{Point{2,Float64}}
    exterior_ring::IGeometry{wkbLineString} = GI.getgeom(polygon, 1)
    pts::Vector{Point{2,Float64}} = Point{2,Float64}.(GI.coordinates(exterior_ring))

    return unique(pts)
end
