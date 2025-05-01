
"""
    get_feasible_matrix(
        nodes::Vector{Point{2, Float64}},
        exclusions::DataFrame
    )::Tuple{Matrix{Float64}, Matrix{Vector{LineString{2, Float64}}}}

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
    nodes::Vector{Point{2, Float64}},
    exclusions::DataFrame
)::Tuple{Matrix{Float64}, Matrix{Vector{LineString{2, Float64}}}}
    n_points = length(nodes)
    dist_matrix = zeros(Float64, n_points, n_points)
    path_matrix = fill(Vector{LineString{2, Float64}}(), n_points, n_points)

    for point_j_idx in 1:n_points
        for point_i_idx in 1:point_j_idx-1
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
    end

    return dist_matrix, path_matrix
end

"""
    get_feasible_vector(
        nodes::Vector{Point{2, Float64}}, exclusions::DataFrame
    )::Tuple{Vector{Float64}, Vector{Vector{LineString{2, Float64}}}}

Create vectors of feasible:
- distances, and
- paths
between sequential nodes, avoiding exclusions.

# Arguments
- `nodes`: Vector of lat long tuples.
- `exclusions`: DataFrame containing exclusion zones representing given vehicle's cumulative
environmental constraints.

# Returns
- `dist_vector`: A vector of distances (in metres) between waypoints.
- `path_vector`: A vector of paths between waypoints, represented as LineStrings.
"""
function get_feasible_vector(
    nodes::Vector{Point{2, Float64}}, exclusions::DataFrame
)::Tuple{Vector{Float64}, Vector{Vector{LineString{2, Float64}}}}
    n_points = length(nodes)-1
    dist_vector = zeros(Float64, n_points)
    path_vector = fill(Vector{LineString{2, Float64}}(), n_points)

    for point_i_idx in 1:n_points
        point_nodes = getindex(nodes, [point_i_idx, point_i_idx+1])
        if nodes[point_i_idx] != nodes[point_i_idx+1]
            # TODO: Process elsewhere
            # Check if any of the points are within an exclusion zone
            if any(point_in_exclusion.(point_nodes, [exclusions]))
                dist_vector[point_i_idx] = Inf
            else
                dist_vector[point_i_idx], path_vector[point_i_idx] =
                shortest_feasible_path(
                    point_nodes[1], point_nodes[2], exclusions
                )
            end
        end
    end

    return dist_vector, path_vector
end

"""
    get_feasible_distances(
        current_location::Point{2, Float64},
        targets::Vector{Point{2, Float64}},
        exclusions::DataFrame
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
    current_location::Point{2, Float64},
    targets::Vector{Point{2, Float64}},
    exclusions::DataFrame
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
    point_in_exclusion(point::Point{2, Float64}, exclusions::DataFrame)::Bool
    point_in_exclusion(point::Point{2, Float64}, exclusion::AG.IGeometry)::Bool

Check if a point is within an exclusion zone.

# Arguments
- `point`: Point to check.
- `exclusions::DataFrame`: A DataFrame containing exclusion zone polygons.
- `exclusion::AG.IGeometry`: One exclusion zone polygon/geometry from a DataFrame.

# Returns
- `true` if point is within an exclusion zone, `false` otherwise.
"""
function point_in_exclusion(point::Point{2, Float64}, exclusions::DataFrame)::Bool
    point_ag = AG.createpoint(point[1], point[2])
    return any(AG.contains.(exclusions.geometry, [point_ag]))
end
function point_in_exclusion(point::Point{2, Float64}, exclusion::AG.IGeometry)::Bool
return any(AG.contains(exclusion, AG.createpoint(point[1], point[2])))
end

"""
    shortest_feasible_path(
        initial_point::Point{2, Float64},
        final_point::Point{2, Float64},
        exclusions::DataFrame
    )::Tuple{Float64, Vector{LineString{2, Float64}}}

Find the shortest feasible path between two points.
- Build network of all polygon vertices that recursively intersect with line to finish, from
start pt and any other intersecting polygons.
- Build graph from network of points.
- Use A* algorithm to find shortest path.

# Arguments
- `initial_point`: Starting point of path.
- `final_point`: Ending point of path.
- `exclusions`: A DataFrame containing exclusion zone polygons.

# Returns
- The distance of the shortest feasible path (in metres).
- The shortest feasible path as a vector of LineStrings.
"""
function shortest_feasible_path(
    initial_point::Point{2, Float64},
    final_point::Point{2, Float64},
    exclusions::DataFrame,
)::Tuple{Float64, Vector{LineString{2, Float64}}}
    final_exclusion_idx = point_in_convexhull(final_point, exclusions)
    initial_exclusion_idx = point_in_convexhull(initial_point, exclusions)

    points_from = Point{2,Float64}[]
    points_to = Point{2,Float64}[]
    exclusion_idxs = Int[]
    # If initial point is not within an exclusion zone
    if iszero(initial_exclusion_idx)
        build_network!(
            points_from,
            points_to,
            exclusion_idxs,
            initial_point,
            final_point,
            exclusions,
            final_exclusion_idx
        )
    else
        # Collect all visible polygon vertices
        initial_polygon::AG.IGeometry{AG.wkbPolygon} = exclusions[
            initial_exclusion_idx, :geometry
        ]
        poly_vertices::Vector{Point{2, Float64}} = collect_polygon_vertices(initial_polygon)
        visible_vertices = poly_vertices[
            is_visible.(
                Ref(initial_point),
                poly_vertices,
                Ref(exclusions)
            )
        ]

        if !isempty(visible_vertices)
            # Connect initial point to every visible polygon vertex
            n_vis_verts = length(visible_vertices)

            append!(points_from, fill(initial_point, n_vis_verts))
            append!(points_to, visible_vertices)
            append!(exclusion_idxs, fill(initial_exclusion_idx, n_vis_verts))

            # For each polygon vertex, add edges to all visible vertices
            for i in poly_vertices
                # Check visibility of all polygon vertices from the current vertex `i`
                visibility_mask = is_visible.(Ref(i), poly_vertices, Ref(exclusions))

                # Get all visible vertices from `i`
                vis_pts_from_i = poly_vertices[visibility_mask]
                n_vis_pts_from_i = length(vis_pts_from_i)

                if !isempty(vis_pts_from_i)
                    append!(points_from, fill(i, n_vis_pts_from_i))
                    append!(points_to, vis_pts_from_i)
                    append!(exclusion_idxs, fill(initial_exclusion_idx, n_vis_pts_from_i))
                end
            end

        end

        # Get widest points to final point on polygon contianing initial point
        widest_verts = find_widest_points(
            final_point,
            initial_point,
            DataFrame(exclusions[initial_exclusion_idx,:])
        )[1]

        # Continue building network from each of widest vertices to final point
        build_network!.(
            [points_from],
            [points_to],
            [exclusion_idxs],
            widest_verts,
            [final_point],
            [exclusions],
            [final_exclusion_idx]
        )
    end

    # Build graph from network of points
    graph, idx_to_point, initial_point_idx, final_point_idx = build_graph(
        points_from,
        points_to,
        exclusions,
        iszero(final_exclusion_idx) ?
            AG.creategeom(AG.wkbPolygon) : exclusions[final_exclusion_idx, :geometry],
        initial_point,
        final_point
    )

    # Use A* algorithm to find shortest path
    path = a_star(graph, initial_point_idx, final_point_idx, graph.weights)
    dist = sum(graph.weights[p.src, p.dst] for p in path)

    linestring_path::Vector{LineString} = (
        s -> LineString([idx_to_point[s.src], idx_to_point[s.dst]])
    ).(path)

    return dist, linestring_path
end

"""
    point_in_convexhull(
        point::Point{2, Float64},
        exclusions::DataFrame
    )::Int

Check if a point is within a convex hull of exclusion zones.

# Arguments
- `point`: Point to check.
- `exclusions`: Exclusion zones.

# Returns
- Index of exclusion zone if point is within a convex hull, 0 otherwise.
"""
function point_in_convexhull(
    point::Point{2, Float64},
    exclusions::DataFrame
)::Int
    point_ag = AG.createpoint(point[1], point[2])
    convex_exclusions_ag = AG.convexhull.(exclusions.geometry)

    point_in_exclusion_zone = AG.contains.(convex_exclusions_ag, [point_ag])

    # If point is within an exclusion zone, add all polygon vertices to graph
    # ? Consider cases where point is within the convex hull of multiple polygons
    # findall(final_point_in_exclusion_zone)
    exclusion_index = findfirst(point_in_exclusion_zone)

    return isnothing(exclusion_index) ? 0 : exclusion_index
end

"""
    build_network!(
        points_from::Vector{Point{2, Float64}},
        points_to::Vector{Point{2, Float64}},
        exclusion_idxs::Vector{Int},
        current_point::Point{2, Float64},
        final_point::Point{2, Float64},
        exclusions::DataFrame,
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
    points_from::Vector{Point{2, Float64}},
    points_to::Vector{Point{2, Float64}},
    exclusion_idxs::Vector{Int},
    current_point::Point{2, Float64},
    final_point::Point{2, Float64},
    exclusions::DataFrame,
    final_exclusion_idx::Int
)::Nothing
    if current_point == final_point
        return
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

    for (vertex, next_exclusion_idx) in zip(candidates, next_exclusion_idxs)

        if isnothing(vertex) || vertex == current_point
            continue
        end

        # If edge already exists in network, skip
        if (current_point, vertex) ∈ zip(points_from, points_to)
            continue
        end

        # Record new point/edge
        push!(points_from, current_point)
        push!(points_to, vertex)
        push!(exclusion_idxs, next_exclusion_idx)

        # If verticex already visited/explored, skip
        # If final exclusion zone reached, vertices added to graph later in build_graph()
        if vertex ∈ points_from ||
            (next_exclusion_idx == final_exclusion_idx && !isnothing(next_exclusion_idx))
            continue
        end

        # Continue building network from this vertex to final point
        build_network!(
            points_from,
            points_to,
            exclusion_idxs,
            vertex,
            final_point,
            exclusions,
            final_exclusion_idx
        )
    end
    # TODO: check if path from current_point to next/intermediate point is feasible (no intersecting polygons in between)
end

"""
    build_graph(
        points_from::Vector{Point{2, Float64}},
        points_to::Vector{Point{2, Float64}},
        exclusions::DataFrame,
        final_polygon::AG.IGeometry{AG.wkbPolygon},
        initial_point::Point{2, Float64},
        final_point::Point{2, Float64}
    )::Tuple{SimpleWeightedGraph{Int64, Float64}, Vector{Point{2, Float64}}, Int64, Int64}

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
    points_from::Vector{Point{2, Float64}},
    points_to::Vector{Point{2, Float64}},
    exclusions::DataFrame,
    final_polygon::AG.IGeometry{AG.wkbPolygon},
    initial_point::Point{2, Float64},
    final_point::Point{2, Float64}
)::Tuple{SimpleWeightedGraph{Int64, Float64}, Vector{Point{2, Float64}}, Int64, Int64}
    is_polygon = !AG.isempty(final_polygon)

    poly_vertices::Vector{Point{2, Float64}} = is_polygon ?
        collect_polygon_vertices(final_polygon) :
        Point{2, Float64}[]

    # Collect connected points from network
    #? Use `Set`?
    chain_points::Vector{Point{2, Float64}} = [Point{2,Float64}(p[1], p[2])
        for p in unique(vcat(
            points_from, points_to, [final_point]
        ))
    ]

    unique_points::Vector{Point{2, Float64}} = unique(vcat(chain_points, poly_vertices))
    n_points::Int = length(unique_points)

    # Create the graph with one vertex per unique point.
    graph = SimpleWeightedGraph(n_points)

    # Map points to indices
    pt_to_idx = Dict{Point{2, Float64}, Int}(unique_points .=> 1:n_points)

    # Add candidate (chain) edges.
    add_edge!.(Ref(graph),
        map(p -> pt_to_idx[p], points_from),
        map(q -> pt_to_idx[q], points_to),
        haversine.(points_from, points_to)
    )

    # Only add polygon edges and extra visible connections if:
    # 1. A polygon is provided, AND
    # 2. Direct line/path from initial_point -> final_point is NOT visible (i.e. obstructed)
    if is_polygon
        if !is_visible(initial_point, final_point, final_polygon)

            # Add edges connecting any visible polygon vertices to other polyogn vertices
            for i in poly_vertices #i_vertices
                visibility_mask = (poly_vertices .!= i) .&
                    is_visible.(Ref(i), poly_vertices, Ref(exclusions))

                visible_points = poly_vertices[visibility_mask]

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

            # Connect all visible points to polygon vertices
            for i in poly_vertices
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
        if poly_vertices[1] != poly_vertices[end]
            add_edge!(
                graph,
                pt_to_idx[poly_vertices[end]],
                pt_to_idx[poly_vertices[1]],
                    haversine(poly_vertices[end], poly_vertices[1])
                )
            end
        else
            # Connect initial point to final point if visible
            add_edge!(
                graph,
                pt_to_idx[initial_point],
                pt_to_idx[final_point],
                haversine(initial_point, final_point)
            )
        end
    end

    return graph, unique_points, pt_to_idx[initial_point], pt_to_idx[final_point]
end

"""
    collect_polygon_vertices(polygon::AG.IGeometry{AG.wkbPolygon}
    )::Vector{Point{2, Float64}}

Collect all vertices of a polygon.

# Arguments
- `polygon`: Polygon to collect vertices from.

# Returns
- Vector of polygon vertices.
"""
function collect_polygon_vertices(polygon::AG.IGeometry{AG.wkbPolygon}
)::Vector{Point{2, Float64}}
    exterior_ring = AG.getgeom(polygon, 0)
    n_pts = AG.ngeom(exterior_ring)
    pts = Vector{Point{2,Float64}}(undef, n_pts)

    # iterate safely within all points in the exterior ring to populate pts with vertices
    @inbounds for i in 1:n_pts
        pt = AG.getpoint(exterior_ring, i - 1) # 0-indexed
        pts[i] = Point{2,Float64}(pt[1], pt[2])
    end

    # convex_hull = AG.convexhull(polygon)
    # #! Select points that are either outside or touching the convex hull
    # points_ag = AG.createpoint.(points_x, points_y)
    # points_outside_convex_hull = .!AG.contains.([convex_hull], points_ag)
    # points_touching_convex_hull = AG.touches.([convex_hull], points_ag)
    # target_points_mask::BitVector = points_outside_convex_hull .|| points_touching_convex_hull
    # poly_vertices = Point{2, Float64}.(points_x[target_points_mask], points_y[target_points_mask])
    # n_pts = length(poly_vertices)

    return unique(pts)
end
