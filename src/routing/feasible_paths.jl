
"""
    get_feasible_matrix(
        nodes::Vector{Point{2, Float64}},
        exclusions::DataFrame
    )::Tuple{Matrix{Float64}, Matrix{Vector{LineString{2, Float64}}}}

Create matrices of distances and paths for feasible routes between waypoints accounting for
    (avoiding) environmental constraints.

# Arguments
- `nodes::Vector{Point{2, Float64}}` : Vector of lat long tuples.
- `exclusions::DataFrame` :Exclusion zones representing vehicle's environmental constraints.

# Returns
- `dist_matrix::Matrix{Float64}` : A matrix of distances between waypoints.
- `path_matrix::Matrix{Vector{LineString{2, Float64}}}` : A vector of paths between
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
                            HierarchicalRouting.shortest_feasible_path(
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
    get_feasible_vector(nodes::Vector{Point{2, Float64}}, exclusions::DataFrame)

Create vectors of:
- distances, and
- feasible paths

between sequential nodes, avoiding exclusions.

# Arguments
- `nodes::Vector{Point{2, Float64}}` : Vector of lat long tuples.
- `exclusions::DataFrame` : DataFrame containing exclusion zones representing given
vehicle's cumulative environmental constraints.

# Returns
- `dist_vector::Vector{Float64}` : A vector of distances between waypoints.
- `path_vector` : A vector of paths between waypoints, represented as LineStrings.
"""
function get_feasible_vector(nodes::Vector{Point{2, Float64}}, exclusions::DataFrame)
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
                HierarchicalRouting.shortest_feasible_path(
                    point_nodes[1], point_nodes[2], exclusions
                )
            end
        end
    end

    return dist_vector, path_vector
end

"""
    point_in_exclusion(point::Point{2, Float64}, exclusions::DataFrame)
    point_in_exclusion(point::Point{2, Float64}, exclusion::AG.IGeometry)

Check if a point is within an exclusion zone.

# Arguments
- `point::Point{2, Float64}`: Point to check.
- `exclusions::DataFrame`: A DataFrame containing exclusion zone polygons.
- `exclusion::AG.IGeometry`: One exclusion zone polygon from exclusions DataFrame.

# Returns
- `true` if point is within an exclusion zone, `false` otherwise.
"""
function point_in_exclusion(point::Point{2,Float64}, exclusions::DataFrame)
    point_ag = AG.createpoint(point[1], point[2])
    return any(AG.contains.(exclusions.geometry, [point_ag]))
end
function point_in_exclusion(point::Point{2,Float64}, exclusion::AG.IGeometry)
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
- `initial_point::Point{2, Float64}`: Starting point of path.
- `final_point::Point{2, Float64}`: Ending point of path.
- `exclusions::DataFrame`: A DataFrame containing exclusion zone polygons.

# Returns
- The distance of the shortest feasible path.
- The shortest feasible path as a vector of LineStrings.
"""
function shortest_feasible_path(
    initial_point::Point{2, Float64},
    final_point::Point{2, Float64},
    exclusions::DataFrame,
)::Tuple{Float64, Vector{LineString{2, Float64}}}
    final_exclusion_idx = HierarchicalRouting.point_in_convexhull(final_point, exclusions)
    initial_exclusion_idx = HierarchicalRouting.point_in_convexhull(initial_point, exclusions)

    points_from = Point{2,Float64}[]
    points_to = Point{2,Float64}[]
    exclusion_idxs = Int[]

    if iszero(initial_exclusion_idx)
        build_network!(
            points_from,
            points_to,
            exclusion_idxs,
            initial_point,
            final_point,
            exclusions,
            0,
            final_exclusion_idx
        )
    else
        initial_polygon = exclusions[initial_exclusion_idx, :geometry]
        poly_vertices = collect_polygon_vertices(initial_polygon)
        visible_vertices = poly_vertices[is_visible.([initial_point], poly_vertices, [exclusions])]

        # Connect each polygon vertex to the initial point if visible.
        n_verts = length(visible_vertices)
        append!(points_from, fill(initial_point, n_verts))
        append!(points_to, visible_vertices)
        append!(exclusion_idxs, fill(initial_exclusion_idx, n_verts))

        # Add all polygon vertices to graph
        n_vertices = length(poly_vertices)
        q = poly_vertices[mod1.((1:n_vertices) .+ 1, n_vertices)]

        append!(points_from, poly_vertices)
        append!(points_to, q)
        append!(exclusion_idxs, fill(initial_exclusion_idx, n_vertices))

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
            [initial_exclusion_idx],
            [final_exclusion_idx]
        )
    end

    graph, idx_to_point, initial_point_idx, final_point_idx = build_graph(
        points_from,
        points_to,
        exclusions,
        iszero(final_exclusion_idx) ? nothing : exclusions[final_exclusion_idx,:geometry],
        initial_point,
        final_point
    )

    path = a_star(graph, initial_point_idx, final_point_idx, graph.weights)
    dist = sum(graph.weights[p.src, p.dst] for p in path)

    linestring_path = [LineString([idx_to_point[segment.src], idx_to_point[segment.dst]]) for segment in path]

    return dist, linestring_path
end

"""
    point_in_convexhull(
        point::Point{2, Float64},
        exclusions::DataFrame
    )::Int

Check if a point is within a convex hull of exclusion zones.

# Arguments
- `point::Point{2, Float64}`: Point to check.
- `exclusions::DataFrame`: Exclusion zones.

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
        current_exclusion_idx::Int,
        final_exclusion_idx::Int
    )

Build a network of points to connect to each other.
Vectors `points_from`, `points_to`, and `exclusion_idxs` are modified in place.
"""
function build_network!(
    points_from::Vector{Point{2, Float64}},
    points_to::Vector{Point{2, Float64}},
    exclusion_idxs::Vector{Int},
    current_point::Point{2, Float64},
    final_point::Point{2, Float64},
    exclusions::DataFrame,
    current_exclusion_idx::Int,
    final_exclusion_idx::Int
)

    if current_point == final_point # || (current_exclusion_idx == final_exclusion_idx && !isnothing(current_exclusion_idx))
        return
    end

    candidates, next_exclusion_idxs = HierarchicalRouting.find_widest_points(
        current_point,
        final_point,
        exclusions
    )

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
            next_exclusion_idx,
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
        final_polygon::Union{Nothing, AG.IGeometry},
        initial_point::Point{2, Float64},
        final_point::Point{2, Float64}
    )

Build a graph network from a set of points.
If polygon is provided, connect it to visible points in graph and visible vertices.

# Arguments
- `points_from::Vector{Point{2, Float64}}`: Vector of points to connect from.
- `points_to::Vector{Point{2, Float64}}`: Vector of points to connect to.
- `exclusions::DataFrame`: DataFrame containing exclusion zones.
- `polygon::Union{Nothing, AG.IGeometry}`: Polygon to connect to.
- `initial_point::Point{2, Float64}`: Initial point to connect from.
- `final_point::Point{2, Float64}`: Final point to connect to.

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
    final_polygon::Union{Nothing, AG.IGeometry},
    initial_point::Point{2, Float64},
    final_point::Point{2, Float64}
)
    is_polygon = !isnothing(final_polygon)

    poly_vertices = is_polygon ?
        collect_polygon_vertices(final_polygon) :
        Point{2, Float64}[]

    # Collect connected points from network
    chain_points = [Point{2,Float64}(p[1], p[2])
        for p in unique(vcat(
            points_from, points_to, [final_point]
        ))
    ]
    all_points = vcat(chain_points, poly_vertices)
    unique_points = unique(all_points)
    n_points = length(unique_points)

    # Create the graph with one vertex per unique point.
    graph = SimpleWeightedGraph(n_points)

    # Map points to indices
    pt_to_idx = Dict{Point{2, Float64}, Int}(unique_points .=> collect(1:n_points))

    # Add candidate (chain) edges.
    for (p,q) in zip(points_from, points_to)
        add_edge!(graph, pt_to_idx[p], pt_to_idx[q], euclidean(p, q))
    end

    # Add polygon edges if a polygon was provided.
    if is_polygon
        n_vertices = length(poly_vertices)
        for i in 1:n_vertices
            p = poly_vertices[i]
            q = poly_vertices[mod1(i + 1, n_vertices)]

            # Connect adjacent polygon vertices, and wrap around.
            i_idx = pt_to_idx[p]
            j_idx = pt_to_idx[q]
            add_edge!(graph, i_idx, j_idx, euclidean(p, q))

            # connect to visible points in the graph
            for point in chain_points
                if is_visible(p, point, exclusions)
                    j_idx = pt_to_idx[point]
                    add_edge!(graph, i_idx, j_idx, euclidean(p, point))
                end
            end
        end
    end

    return graph, unique_points, pt_to_idx[initial_point], pt_to_idx[final_point]
end

function collect_polygon_vertices(polygon::AG.IGeometry)
    exterior_ring = AG.getgeom(polygon, 0)
    n_pts = AG.ngeom(exterior_ring)

    vertices = AG.getpoint.([exterior_ring], 0:n_pts-1)
    points_x, points_y = getindex.(vertices, 1), getindex.(vertices, 2)

    # convex_hull = AG.convexhull(polygon)
    # # Select points that are either outside or touching the convex hull
    # points_ag = AG.createpoint.(points_x, points_y)
    # points_outside_convex_hull = .!AG.contains.([convex_hull], points_ag)
    # points_touching_convex_hull = AG.touches.([convex_hull], points_ag)
    # target_points_mask::BitVector = points_outside_convex_hull .|| points_touching_convex_hull
    # poly_vertices = Point{2, Float64}.(points_x[target_points_mask], points_y[target_points_mask])
    # n_pts = length(poly_vertices)

    return Point{2, Float64}.(points_x, points_y)
end
