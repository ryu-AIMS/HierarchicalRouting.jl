
"""
    get_feasible_matrix(nodes::Vector{Point{2, Float64}}, exclusions::DataFrame)

Create a matrix of distances of feasible paths between waypoints accounting for (avoiding) environmental constraints.

# Arguments
- `nodes::Vector{Point{2, Float64}}` : Vector of lat long tuples.
- `exclusions::DataFrame` : DataFrame containing exclusion zones representing given vehicle's cumulative environmental constraints.

# Returns
- `dist_matrix::Matrix{Float64}` : A matrix of distances between waypoints.
- `path_matrix` : A vector of tuples containing the graph, point to index mapping, and edges for each pair of waypoints.
"""
function get_feasible_matrix(nodes::Vector{Point{2, Float64}}, exclusions::DataFrame)
    n_points = length(nodes)
    dist_matrix = zeros(Float64, n_points, n_points)
    path_matrix = fill(Vector{LineString{2, Float64}}(), n_points, n_points)

    for j in 1:n_points
        for i in 1:j-1
            if nodes[i] != nodes[j]
                # TODO: Process elsewhere
                # Check if any of the points are within an exclusion zone
                if any(point_in_exclusion.([nodes[i], nodes[j]], [exclusions]))
                    dist_matrix[i, j] = dist_matrix[j, i] = Inf
                else
                    dist_matrix[i, j], path_matrix[i, j] = HierarchicalRouting.shortest_feasible_path(nodes[i], nodes[j], exclusions)
                    dist_matrix[j, i] = dist_matrix[i, j]
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
- `exclusions::DataFrame` : DataFrame containing exclusion zones representing given vehicle's cumulative environmental constraints.

# Returns
- `dist_vector::Vector{Float64}` : A vector of distances between waypoints.
- `path_vector` : A vector of paths between waypoints, represented as LineStrings.
"""
function get_feasible_vector(nodes::Vector{Point{2, Float64}}, exclusions::DataFrame)
    n_points = length(nodes)-1
    dist_vector = zeros(Float64, n_points)
    path_vector = fill(Vector{LineString{2, Float64}}(), n_points)

    for i in 1:n_points
            if nodes[i] != nodes[i+1]
                # TODO: Process elsewhere
                # Check if any of the points are within an exclusion zone
                if any(point_in_exclusion.([nodes[i], nodes[i+1]], [exclusions]))
                    dist_vector[i, i+1] = Inf
                else
                    dist_vector[i], path_vector[i] = HierarchicalRouting.shortest_feasible_path(nodes[i], nodes[i+1], exclusions)
                end
            end
    end

    return dist_vector, path_vector
end

function point_in_exclusion(point::Point{2,Float64}, exclusions::DataFrame)
    point_ag = AG.createpoint(point[1], point[2])
    return any(AG.contains.(exclusions.geometry, [point_ag]))
end

"""
    shortest_feasible_path(initial_point::Point{2, Float64}, final_point::Point{2, Float64}, exclusions::DataFrame)

Find the shortest feasible path between two points.
Use A* between all vertices on polygons that intersect with straight line to finish, from start pt and any other intersecting polygons.

# Arguments
- `initial_point::Point{2, Float64}`: Starting point of path.
- `final_point::Point{2, Float64}`: Ending point of path.
- `exclusions::DataFrame`: A DataFrame containing exclusion zone polygons.

# Returns
- `dist::Float64`: The distance of the shortest feasible path.
- `path::Vector{SimpleWeightedGraph{Int64, Float64}.Edge}`: The shortest feasible path as a vector of edges.
"""
function shortest_feasible_path(initial_point::Point{2, Float64}, final_point::Point{2, Float64}, exclusions::DataFrame,
    )

    # ! Under construction
    final_exclusion_idx = HierarchicalRouting.point_in_convexhull(final_point, exclusions)
    # initial_exclusion_idx = HierarchicalRouting.point_in_convexhull(initial_point, exclusions)

    # # ! Workaround for now
    # # If initial point is within an exclusion zone, reverse route and add all vertices of 'final' exclusion polygon to graph
    # if !iszero(initial_point_in_exclusion_zone)
    #     final_exclusion_idx = initial_exclusion_idx
    #     initial_point, final_point = final_point, initial_point
    # end

    points_from = Point{2,Float64}[]
    points_to = Point{2,Float64}[]
    exclusion_idx = Int[]

    # TODO: check if path from current_point to (left_point, right_point) is feasible (no intersecting polygons in between)
    build_network!(
        points_from,
        points_to,
        exclusion_idx,
        initial_point,
        final_point,
        exclusions,
        0,
        final_exclusion_idx
    )

    graph, idx_to_point, final_point_idx = build_graph(
        points_from,
        points_to,
        iszero(final_exclusion_idx) ? nothing : exclusions[final_exclusion_idx,:geometry],
        final_point
    )

    path = a_star(graph, 1, final_point_idx, graph.weights)
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
    exclusion_idx::Vector{Union{Int,Nothing}},
    current_point::Point{2, Float64},
    final_point::Point{2, Float64},
    exclusions::DataFrame,
    current_exclusion::Union{Int,Nothing} = nothing,
    final_exclusion_idx::Union{Int,Nothing} = nothing
)

Build a network of points to connect to each other.
Vectors `points_from`, `points_to`, and `exclusion_idx` are modified in place.

# Returns
- `nothing`
"""
function build_network!(
    points_from::Vector{Point{2, Float64}},
    points_to::Vector{Point{2, Float64}},
    exclusion_idx::Vector{Int},
    current_point::Point{2, Float64},
    final_point::Point{2, Float64},
    exclusions::DataFrame,
    current_exclusion::Int,
    final_exclusion_idx::Int
)

    if current_point == final_point # || (current_exclusion == final_exclusion_idx && !isnothing(current_exclusion))
        return
    end

    candidates, next_exclusion_idx = HierarchicalRouting.find_widest_points(
        current_point,
        final_point,
        exclusions,
        current_exclusion
    )

    for vertex in candidates

        # Record new point/edge
        push!(points_from, current_point)
        push!(points_to, vertex)
        push!(exclusion_idx, next_exclusion_idx)

        if (next_exclusion_idx == final_exclusion_idx && !isnothing(next_exclusion_idx))
            continue
        end
        # Skip vertices already visited
        if vertex ∈ points_from #|| isnothing(vertex) || vertex == current_point
            continue
        end

        # Continue building network from this vertex to final point
        build_network!(
            points_from,
            points_to,
            exclusion_idx,
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
        polygon::Union{Nothing, AG.IGeometry},
        final_point::Point{2, Float64}
    )

Build a graph network from a set of points.
If polygon is provided, connect it to visible points in graph and visible vertices.

# Arguments
- `points_from::Vector{Point{2, Float64}}`: Vector of points to connect from.
- `points_to::Vector{Point{2, Float64}}`: Vector of points to connect to.
- `polygon::Union{Nothing, AG.IGeometry}`: Polygon to connect to.
- `final_point::Point{2, Float64}`: Final point to connect to.

# Returns
- `graph::SimpleWeightedGraph`: Graph object.
- `idx_to_point::Vector{Point{2, Float64}}`: Mapping of indices to points.
- `final_pt_idx::Int`: Index of final point. Used in A* to specify target point.
"""
function build_graph(
    points_from::Vector{Point{2, Float64}},
    points_to::Vector{Point{2, Float64}},
    polygon::Union{Nothing, AG.IGeometry},
    final_point::Point{2, Float64}
)
    # Collect connected points from network
    chain_points = [Point{2,Float64}(p[1], p[2]) for p in vcat(points_from, points_to, [final_point])]

    is_polygon = !isnothing(polygon)
    poly_vertices = Point{2, Float64}[]

    if is_polygon
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

        poly_vertices = Point{2, Float64}.(points_x, points_y)
    end

    all_points = vcat(chain_points, poly_vertices)
    unique_points = unique(all_points)
    n_vertices = length(unique_points)

    # Create the graph with one vertex per unique point.
    graph = SimpleWeightedGraph(n_vertices)

    # Map points to indices
    pt_to_idx = Dict{Point{2, Float64}, Int}(unique_points .=> collect(1:n_vertices))

    # Add candidate (chain) edges.
    for (p,q) in zip(points_from, points_to)
        add_edge!(graph, pt_to_idx[p], pt_to_idx[q], euclidean(p, q))
    end

    # Add polygon edges if a polygon was provided.
    if is_polygon # && !isempty(poly_vertices)
        n_poly = length(poly_vertices)
        # Connect adjacent polygon vertices, and wrap around.
        for i in 1:n_poly
            p = poly_vertices[i]
            q = poly_vertices[mod1(i + 1, n_poly)]

            if is_visible(p, q, polygon)
                i_idx = pt_to_idx[p]
                j_idx = pt_to_idx[q]
                add_edge!(graph, i_idx, j_idx, euclidean(p, q))
            end
        end
        # Connect each polygon vertex to the final point if visible.
        for p in poly_vertices
            if is_visible(p, final_point, polygon)
                i_idx = pt_to_idx[p]
                j_idx = pt_to_idx[final_point]
                add_edge!(graph, i_idx, j_idx, euclidean(p, final_point))
            end
        end
    end

    return graph, unique_points, pt_to_idx[final_point]
end
