
"""
    get_feasible_matrix(waypoints::Vector{Point{2, Float64}}, exclusions::DataFrame)

Create a matrix of distances of feasible paths between waypoints accounting for (avoiding) environmental constraints.

# Arguments
- `points::Vector{Point{2, Float64}` : Vector of lat long tuples.
- `exclusions::DataFrame` : DataFrame containing exclusion zones representing given vehicle's cumulative environmental constraints.

# Returns
- `feasible_matrix::Matrix{Float64}` : A matrix of distances between waypoints.
- `feasible_path` : A vector of tuples containing the graph, point to index mapping, and edges for each pair of waypoints.
"""
function get_feasible_matrix(points::Vector{Point{2, Float64}}, exclusions::DataFrame)
    n_points = length(points)
    feasible_matrix = zeros(Float64, n_points, n_points)
    feasible_path = fill((Dict{Int64, Point{2, Float64}}(), Vector{SimpleWeightedGraphs.SimpleWeightedEdge{Int64, Float64}}()), n_points, n_points)

    for j in 1:n_points
        for i in 1:j-1
            if points[i] != points[j]
                feasible_matrix[i, j], feasible_path[i,j] = shortest_feasible_path(points[i], points[j], exclusions)
                feasible_matrix[j, i] = feasible_matrix[i, j]
            end
        end
    end

    return feasible_matrix, feasible_path
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
function shortest_feasible_path(initial_point::Point{2, Float64}, final_point::Point{2, Float64}, exclusions::DataFrame)
    points = [initial_point]
    parent_points = [initial_point]
    current_node_idx = 1
    exclusion_idx = [0]

    function process_point(
        candidate_point::Point{2, Float64},
        target_point::Union{Nothing, Point{2, Float64}},
        exclusions::DataFrame,
        ignore_exclusion_indices::Vector{Int})

        if target_point !== nothing
            (left_point, right_point), new_exclusion_idx = HierarchicalRouting.find_next_points(candidate_point, target_point, exclusions, ignore_exclusion_indices)

            if left_point == target_point || right_point == target_point
                push!(parent_points, candidate_point)
                push!(points, target_point)
                push!(exclusion_idx, ignore_exclusion_indices[1])
                return
            end

            for intermediate_point in [left_point, right_point]
                if intermediate_point !== nothing && intermediate_point !== target_point
                    push!(parent_points, candidate_point)
                    push!(points, intermediate_point)
                    push!(exclusion_idx, new_exclusion_idx)
                end
            end
        end
    end

    while current_node_idx <= length(points)
        current_point = points[current_node_idx]
        if current_point == final_point || current_point in points[1:current_node_idx - 1]
            current_node_idx += 1
            continue
        end

        # Find edge points of next intersecting exclusion polygon, and polygon index
        (left_point, right_point), current_exclusion_idx = HierarchicalRouting.find_next_points(current_point, final_point, exclusions, [exclusion_idx[current_node_idx]])

        # check if current_point to (left_point, right_point) is feasible or if there are intersecting polygons in between
        map(
            point -> process_point(current_point, point, exclusions, [current_exclusion_idx, exclusion_idx[current_node_idx]]),
            filter(x -> x !== nothing, [left_point, right_point])
        )

        current_node_idx += 1
    end

    g, point_to_idx, idx_to_point = build_graph(points, parent_points)

    path = a_star(g, 1, point_to_idx[final_point], g.weights)
    dist = sum(g.weights[p.src, p.dst] for p in path)

    return dist, (idx_to_point, path)
end

"""
    find_next_points(current_point::Point{2, Float64}, final_point::Point{2, Float64}, exclusions::DataFrame, current_exclusions_idx::Vector{Int}) #::Union{Nothing, Tuple}

Find the widest vertices on each (left/right) side of the line to the final point.

# Arguments
- `current_point::Point{2, Float64}`: The current point marking start of line.
- `final_point::Point{2, Float64}`: The final point marking end of line.
- `exclusions::DataFrame`: The dataframe containing the polygon exclusions.
- `current_exclusions_idx::Vector{Int}`: The indices of the exclusion polygons to disregard for crossings.

# Returns
- `furthest_vert_L::Union{Nothing, Point{2, Float64}}`: The furthest point on the left side of the line.
- `furthest_vert_R::Union{Nothing, Point{2, Float64}}`: The furthest point on the right side of the line.
- `polygon_idx::Int`: The index of the polygon crossed.
"""
function find_next_points(current_point::Point{2, Float64}, final_point::Point{2, Float64}, exclusions::DataFrame, current_exclusions_idx::Vector{Int}) #::Union{Nothing, Tuple}
    max_dist_L, max_dist_R = -Inf32, -Inf32
    furthest_vert_L, furthest_vert_R = nothing, nothing # Set{Point{2, Float64}}()

    # TODO: Allow for multiple equidistant furthest points

    # Next polygon crossed by line to end point
    polygon, polygon_idx = closest_crossed_polygon(current_point, final_point, exclusions, current_exclusions_idx)
    if polygon === nothing
        return (final_point, nothing), 0
    end

    # Geometry of polygon and number of vertices
    exterior_ring = AG.getgeom(polygon, 0)
    n_pts = AG.ngeom(exterior_ring)

    # For each polygon vertex, find the furthest points on the left and right side of line
    for i in 0:n_pts - 1
        x, y, _ = AG.getpoint(exterior_ring, i)
        pt = Point{2, Float64}(x, y)

        # Perp dist of point to line
        dist = GO.distance(pt, LineString([current_point, final_point]))
        side = (final_point[1] - current_point[1]) * (pt[2] - current_point[2]) - (final_point[2] - current_point[2]) * (pt[1] - current_point[1])

        # Check side (L/R) and update furthest points
        if side > 0 && dist > max_dist_L
                max_dist_L = dist
                furthest_vert_L = pt #Set([pt])
        elseif side < 0 && dist > max_dist_R
                max_dist_R = dist
                furthest_vert_R = pt #Set([pt])
        end
    end

    return (furthest_vert_L, furthest_vert_R), polygon_idx
end

"""
    closest_crossed_polygon(current_point::Point{2, Float64}, final_point::Point{2, Float64}, exclusions::DataFrame, current_exclusions_idx::Vector{Int})

Find polygons that intersect with a line segment.

# Arguments
- `current_point::Point{2, Float64}`: The current point marking start of line.
- `final_point::Point{2, Float64}`: The final point marking end of line.
- `exclusions::DataFrame`: The dataframe containing the polygon exclusions.
- `current_exclusions_idx::Vector{Int}`: The indices of the exclusion polygons to disregard.

# Returns
- `closest_polygon`: The first/closest polygon that intersects with the line segment.
- `polygon_idx`: The index of the(first) polygon crossed.
"""
function closest_crossed_polygon(current_point::Point{2, Float64}, final_point::Point{2, Float64}, exclusions::DataFrame, current_exclusions_idx::Vector{Int})
    closest_polygon = nothing
    min_dist = Inf
    polygon_idx = nothing

    line = AG.createlinestring([current_point[1], final_point[1]], [current_point[2], final_point[2]])

    # Define bounding box for line to exclude polygons that do not intersect
    line_min_x, line_max_x = min(current_point[1], final_point[1]), max(current_point[1], final_point[1])
    line_min_y, line_max_y = min(current_point[2], final_point[2]), max(current_point[2], final_point[2])

    for (i, row) in enumerate(eachrow(exclusions))
        # Skip current polygon
        if i in current_exclusions_idx
            continue
        end

        # Check if any vertices of the polygon are inside the bounding box
        exterior_ring = AG.getgeom(row.geometry, 0)
        n_pts = AG.ngeom(exterior_ring)

        # Check if any vertices of the polygon are inside the bounding box of line
        has_vertex_in_bbox = any(
            line_min_x <= AG.getpoint(exterior_ring, j)[1] <= line_max_x &&
            line_min_y <= AG.getpoint(exterior_ring, j)[2] <= line_max_y
            for j in 0:n_pts - 1
        )

        crosses_bbox = false
        if !has_vertex_in_bbox
            # Check if line crosses bounding box
            poly_xs, poly_ys, _ = [AG.getpoint(exterior_ring, j) for j in 0:n_pts - 1]

            crosses_bbox = (
                line_min_x <= maximum(poly_xs) &&
                line_max_x >= minimum(poly_xs) &&
                line_min_y <= maximum(poly_ys) &&
                line_max_y >= minimum(poly_ys)
            )
        end

        # Skip polygons with no vertices in or crossing bounding box
        if !has_vertex_in_bbox && !crosses_bbox
            continue
        end

        if AG.crosses(line, row.geometry)
            dist = GO.distance(current_point, row.geometry)

            if dist < min_dist
                min_dist = dist
                closest_polygon = row.geometry
                polygon_idx = i
            end
        end
    end

    return closest_polygon, polygon_idx
end

"""
    build_graph(pts::Vector{Point{2, Float64}}, exclusions::DataFrame)::SimpleWeightedGraph{Int64, Float64}

Construct a simple weighted graph between given points that do not intersect exclusions.

# Arguments
- `points::Vector{Point{2, Float64}`: A vector of points respresenting ordered end points.
- `parent_points::Vector{Point{2, Float64}`: A vector of points representing ordered start points.

# Returns
Simple weighted graph with distances between points.
"""
function build_graph(points::Vector{Point{2, Float64}}, parent_points::Vector{Point{2, Float64}})#::SimpleWeightedGraph{Int64, Float64}
    # Dictionaries to map unique points and their indices
    point_to_idx = Dict{Point{2, Float64}, Int64}()
    idx_to_point = Dict{Int64, Point{2, Float64}}()
    idx_counter = 1

    for pt in points
        if !haskey(point_to_idx, pt)
            point_to_idx[pt] = idx_counter
            idx_to_point[idx_counter] = pt
            idx_counter += 1
        end
    end

    g = SimpleWeightedGraph(idx_counter - 1)

    # Add edges between points & parents (from 2 because first point has no parent)
    for i in 2:length(points)
        pt_i = points[i]
        parent_pt = parent_points[i]

        idx_pt = point_to_idx[pt_i]
        idx_parent = point_to_idx[parent_pt]

        add_edge!(g, idx_parent, idx_pt, euclidean(pt_i, parent_pt)) # haversine
    end

    return g, point_to_idx, idx_to_point
end
