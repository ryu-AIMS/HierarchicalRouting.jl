
"""
    get_feasible_matrix(nodes::Vector{Point{2, Float64}}, exclusions::DataFrame)

Create a matrix of distances of feasible paths between waypoints accounting for (avoiding) environmental constraints.

# Arguments
- `nodes::Vector{Point{2, Float64}}` : Vector of lat long tuples.
- `exclusions::DataFrame` : DataFrame containing exclusion zones representing given vehicle's cumulative environmental constraints.
- `ignore_exclusions_flag::Bool` : Flag to ignore exclusions. Default is `true`.

# Returns
- `feasible_matrix::Matrix{Float64}` : A matrix of distances between waypoints.
- `feasible_path` : A vector of tuples containing the graph, point to index mapping, and edges for each pair of waypoints.
"""
function get_feasible_matrix(nodes::Vector{Point{2, Float64}}, exclusions::DataFrame,
    ignore_exclusions_flag::Bool = true
    )
    n_points = length(nodes)
    feasible_matrix = zeros(Float64, n_points, n_points)
    feasible_path = fill((Dict{Int64, Point{2, Float64}}(), Vector{SimpleWeightedGraphs.SimpleWeightedEdge{Int64, Float64}}()), n_points, n_points)

    for j in 1:n_points
        for i in 1:j-1
            if nodes[i] != nodes[j]
                # TODO: Process elsewhere
                # Check if any of the points are within an exclusion zone
                if any(
                    AG.contains.(
                        exclusions.geometry, Ref(AG.createpoint(nodes[j][1], nodes[j][2])))
                ) ||
                any(
                    AG.contains.(
                        exclusions.geometry, Ref(AG.createpoint(nodes[i][1], nodes[i][2])))
                )
                    feasible_matrix[i, j] = feasible_matrix[j, i] = Inf
                else
                    feasible_matrix[i, j], feasible_path[i,j] = HierarchicalRouting.shortest_feasible_path(nodes[i], nodes[j], exclusions, ignore_exclusions_flag)
                    feasible_matrix[j, i] = feasible_matrix[i, j]
                end
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
- `ignore_exclusions_flag::Bool`: Flag to ignore exclusions.

# Returns
- `dist::Float64`: The distance of the shortest feasible path.
- `path::Vector{SimpleWeightedGraph{Int64, Float64}.Edge}`: The shortest feasible path as a vector of edges.
"""
function shortest_feasible_path(initial_point::Point{2, Float64}, final_point::Point{2, Float64}, exclusions::DataFrame,
    ignore_exclusions_flag::Bool
    )

    final_exclusion_idx = nothing
    if ignore_exclusions_flag
        for (i, exclusion) in enumerate(eachrow(exclusions))
            # If final point is within an exclusion zone, add all vertices of 'final' exclusion polygon to graph
            if AG.contains(AG.convexhull(exclusion.geometry), AG.createpoint(final_point[1], final_point[2]))
                final_exclusion_idx = i
                break
            end

            # If final point is within an exclusion zone, reverse route and add all vertices of 'final' exclusion polygon to graph
            if AG.contains(AG.convexhull(exclusion.geometry), AG.createpoint(initial_point[1], initial_point[2]))
                final_exclusion_idx = i
                temp_point = initial_point
                initial_point = final_point
                final_point = temp_point
                points = [initial_point]
                break
            end
        end
    end

    points = [initial_point]
    parent_points = [initial_point]
    current_node_idx = 1
    exclusion_idx = [0]

    """
        process_point(
        candidate_point::Point{2, Float64},
        target_point::Union{Nothing, Point{2, Float64}},
        exclusions::DataFrame,
        ignore_exclusion_indices::Vector{Int}
        )

    Check if the line between two points is feasible, and add intermediate points to the path.

    # Arguments
    - `candidate_point::Point{2, Float64}`: The start of the line.
    - `target_point::Union{Nothing, Point{2, Float64}}`: The end of the line.
    - `exclusions::DataFrame`: The dataframe containing the polygon exclusions.
    - `ignore_exclusion_indices::Vector{Int}`: The indices of the exclusion polygons to disregard for crossings - the polygons associated with candidate and target points.
    """
    function process_point(
        candidate_point::Point{2, Float64},
        target_point::Union{Nothing, Point{2, Float64}},
        exclusions::DataFrame,
        ignore_exclusion_indices::Vector{Int}
    )
        if target_point !== nothing
            new_vertices, new_exclusion_idx = HierarchicalRouting.find_widest_points(candidate_point, target_point, exclusions, ignore_exclusions_flag ? ignore_exclusion_indices : [0])

            if target_point in new_vertices
                push!(parent_points, candidate_point)
                push!(points, target_point)
                push!(exclusion_idx, ignore_exclusion_indices[1])
                return
            end

            for intermediate_point in new_vertices
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
        if current_point == final_point || current_point in points[1:current_node_idx - 1] || exclusion_idx[current_node_idx] == final_exclusion_idx
            current_node_idx += 1
            continue
        end

        # Find edge points of next intersecting exclusion polygon, and polygon index
        vertices, current_exclusion_idx = HierarchicalRouting.find_widest_points(current_point, final_point, exclusions, ignore_exclusions_flag ? [exclusion_idx[current_node_idx]] : [0])

        # TODO: check if path from current_point to (left_point, right_point) is feasible (no intersecting polygons in between)
        process_point.(
            Ref(current_point),
            vertices,
            Ref(exclusions),
            Ref([current_exclusion_idx, exclusion_idx[current_node_idx]])
        )

        current_node_idx += 1
    end

    g, point_to_idx, idx_to_point = build_graph(
        points,
        parent_points,
        isnothing(final_exclusion_idx) ? nothing : exclusions[final_exclusion_idx,:geometry],
        final_point
    )

    path = a_star(g, 1, point_to_idx[final_point], g.weights)
    dist = sum(g.weights[p.src, p.dst] for p in path)

    return dist, (idx_to_point, path)
end

"""
    build_graph(pts::Vector{Point{2, Float64}}, exclusions::DataFrame)::SimpleWeightedGraph{Int64, Float64}

Construct a simple weighted graph between given points that do not intersect exclusions.

# Arguments
- `points::Vector{Point{2, Float64}`: A vector of points respresenting ordered end points.
- `parent_points::Vector{Point{2, Float64}`: A vector of points representing ordered start points.
- `polygon` : Polygon provided if it's convex hull contains final point.
    Include all vertices of exclusion polygon in graph.
- `final_point` : Include final point in graph.
    Provided if final point is within the convex hull of an exclusion zone.

# Returns
- `g::SimpleWeightedGraph{Int64, Float64}`: A simple weighted graph of all points/edges.
- `point_to_idx::Dict{Point{2, Float64}, Int64}`: A dictionary mapping points to indices.
- `idx_to_point::Dict{Int64, Point{2, Float64}`: A dictionary mapping indices to points.
"""
function build_graph(
    points::Vector{Point{2, Float64}},
    parent_points::Vector{Point{2, Float64}},
    polygon,
    final_point
    )
    # Dictionaries to map unique points and their indices
    point_to_idx = Dict{Point{2, Float64}, Int64}()
    idx_to_point = Dict{Int64, Point{2, Float64}}()
    idx_counter = 0

    for pt in points
        if !haskey(point_to_idx, pt)
            idx_counter += 1
            point_to_idx[pt] = idx_counter
            idx_to_point[idx_counter] = pt
        end
    end

    poly_vertices = []
    n_pts = 0
    if !isnothing(polygon)
        # Add all polygon vertices/edges to graph
        exterior_ring = AG.getgeom(polygon, 0)
        n_pts = AG.ngeom(exterior_ring)
        poly_vertices = [
            Point{2,Float64}(
                AG.getpoint(exterior_ring, i)[1],
                AG.getpoint(exterior_ring, i)[2]
            ) for i in 0:n_pts-1
        ]

        for v in poly_vertices
            if !haskey(point_to_idx, v)
                idx_counter += 1
                point_to_idx[v] = idx_counter
                idx_to_point[idx_counter] = v
            end
        end

        # Ensure final_point is in dictionaries
        if !haskey(point_to_idx, final_point)
            idx_counter += 1
            point_to_idx[final_point] = idx_counter
            idx_to_point[idx_counter] = final_point
        end
    end

    g = SimpleWeightedGraph(idx_counter)

    # Add edges between points & parents (points[1] has no parent)
    for i in 2:length(points)
        pt_i = points[i]
        parent_pt = parent_points[i]

        idx_pt = point_to_idx[pt_i]
        idx_parent = point_to_idx[parent_pt]

        add_edge!(g, idx_parent, idx_pt, euclidean(pt_i, parent_pt)) # haversine
    end

    # If `polygon`: Add edges between polygon vertices and final point if is_visible
    if !isnothing(polygon)
        final_pt_idx = point_to_idx[final_point]
        for i in 1:n_pts
            pt_i = poly_vertices[i]
            idx_i = point_to_idx[pt_i]

            # Connect adjacent vertex (wrapping around)
            j = (i % n_pts) + 1

            pt_j = poly_vertices[j]
            idx_j = point_to_idx[pt_j]

            add_edge!(g, idx_i, idx_j, euclidean(pt_i, pt_j))

            # Add edge between polygon vertex and final_point if is_visible
            if is_visible(pt_i, final_point, polygon)
                add_edge!(g, idx_i, final_pt_idx, euclidean(pt_i, final_point))
            end
        end
    end

    return g, point_to_idx, idx_to_point
end
