
"""
    get_feasible_matrix(nodes::Vector{Point{2, Float64}}, exclusions::DataFrame)

Create a matrix of distances of feasible paths between waypoints accounting for (avoiding) environmental constraints.

# Arguments
- `nodes::Vector{Point{2, Float64}}` : Vector of lat long tuples.
- `exclusions::DataFrame` : DataFrame containing exclusion zones representing given vehicle's cumulative environmental constraints.

# Returns
- `dist_matrix::Matrix{Float64}` : A matrix of distances between waypoints.
- `path_matrix` : A matrix of paths between waypoints, represented as LineStrings.
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
                if point_in_exclusion([nodes[i], nodes[j]], exclusions)
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

function point_in_exclusion(point::Point{2,Float64}, exclusions::DataFrame)
    point_ag = AG.createpoint(point[1], point[2])
    return any(AG.contains.(exclusions.geometry, [point_ag]))
end
function point_in_exclusion(points::Vector{Point{2,Float64}}, exclusions::DataFrame)
    return any(point_in_exclusion.(points, [exclusions]))
end

"""
    shortest_feasible_path(
    initial_point::Point{2, Float64},
    final_point::Point{2, Float64},
    exclusions::DataFrame
)

Find the shortest feasible path between two points.
Use A* between all vertices on polygons that intersect with straight line to finish, from start pt and any other intersecting polygons.

# Arguments
- `initial_point::Point{2, Float64}`: Starting point of path.
- `final_point::Point{2, Float64}`: Ending point of path.
- `exclusions::DataFrame`: A DataFrame containing exclusion zone polygons.

# Returns
- `dist::Float64`: The distance of the shortest feasible path.
- `linestring_path::Vector{LineString{2, Float64}}`: The path as a vector of LineStrings.
"""
function shortest_feasible_path(
    initial_point::Point{2, Float64},
    final_point::Point{2, Float64},
    exclusions::DataFrame,
)
    final_point_ag = AG.createpoint(final_point[1], final_point[2])
    initial_point_ag = AG.createpoint(initial_point[1], initial_point[2])

    convex_exclusions_ag = AG.convexhull.(exclusions.geometry)
    final_point_in_exclusion_zone = AG.contains.(convex_exclusions_ag, [final_point_ag])
    initial_point_in_exclusion_zone = AG.contains.(convex_exclusions_ag, [initial_point_ag])

    # If final point is within an exclusion zone, add all vertices of 'final' exclusion polygon to graph
    # ? Consider cases where point is within the convex hull of multiple polygons
    final_exclusion_idx = findfirst(final_point_in_exclusion_zone) # findall(final_point_in_exclusion_zone)
    initial_exclusion_idx = findfirst(initial_point_in_exclusion_zone) # findall(initial_point_in_exclusion_zone)

    exclusions_containing_points = unique(filter(x -> !isnothing(x), [initial_exclusion_idx, final_exclusion_idx]))

    # ! Workaround for now
    # If initial point is within an exclusion zone, reverse route and add all vertices of 'final' exclusion polygon to graph
    if any(initial_point_in_exclusion_zone)
        final_exclusion_idx = initial_exclusion_idx
        initial_point, final_point = final_point, initial_point
    end

    points_from = Point{2,Float64}[]
    points_to = Point{2,Float64}[]
    exclusion_idx = Union{Int,Nothing}[]

    # recursive function to build list of points to add to graph
    build_network!(
        points_from,
        points_to,
        exclusion_idx,
        initial_point,
        final_point,
        exclusions,
        nothing, #exclusions_containing_points, #
        final_exclusion_idx
    )

    graph, idx_to_pt, final_pt_idx = build_graph(
        points_from,
        points_to,
        isnothing(final_exclusion_idx) ? nothing : exclusions[final_exclusion_idx,:geometry],
        final_point
    )

    path = a_star(graph, 1, final_pt_idx, graph.weights)
    dist = sum(graph.weights[p.src, p.dst] for p in path)

    linestring_path = [LineString([idx_to_pt[segment.src], idx_to_pt[segment.dst]]) for segment in path]

    return dist, linestring_path
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
    exclusion_idx::Vector{Union{Int,Nothing}},
    current_point::Point{2, Float64},
    final_point::Point{2, Float64},
    exclusions::DataFrame,
    current_exclusion::Union{Int,Nothing} = nothing,
    final_exclusion_idx::Union{Int,Nothing} = nothing
)

    if current_point == final_point # || (current_exclusion == final_exclusion_idx && !isnothing(current_exclusion))
        return
    end

    candidates, next_exclusion_idx = HierarchicalRouting.find_widest_points(
        current_point,
        final_point,
        exclusions,
        isnothing(current_exclusion) ? [0] : [current_exclusion]
    )

    for vertex in candidates
        # Skip vertices already visited
        if vertex ∈ points_from
            continue
        end

        # Record new point/edge
        push!(points_from, current_point)
        push!(points_to, vertex)
        push!(exclusion_idx, next_exclusion_idx)

        if (next_exclusion_idx == final_exclusion_idx && !isnothing(next_exclusion_idx))
            continue
        end

        # Recursively process this candidate vertex.
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
