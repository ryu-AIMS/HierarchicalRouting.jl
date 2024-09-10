
"""
    get_feasible_matrix(waypoints::Vector{Point{2, Float32}}, exclusions::DataFrame)

Create a distance matrix between waypoints accounting for environmental constraints.

# Arguments
- `waypoint::Vector{Point{2, Float32}s` : Vector of lat long tuples.
- `exclusions::DataFrame` : DataFrame containing exclusion zones representing given vehicle's cumulative environmental constraints.

# Returns
Feasible distance matrix between waypoints.
"""
function get_feasible_matrix(waypoints::Vector{Point{2, Float32}}, exclusions::DataFrame)::Matrix{Float64}
    n_waypoints = length(waypoints) - 1
    feasible_matrix = zeros(Float64, n_waypoints, n_waypoints)

    for j in 1:n_waypoints
        for i in 1:j-1
            feasible_matrix[i, j] = shortest_feasible_path((waypoints[i], waypoints[j]), exclusions)[1]
            feasible_matrix[j, i] = feasible_matrix[i, j]
        end
    end

    return feasible_matrix
end

"""
    shortest_feasible_path(line_pts::Tuple{Point{2, Float32}, Point{2, Float32}}, exclusions::DataFrame)

Find the shortest feasible path between two points.
Use A* between all vertices on polygons that intersect with straight line to finish, from start pt and any other intersecting polygons.

# Arguments
- `line_pts::Tuple{Point{2, Float32}, Point{2, Float32}}`: A tuple containing two points on a line.
- `exclusions::DataFrame`: A DataFrame containing exclusions.

# Returns
- `dist::Float64`: The distance of the shortest feasible path.
- `path::Vector{SimpleWeightedGraph{Int64, Float64}.Edge}`: The shortest feasible path as a vector of edges.
"""
function shortest_feasible_path(line_pts::Tuple{Point{2, Float32}, Point{2, Float32}}, exclusions::DataFrame)
    pts = Set([line_pts[1]])

    new_pts = extract_unique_vertices(line_pts, exclusions)

    while new_pts !== nothing
        union!(pts, new_pts)

        extracted_pts = unique(vcat([extract_unique_vertices((p, line_pts[2]), exclusions) for p in new_pts]...))
        new_pts = [pt for pt in extracted_pts if !(pt in pts) && pt !== nothing]
    end

    pts = collect(pts)
    push!(pts, line_pts[2])

    g = build_graph(pts, exclusions)

    path = a_star(g, 1, length(pts), weights(g))
    dist = sum(g.weights[p.src, p.dst] for p in path)

    return dist, path
end

"""
    extract_unique_vertices(exclusions::DataFrame)::Vector{Point{2, Float32}}

Extracts unique vertices from the given DataFrame of exclusions.

# Arguments
- `line_pts::Tuple{Point{2, Float32}, Point{2, Float32}}`: A tuple containing start/emd points of a line.
- `exclusions::DataFrame`: The DataFrame containing exclusions.

# Returns
- `Vector{Point{2, Float32}}`: A vector of unique vertices.
"""
function extract_unique_vertices(line_pts::Tuple{Point{2, Float32}, Point{2, Float32}}, exclusions::DataFrame)::Union{Nothing, Vector{Point{2, Float32}}}
    unique_vertices = Set{Point{2, Float32}}()

    polygons = intersecting_polygons(line_pts, exclusions)

    for polygon in polygons

        exterior_ring = AG.getgeom(polygon, 0)
        n_pts = AG.ngeom(exterior_ring)

        for i in 0:n_pts - 1
            x, y, _ = AG.getpoint(exterior_ring, i)
            point = Point{2, Float32}(x, y)

            if !(point in unique_vertices)
                push!(unique_vertices, point)
            end
        end
    end

    if isempty(unique_vertices)
        return nothing
    else
    return collect(unique_vertices)
    end
end

"""
    intersecting_polygons(pt_a::Point{2, Float32}, pt_b::Point{2, Float32}, exclusions::DataFrame)

Find polygons that intersect with a line segment.

# Arguments
-
- `exclusions::DataFrame`: The dataframe containing the polygon exclusions.

# Returns
- `crossed_polygons`: A list of polygons that intersect with the line segment.
"""
function intersecting_polygons(line_pts::Tuple{Point{2, Float32}, Point{2, Float32}}, exclusions::DataFrame)
    crossed_polygons = []

    line = AG.createlinestring([[line_pts[1][1], line_pts[1][2]], [line_pts[2][1], line_pts[2][2]]])

    for (i, row) in enumerate(eachrow(exclusions))
        if AG.crosses(line, row.geometry)
            push!(crossed_polygons, row.geometry)
        end
    end

    return crossed_polygons
end

"""
    build_graph(pts::Vector{Point{2, Float32}}, exclusions::DataFrame)::SimpleWeightedGraph{Int64, Float64}

Construct a simple weighted graph between given points that do not intersect exclusions.

# Arguments
- `pts::Vector{Point{2, Float32}`: A vector of 2D points.
- `exclusions::DataFrame`: A DataFrame containing exclusions.

# Returns
Simple weighted graph with distances between points.
"""
function build_graph(pts::Vector{Point{2, Float32}}, exclusions::DataFrame)::SimpleWeightedGraph{Int64, Float64}
    g = SimpleWeightedGraph(length(pts))
    local_graphs = [SimpleWeightedGraph{Int64, Float64}(length(pts)) for _ in 1:nthreads()]

    @floop for j in 1:length(pts)
        local_g = local_graphs[threadid()]

        for i in 1:j-1
            if !intersects_polygon(pts[i], pts[j], exclusions)
                add_edge!(local_g, i, j, haversine(pts[i], pts[j]))
            end
        end
    end

    for local_g in local_graphs
        for e in edges(local_g)
            add_edge!(g, src(e), dst(e), weight(e))
        end
    end

    return g
end

"""
    intersects_polygon(line::LineString{2, Float32}, exclusions::DataFrame)::Bool

Check if a line intersects with polygons in a dataframe.

# Arguments
- `pt_a::Point{2, Float32}`: The start point of the line.
- `pt_b::Point{2, Float32}`: The end point of the line.
- `exclusions::DataFrame`: The dataframe containing the polygon exclusions.

# Returns
- `Bool`: `true` if the line intersects with any polygon in the dataframe, `false` otherwise.
"""
function intersects_polygon(pt_a::Point{2, Float32}, pt_b::Point{2, Float32}, exclusions::DataFrame)::Bool
    # TODO: Parallelise loop?
    for row in eachrow(exclusions)
        if AG.crosses(AG.createlinestring([[pt_a[1], pt_a[2]], [pt_b[1], pt_b[2]]]), row.geometry)
            return true
        end
    end
    return false
end
