
"""
    find_widest_points(
        current_point::Point{2, Float64},
        final_point::Point{2, Float64},
        exclusions::DataFrame,
        visited_exclusion_idxs::Vector{Int} = [0]
    )::Tuple{Vector{Point{2, Float64}}, Vector{Int}}

Find the intermediate vertices between two points by using the widest visible vertices on
each (left/right) side of the direct-line path by closest crossed polygon.
If these intermediate vertices are not visible, the function is called recursively to find
the widest visible vertices.

# Arguments
- `current_point`: The current point marking start of line.
- `final_point`: The final point marking end of line.
- `exclusions`: The dataframe containing the polygon exclusions.
- `visited_exclusion_idxs`: Already-crossed exclusion zone indices.

# Returns
- A vector of points representing the widest visible points on each side of the line.
- A vector of polygon indices crossed by the line.
"""
function find_widest_points(
    current_point::Point{2, Float64},
    final_point::Point{2, Float64},
    exclusions::DataFrame,
    visited_exclusion_idxs::Vector{Int} = [0]
)::Tuple{Vector{Point{2, Float64}}, Vector{Int}}

    if is_visible(current_point, final_point, exclusions)
        return [final_point], [0]
    end

    # Find the first polygon that the line (current_point, final_point) crosses,
    # ignoring any polygons already visited.
    (polygon, exterior_ring, n_pts), polygon_idx = closest_crossed_polygon(
        current_point, final_point, exclusions, visited_exclusion_idxs
    )

    if AG.isempty(polygon)
        return [final_point], [0]
    end

    # Search for the widest vertices (left and right) along the polygon's exterior.
    widest_vertices::NTuple{2, Union{Point{2, Float64}, Nothing}} = search_widest_points(
        current_point, final_point, exterior_ring, n_pts
    )

    candidates = Point{2,Float64}[]
    poly_indices = Int[]
    # For each vertex: record if visible, else recursively find widest visible vertices.
    for vertex::Point{2,Float64} ∈ filter(!isnothing, (widest_vertices))
        if is_visible(current_point, vertex, exclusions, visited_exclusion_idxs)
            push!(candidates, vertex)
            push!(poly_indices, polygon_idx)
        else
            new_pts, new_poly_idxs = find_widest_points(
                current_point,
                vertex,
                exclusions,
                union(visited_exclusion_idxs, [polygon_idx])
            )
            append!(candidates, new_pts)
            append!(poly_indices, new_poly_idxs)
        end
    end

    return candidates, poly_indices
end

"""
    vector_angle(
        base_vector::Vector{Float64}, candidate_vector::Vector{Float64}
    )::Float64
    vector_angle(
        base_vector::Vector{Float64}, candidate_tuple::Tuple{Float64, Float64}
    )::Float64

Calculate the signed angle (radians) of a candidate vector relative to a base vector.

# Arguments
- `base_vector`: The base vector.
- `candidate_vector`: The candidate vector.
- `candidate_tuple`: The candidate vector as a tuple.

# Returns
- The signed angle between the two vectors in radians.
"""
function vector_angle(
    base_vector::Vector{Float64}, candidate_tuple::Tuple{Float64, Float64}
)::Float64
    return vector_angle(base_vector, collect(candidate_tuple))
end
function vector_angle(
    base_vector::Vector{Float64}, candidate_vector::Vector{Float64}
)::Float64
    # Cross product tells us the direction of the angle
    cross = base_vector[1]*candidate_vector[2] - base_vector[2]*candidate_vector[1]
    # Dot product tells us the magnitude of the angle
    dot_val = base_vector[1]*candidate_vector[1] + base_vector[2]*candidate_vector[2]

    return atan(cross, dot_val)
end

"""
    search_widest_points(
        current_point::Point{2, Float64},
        final_point::Point{2, Float64},
        exterior_ring::AG.IGeometry{AG.wkbLineString},
        n_pts::Int32,
    )::NTuple{2, Union{Point{2, Float64}, Nothing}}

Search for the widest visible vertices on each side of the base vector.

# Arguments
- `current_point`: The current point marking start of line.
- `final_point`: The final point marking end of line.
- `exterior_ring`: The exterior ring of the polygon.
- `n_pts`: The number of vertices in the polygon.

# Returns
The widest visible polygon vertex on each (L/R) side of the base vector/line.
"""
function search_widest_points(
    current_point::Point{2, Float64},
    final_point::Point{2, Float64},
    exterior_ring::AG.IGeometry{AG.wkbLineString},
    n_pts::Int32,
)::NTuple{2, Union{Point{2, Float64}, Nothing}}
    max_angle_L, max_angle_R = 0.0, 0.0
    furthest_vert_L, furthest_vert_R = nothing, nothing
    # Base vector from current_point to final_point to calc angle from
    base_vector = [
        final_point[1] - current_point[1],
        final_point[2] - current_point[2]
    ]

    getpoint = AG.getpoint
    x::Float64, y::Float64 = 0.0, 0.0
    pt = Point{2, Float64}(x, y)
    candidate_tuple::Tuple{Float64, Float64} = (0.0, 0.0)
    angle::Float64 = 0.0
    for i in 0:n_pts - 1
        x, y, _ = getpoint(exterior_ring, i)
        pt = Point{2, Float64}(x, y)

        # Tuple from current_point to this vertex.
        candidate_tuple = ((pt[1] - current_point[1]), (pt[2] - current_point[2]))

        # Signed angle between base vector and vector to vertex
        angle = vector_angle(base_vector, candidate_tuple)

        if pt != current_point
            if angle > 0 && angle > max_angle_L
                max_angle_L = angle
                furthest_vert_L = pt
            elseif angle < 0 && abs(angle) > max_angle_R
                max_angle_R = abs(angle)
                furthest_vert_R = pt
            end
        end
    end
    return (furthest_vert_L, furthest_vert_R)
end

"""
    closest_crossed_polygon(
        current_point::Point{2, Float64},
        final_point::Point{2, Float64},
        exclusions::DataFrame,
        visited_exclusion_idxs::Vector{Int}
    )::Tuple{
        Tuple{
            AG.IGeometry{AG.wkbPolygon},
            AG.IGeometry{AG.wkbLineString},
            Int32
        },
        Int64
    }

Find polygons that intersect with a line segment.

# Arguments
- `current_point`: The current point marking start of line.
- `final_point`: The final point marking end of line.
- `exclusions`: The dataframe containing the polygon exclusions.
- `visited_exclusion_idxs`: Indices of exclusion zones already crossed.

# Returns
- `closest_polygon`: A tuple containing the:
    - `geometry`: polygon geometry,
    - `exterior_ring`: LineString of the polygon's exterior, and
    - `n_pts`: number of vertices in the polygon,\n
    of the first/closest polygon intersecting with the line segment.
- `polygon_idx`: The index of the (first) polygon crossed.
"""
function closest_crossed_polygon(
    current_point::Point{2, Float64},
    final_point::Point{2, Float64},
    exclusions::DataFrame,
    visited_exclusion_idxs::Vector{Int64}
)::Tuple{
    Tuple{
        AG.IGeometry{AG.wkbPolygon},
        AG.IGeometry{AG.wkbLineString},
        Int32
    },
    Int64
}
    closest_polygon = (AG.creategeom(AG.wkbPolygon), AG.creategeom(AG.wkbLineString), 0)
    min_dist = Inf
    polygon_idx = 0

    line::AG.IGeometry{AG.wkbLineString} = AG.createlinestring(
        [current_point[1], final_point[1]],
        [current_point[2], final_point[2]]
    )

    # Define bounding box for line to exclude polygons that do not intersect
    line_min_x::Float64 = min(current_point[1], final_point[1])
    line_max_x::Float64 = max(current_point[1], final_point[1])
    line_min_y::Float64 = min(current_point[2], final_point[2])
    line_max_y::Float64 = max(current_point[2], final_point[2])

    for (i, geom::AG.IGeometry{AG.wkbPolygon}) ∈ enumerate(exclusions.geometry)
        if i in visited_exclusion_idxs
            continue
        end

        exterior_ring::AG.IGeometry{AG.wkbLineString} = AG.getgeom(geom, 0)
        n_pts::Int32 = AG.ngeom(exterior_ring)

        # Check if any polygon vertices are inside the line's bounding box
        vertex_in_line_bbox = any(
            line_min_x <= AG.getpoint(exterior_ring, j)[1] <= line_max_x &&
            line_min_y <= AG.getpoint(exterior_ring, j)[2] <= line_max_y
            for j in 0:n_pts - 1
        )

        # Check if line crosses bounding box
        line_in_polygon_bbox = false
        if !vertex_in_line_bbox
            poly_xs, poly_ys, _ = [AG.getpoint(exterior_ring, j) for j in 0:n_pts - 1]

            line_in_polygon_bbox = (
                line_min_x <= maximum(view(poly_xs,1:n_pts)) &&
                line_max_x >= minimum(view(poly_xs,1:n_pts)) &&
                line_min_y <= maximum(view(poly_ys,1:n_pts)) &&
                line_max_y >= minimum(view(poly_ys,1:n_pts))
            )

        end

        # Skip polygons with no vertices in or crossing bounding box
        if !vertex_in_line_bbox && !line_in_polygon_bbox
            continue
        end

        if AG.crosses(line, geom) ||
            (
                AG.touches(AG.createpoint(current_point[1], current_point[2]), geom) &&
                AG.touches(AG.createpoint(final_point[1], final_point[2]), geom)
            )

            intersections::AG.IGeometry = AG.intersection(line, geom)

            pts::Vector{AG.IGeometry} = get_pts(intersections)

            # Find distance to polygon
            dist::Float64 = minimum(
                GO.distance.(
                    [current_point]::Vector{Point{2, Float64}},
                    pts
                )
            )

            # If closer than current closest polygon, update closest polygon
            if dist < min_dist
                min_dist = dist
                closest_polygon = (geom, exterior_ring, n_pts)
                polygon_idx = i
            end
        end
    end
    return closest_polygon, polygon_idx
end

function get_pts(intersections::AG.IGeometry{AG.wkbPoint})
    return [intersections]
end
function get_pts(intersections::AG.IGeometry)
    n = AG.ngeom(intersections)
    return AG.getgeom.(Ref(intersections), 0:n-1)
end


"""
    is_visible(
        current_point::Point{2, Float64},
        final_point::Point{2, Float64},
        exclusion_poly::AG.IGeometry{AG.wkbPolygon}
    )::Bool
    is_visible(
        current_point::Point{2, Float64},
        final_point::Point{2, Float64},
        exclusions::DataFrame,
        current_exclusions_idx::Vector{Int} = [0]
    )::Bool

Check if a point is visible from another point, given a vector of, or single exclusion
    polygon/s. i.e. no exclusion zone intersects the straight line between them.

# Arguments
- `current_point`: The current point marking start of line.
- `final_point`: The final point marking end of line.
- `exclusion_poly`: The exclusion polygon.
- `exclusions`: The dataframe containing the polygon exclusions.
- `current_exclusions_idx`: The indices of the exclusion zones that have already been crossed.

# Returns
`true` if the line between the two points is visible, `false` otherwise.
"""
function is_visible(
    current_point::Point{2, Float64},
    final_point::Point{2, Float64},
    exclusion_poly::AG.IGeometry{AG.wkbPolygon}
)::Bool
    line_to_point::AG.IGeometry{AG.wkbLineString} = AG.createlinestring([
        (current_point[1], current_point[2]),
        (final_point[1], final_point[2])
    ])

    # TODO: CHECK if commented statement below is necessary/correct/overkill
    return !AG.intersects(exclusion_poly, line_to_point) ||
        AG.touches(exclusion_poly, line_to_point)
end
function is_visible(
    current_point::Point{2, Float64},
    final_point::Point{2, Float64},
    exclusions::DataFrame,
    current_exclusions_idx::Vector{Int} = [0]
)::Bool
    #! Broadcast/vectorise function
    current_exclusions = Set{Int}()
    sizehint!(current_exclusions, length(current_exclusions_idx))
    for idx in current_exclusions_idx
        push!(current_exclusions, idx)
    end

    geometries::Vector{AG.IGeometry{AG.wkbPolygon}} = exclusions.geometry
    line_to_point::AG.IGeometry{AG.wkbLineString} = AG.createlinestring([
        (current_point[1], current_point[2]),
        (final_point[1], final_point[2])
    ])

    isempty_geo = AG.isempty

    for (i, geometry) ∈ enumerate(geometries)
        # Skip exclusion zones
        if isempty_geo(geometry) || (i ∈ current_exclusions)
            continue
        end

        if !is_visible(current_point, final_point, geometry_i)
            return false
        end
    end
    return true
end
