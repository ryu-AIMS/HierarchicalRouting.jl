
"""
    find_widest_points(
        current_point::Point{2, Float64},
        final_point::Point{2, Float64},
        exclusions::DataFrame,
        current_exclusion_idx::Int
    )::Tuple{Vector{Point{2, Float64}}, Vector{Int}}

Find the intermediate vertices between two points by using the widest visible vertices on
each (left/right) side of the direct-line path by closest crossed polygon.
If these intermediate vertices are not visible, the function is called recursively to find
the widest visible vertices.

# Arguments
- `current_point::Point{2, Float64}`: The current point marking start of line.
- `final_point::Point{2, Float64}`: The final point marking end of line.
- `exclusions::DataFrame`: The dataframe containing the polygon exclusions.
- `current_exclusion_idx::Int`: The index of the current exclusion zone - to ignore.

# Returns
- A vector of points representing the widest visible points on each side of the line.
- A vector of polygon indices crossed by the line.
"""
function find_widest_points(
    current_point::Point{2, Float64},
    final_point::Point{2, Float64},
    exclusions::DataFrame,
    current_exclusion_idx::Int
)::Tuple{Vector{Point{2, Float64}}, Vector{Int}}

    if is_visible(current_point, final_point, exclusions)
        return [final_point], [0]
    end
    # Find the first/next polygon that the line (current_point, final_point) crosses.
    (polygon, exterior_ring, n_pts), polygon_idx = HierarchicalRouting.
        closest_crossed_polygon(current_point, final_point, exclusions, 0)

    if isnothing(polygon)
        return [final_point], [0]
    end

    max_angle_L, max_angle_R = 0.0, 0.0
    furthest_vert_L, furthest_vert_R = nothing, nothing

    # Base vector from current_point to final_point.
    base_vec = [
        (final_point[1] - current_point[1]),
        (final_point[2] - current_point[2])
    ]

    for i in 0:n_pts - 1
        x, y, _ = AG.getpoint(exterior_ring, i)
        pt = Point{2, Float64}(x, y)

        # Vector from current_point to this vertex.
        vec = [
            (pt[1] - current_point[1]),
            (pt[2] - current_point[2])
        ]

        # Signed angle between base vector and vector to vertex
        angle = HierarchicalRouting.vector_angle(base_vec, vec)

        if HierarchicalRouting.is_visible(current_point, pt, exclusions) &&
            pt != current_point
            # TODO Try to avoid going backwards - i.e. angle > abs(π/2)
            # ! only the case when point is in the convex hull of the exclusion zone
            if angle > 0 && angle > max_angle_L
                max_angle_L = angle
                furthest_vert_L = pt
            elseif angle < 0 && abs(angle) > max_angle_R
                max_angle_R = abs(angle)
                furthest_vert_R = pt
            end
        end
    end

    candidates = Point{2,Float64}[]
    poly_indices = Int[]

    for vertex ∈ [furthest_vert_L, furthest_vert_R]
        if !isnothing(vertex)
            if is_visible(current_point, vertex, exclusions, [current_exclusion_idx])
                push!(candidates, vertex)
                push!(poly_indices, polygon_idx)
            else
                new_pts, new_poly_idx = find_widest_points(
                    current_point,
                    vertex,
                    exclusions,
                    current_exclusion_idx
                )
                append!(candidates, new_pts)
                append!(poly_indices, new_poly_idx)
            end
        end
    end

    return candidates, poly_indices
end

# Compute the signed angle (radian) between base_vec and vec
function vector_angle(base_vec::Vector{Float64}, vec::Vector{Float64})

    # Cross product tells us the direction of the angle
    cross = base_vec[1]*vec[2] - base_vec[2]*vec[1]
    # Dot product tells us the magnitude of the angle
    dot_val = base_vec[1]*vec[1] + base_vec[2]*vec[2]

    return atan(cross, dot_val)
end

"""
    closest_crossed_polygon(
    current_point::Point{2, Float64},
    final_point::Point{2, Float64},
    exclusions::DataFrame,
    current_exclusionn_idx::Int
    )

Find polygons that intersect with a line segment.

# Arguments
- `current_point::Point{2, Float64}`: The current point marking start of line.
- `final_point::Point{2, Float64}`: The final point marking end of line.
- `exclusions::DataFrame`: The dataframe containing the polygon exclusions.
- `current_exclusion_idx::Int`: The index of the current exclusion zone - to ignore.

# Returns
- `closest_polygon`: The polygon, LineString and number of vertices of the first/closest polygon intersecting with the line segment.
- `polygon_idx`: The index of the(first) polygon crossed.
"""
function closest_crossed_polygon(
    current_point::Point{2, Float64},
    final_point::Point{2, Float64},
    exclusions::DataFrame,
    current_exclusion_idx::Int
)
    closest_polygon = (nothing, nothing, 0)
    min_dist = Inf
    polygon_idx = nothing

    line = AG.createlinestring([current_point[1], final_point[1]], [current_point[2], final_point[2]])

    # Define bounding box for line to exclude polygons that do not intersect
    line_min_x, line_max_x = min(current_point[1], final_point[1]), max(current_point[1], final_point[1])
    line_min_y, line_max_y = min(current_point[2], final_point[2]), max(current_point[2], final_point[2])

    for (i, row) in enumerate(eachrow(exclusions))
        if i == current_exclusion_idx
            continue
        end

        geom = row.geometry
        exterior_ring = AG.getgeom(geom, 0)
        n_pts = AG.ngeom(exterior_ring)

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
                line_min_x <= maximum(poly_xs) &&
                line_max_x >= minimum(poly_xs) &&
                line_min_y <= maximum(poly_ys) &&
                line_max_y >= minimum(poly_ys)
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

            # Find distance to polygon
            dist = GO.distance(current_point, geom)

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

"""
    is_visible(
    current_point::Point{2, Float64},
    final_point::Point{2, Float64},
    exclusion_poly
)::Bool

Check if a point is visible from another point, given a single exclusion polygon.
    i.e. no exclusion zone intersects the straight line between them.

# Arguments
- `current_point::Point{2, Float64}`: The current point marking start of line.
- `final_point::Point{2, Float64}`: The final point marking end of line.
- `exclusion_poly`: The exclusion polygon.

# Returns
- `Bool`: True if the line between the two points is visible, false otherwise.
"""
function is_visible(
    current_point::Point{2, Float64},
    final_point::Point{2, Float64},
    exclusion_poly
)::Bool

    line_to_point = AG.createlinestring([
        (current_point[1], current_point[2]),
        (final_point[1], final_point[2])
    ])

    # TODO: CHECK if commented statement below is necessary/correct/overkill
    return !AG.intersects(exclusion_poly, line_to_point) || AG.touches(exclusion_poly, line_to_point)
end

"""
    is_visible(
    current_point::Point{2, Float64},
    final_point::Point{2, Float64},
    exclusions::DataFrame,
    current_exclusions_idx::Vector{Int} = [0]
    )::Bool

Check if a point is visible from another point, considering a vector of polygons.
    i.e. no exclusion zone intersects the straight line between them.

# Arguments
- `current_point::Point{2, Float64}`: The current point marking start of line.
- `final_point::Point{2, Float64}`: The final point marking end of line.
- `exclusions::DataFrame`: The dataframe containing the polygon exclusions.
- `current_exclusions_idx::Vector{Int}`: The indices of the exclusion zones that have already been crossed.

# Returns
- `Bool`: True if the line between the two points is visible, false otherwise.
"""
function is_visible(
    current_point::Point{2, Float64},
    final_point::Point{2, Float64},
    exclusions::DataFrame,
    current_exclusions_idx::Vector{Int} = [0]
)::Bool

    for i in 1:size(exclusions, 1)

        # Skip exclusion zones
        if i in current_exclusions_idx || AG.isempty(exclusions.geometry[i])
            continue
        end

        if !is_visible(current_point, final_point, exclusions.geometry[i])
            return false
        end
    end
    return true
end
