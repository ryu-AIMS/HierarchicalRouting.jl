
"""
    find_next_points(
    current_point::Point{2, Float64},
    final_point::Point{2, Float64},
    exclusions::DataFrame,
    current_exclusions_idx::Vector{Int}
    )

Find the widest vertices on each (left/right) side of the line to the final point.

# Arguments
- `current_point::Point{2, Float64}`: The current point marking start of line.
- `final_point::Point{2, Float64}`: The final point marking end of line.
- `exclusions::DataFrame`: The dataframe containing the polygon exclusions.
- `current_exclusions_idx::Vector{Int}`: The indices of the exclusion polygons to disregard for crossings.

# Returns
- `[furthest_vert_L, furthest_vert_R]`: The furthest left and right polygon vertices of the line to the final point.
- `polygon_idx::Int`: The index of the polygon crossed.
"""
function find_next_points(
    current_point::Point{2, Float64},
    final_point::Point{2, Float64},
    exclusions::DataFrame,
    current_exclusions_idx::Vector{Int}
)
    max_dist_L, max_dist_R = -Inf32, -Inf32
    furthest_vert_L, furthest_vert_R = nothing, nothing # Set{Point{2, Float64}}()

    # TODO: Allow for multiple equidistant furthest points

    # Next polygon crossed by line to end point
    (polygon, exterior_ring, n_pts), polygon_idx = closest_crossed_polygon(current_point, final_point, exclusions, current_exclusions_idx)
    if polygon === nothing
        return [final_point], 0
    end

    """
        perpendicular_distance_line_to_point(line::Line{2, Float64}, point::Point{2, Float64})

    Calculate the perpendicular distance of a point to a line.

    # Arguments
    - `line::Line{2, Float64}`: The line to calculate the distance to.
    - `point::Point{2, Float64}`: The point to calculate the distance from.

    # Returns
    - `dist::Float64`: The perpendicular distance of the point to the line.
    """
    function perpendicular_distance_line_to_point(line::Line{2, Float64}, point::Point{2, Float64})
        p1, p2 = line.points[1], line.points[2]

        # Convert points to vectors
        v1 = Vec(p1)
        v2 = Vec(p2)
        p = Vec(point)

        # Vector projection and the distance
        line_vec = v2 - v1
        point_vec = p - v1
        proj_len = GeometryBasics.dot(point_vec, line_vec) / GeometryBasics.norm(line_vec)
        proj_point = v1 + proj_len * GeometryBasics.normalize(line_vec)

        return GeometryBasics.norm(p - proj_point)
    end

    # For each polygon vertex, find the furthest points on the left and right side of line
    for i in 0:n_pts - 1
        x, y, _ = AG.getpoint(exterior_ring, i)
        pt = Point{2, Float64}(x, y)

        dist = perpendicular_distance_line_to_point(Line(current_point, final_point), pt)

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

    return [furthest_vert_L, furthest_vert_R], polygon_idx
end

"""
    closest_crossed_polygon(
    current_point::Point{2, Float64},
    final_point::Point{2, Float64},
    exclusions::DataFrame,
    current_exclusions_idx::Vector{Int}
    )

Find polygons that intersect with a line segment.

# Arguments
- `current_point::Point{2, Float64}`: The current point marking start of line.
- `final_point::Point{2, Float64}`: The final point marking end of line.
- `exclusions::DataFrame`: The dataframe containing the polygon exclusions.
- `current_exclusions_idx::Vector{Int}`: The indices of the exclusion polygons to disregard.

# Returns
- `closest_polygon`: The polygon, LineString and number of vertices of the first/closest polygon intersecting with the line segment.
- `polygon_idx`: The index of the(first) polygon crossed.
"""
function closest_crossed_polygon(
    current_point::Point{2, Float64},
    final_point::Point{2, Float64},
    exclusions::DataFrame,
    current_exclusions_idx::Vector{Int}
)
    closest_polygon = (nothing, nothing, 0)
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
