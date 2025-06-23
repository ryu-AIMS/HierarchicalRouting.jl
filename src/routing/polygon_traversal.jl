
"""
    find_widest_points(
        current_point::Point{2,Float64},
        final_point::Point{2,Float64},
        exclusions::Vector{IGeometry{wkbPolygon}},
        visited_exclusion_idxs::Vector{Int}=[0]
    )::Tuple{Vector{Point{2,Float64}},Vector{Int}}

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
    current_point::Point{2,Float64},
    final_point::Point{2,Float64},
    exclusions::POLY_VEC,
    visited_exclusion_idxs::Vector{Int}=[0]
)::Tuple{Vector{Point{2,Float64}},Vector{Int}}

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
    widest_vertices::NTuple{2,Union{Point{2,Float64},Nothing}} = search_widest_points(
        current_point, final_point, exterior_ring, n_pts
    )

    candidates = Point{2,Float64}[]
    poly_indices = Int[]
    sizehint!(candidates, length(widest_vertices))
    sizehint!(poly_indices, length(widest_vertices))

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
        base_vector::Vector{Float64}, candidate_tuple::NTuple{2,Float64}
    )::Float64
    vector_angle(
        base_vector::Vector{Float64}, candidate_vector::Vector{Float64}
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
    base_vector::Vector{Float64}, candidate_tuple::NTuple{2,Float64}
)::Float64
    return vector_angle(base_vector, collect(candidate_tuple))
end
function vector_angle(
    base_vector::Vector{Float64}, candidate_vector::Vector{Float64}
)::Float64
    # Cross product tells us the direction of the angle
    cross = base_vector[1] * candidate_vector[2] - base_vector[2] * candidate_vector[1]
    # Dot product tells us the magnitude of the angle
    dot_val = base_vector[1] * candidate_vector[1] + base_vector[2] * candidate_vector[2]

    return atan(cross, dot_val)
end

"""
    search_widest_points(
        current_point::Point{2,Float64},
        final_point::Point{2,Float64},
        exterior_ring::IGeometry{wkbLineString},
        n_pts::Int32,
    )::NTuple{2,Union{Point{2,Float64},Nothing}}

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
    current_point::Point{2,Float64},
    final_point::Point{2,Float64},
    exterior_ring::IGeometry{wkbLineString},
    n_pts::Int32,
)::NTuple{2,Union{Point{2,Float64},Nothing}}
    max_angle_L, max_angle_R = 0.0, 0.0
    furthest_vert_L, furthest_vert_R = nothing, nothing
    # Base vector from current_point to final_point to calc angle from
    base_vector = [
        final_point[1] - current_point[1],
        final_point[2] - current_point[2]
    ]

    x::Float64, y::Float64 = 0.0, 0.0
    pt = Point{2,Float64}(x, y)
    candidate_tuple::Tuple{Float64,Float64} = (0.0, 0.0)
    angle::Float64 = 0.0
    for i in 0:n_pts-1
        x, y, _ = AG.getpoint(exterior_ring, i)
        pt = Point{2,Float64}(x, y)

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
        current_point::Point{2,Float64},
        final_point::Point{2,Float64},
        exclusions::Vector{IGeometry{wkbPolygon}},
        visited_exclusion_idxs::Vector{Int}
    )::Tuple{
        Tuple{
            IGeometry{wkbPolygon},
            IGeometry{wkbLineString},
            Int32
        },
        Int64
    }

Find the first polygon that is intersected by a line between two given points.

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
    current_point::Point{2,Float64},
    final_point::Point{2,Float64},
    exclusions::POLY_VEC,
    visited_exclusion_idxs::Vector{Int64}
)::Tuple{
    Tuple{
        IGeometry{wkbPolygon},
        IGeometry{wkbLineString},
        Int32
    },
    Int64
}
    # Create placeholder variables
    exterior_ring::IGeometry{wkbLineString} = AG.creategeom(wkbLineString)
    n_pts::Int32 = 0
    pt_ref = Ref{NTuple{3,Float64}}()
    vertex_in_line_bbox::Bool = false
    line_in_polygon_bbox::Bool = false

    closest_polygon = (AG.creategeom(wkbPolygon), AG.creategeom(wkbLineString), 0)
    min_dist = Inf
    polygon_idx = 0

    # Set vector lengths for polygon vertices with length of 100
    poly_xs = Vector{Float64}(undef, 100)
    poly_ys = Vector{Float64}(undef, 100)

    # Define bounding box for line to exclude polygons that do not intersect
    line::IGeometry{wkbLineString} = AG.createlinestring(
        Float64[current_point[1], final_point[1]],
        Float64[current_point[2], final_point[2]]
    )

    line_min_x::Float64 = min(current_point[1], final_point[1])
    line_max_x::Float64 = max(current_point[1], final_point[1])
    line_min_y::Float64 = min(current_point[2], final_point[2])
    line_max_y::Float64 = max(current_point[2], final_point[2])

    for (i, geom::IGeometry{wkbPolygon}) ∈ enumerate(exclusions)
        if i in visited_exclusion_idxs
            continue
        end

        exterior_ring = AG.getgeom(geom, 0)
        n_pts = AG.ngeom(exterior_ring)

        # Resize bitvector and arrays if length of polygon vertices exceeds current size
        if n_pts > length(poly_xs)
            size_diff = n_pts - length(poly_xs)
            append!(poly_xs, Vector{Float64}(undef, size_diff))
            append!(poly_ys, Vector{Float64}(undef, size_diff))
        end

        # Check if any polygon vertices are inside the line's bounding box
        vertex_in_line_bbox = false
        @inbounds for i in 1:n_pts
            pt_ref = AG.getpoint(exterior_ring, i - 1)
            x = poly_xs[i] = pt_ref[1]
            y = poly_ys[i] = pt_ref[2]

            if (line_min_x ≤ x ≤ line_max_x) && (line_min_y ≤ y ≤ line_max_y)
                vertex_in_line_bbox = true
                break
            end
        end

        # Check if line crosses bounding box
        line_in_polygon_bbox = !vertex_in_line_bbox ? (
            line_min_x <= maximum(view(poly_xs, 1:n_pts)) &&
            line_max_x >= minimum(view(poly_xs, 1:n_pts)) &&
            line_min_y <= maximum(view(poly_ys, 1:n_pts)) &&
            line_max_y >= minimum(view(poly_ys, 1:n_pts))
        ) : false

        # Skip polygons with no vertices in or crossing bounding box
        if !vertex_in_line_bbox && !line_in_polygon_bbox
            continue
        end

        # Check if line crosses polygon or touches both current and final points
        if AG.crosses(line, geom) ||
           (GO.touches(current_point, geom) && GO.touches(final_point, geom))
            intersections::IGeometry = AG.intersection(line, geom)
            pts::Vector{IGeometry{wkbPoint}} = get_pts(intersections)
            dist::Float64 = minimum(
                GO.distance.(
                    Ref{Point{2,Float64}}(current_point),
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

"""
    get_pts(intersections::IGeometry{AG.wkbPoint})::Vector{IGeometry{wkbPolygon}}
    get_pts(intersections::IGeometry)::Vector{IGeometry{wkbPolygon}}

Get the points from the intersection geometry.

# Arguments
- `intersections`: The intersection geometry.

# Returns
A vector of points from the intersection geometry.
"""
function get_pts(intersections::IGeometry{AG.wkbPoint})::Vector{IGeometry{wkbPoint}}
    return [intersections]
end
function get_pts(intersections::IGeometry{wkbLineString})::Vector{IGeometry{wkbPoint}}
    n = AG.ngeom(intersections)
    return AG.getgeom.(Ref(intersections), 0:n-1)
end
function get_pts(intersections::IGeometry{wkbMultiLineString})::Vector{IGeometry{wkbPoint}}
    n = AG.ngeom(intersections)
    return AG.getgeom.(Ref(intersections), 0:n-1)
end

"""
    is_visible(
        line_to_point::IGeometry{wkbLineString},
        exclusion_poly::IGeometry{wkbPolygon}
    )::Bool
    is_visible(
        current_point::Point{2,Float64},
        final_point::Point{2,Float64},
        exclusion_poly::IGeometry{wkbPolygon}
    )::Bool
    is_visible(
        current_point::Point{2,Float64},
        final_point::Point{2,Float64},
        exclusions::Vector{IGeometry{wkbPolygon}},
        current_exclusions_idx::Vector{Int}=[0]
    )::Bool

Check if a point is visible from another point, given a vector of, or single exclusion
    polygon/s. i.e. no exclusion zone intersects the straight line between them.

# Arguments
- `line_to_point`: The line from the current point to the final point.
- `current_point`: The current point marking start of line.
- `final_point`: The final point marking end of line.
- `exclusion_poly`: The exclusion polygon.
- `exclusions`: The dataframe containing the polygon exclusions.
- `current_exclusions_idx`: The indices of the exclusion zones that have already been crossed.

# Returns
`true` if the line between the two points is visible, `false` otherwise.
"""
function is_visible(
    line_to_point::IGeometry{wkbLineString},
    exclusion_poly::IGeometry{wkbPolygon}
)::Bool
    return !AG.intersects(exclusion_poly, line_to_point) ||
           AG.touches(exclusion_poly, line_to_point)
end
function is_visible(
    current_point::Point{2,Float64},
    final_point::Point{2,Float64},
    exclusion_poly::IGeometry{wkbPolygon}
)::Bool
    line_to_point::IGeometry{wkbLineString} = AG.createlinestring([
        (current_point[1], current_point[2]),
        (final_point[1], final_point[2])
    ])
    return is_visible(line_to_point, exclusion_poly)
end
function is_visible(
    current_point::Point{2,Float64},
    final_point::Point{2,Float64},
    exclusions::DataFrame,
    current_exclusions_idx::Vector{Int64}=[0]
)::Bool
    return is_visible(current_point, final_point, exclusions.geometry, current_exclusions_idx)
end
function is_visible(
    current_point::Point{2,Float64},
    final_point::Point{2,Float64},
    geometries::POLY_VEC,
    current_exclusions_idx::Vector{Int}=[0]
)
    current_exclusions = Set{Int64}(current_exclusions_idx)

    line_to_point::IGeometry{wkbLineString} = AG.createlinestring([
        (current_point[1], current_point[2]),
        (final_point[1], final_point[2])
    ])

    for (i, geometry) ∈ enumerate(geometries)
        # Skip exclusion zones
        if AG.isempty(geometry) || (i ∈ current_exclusions)
            continue
        end

        if !is_visible(line_to_point, geometry)
            return false
        end
    end

    return true
end
