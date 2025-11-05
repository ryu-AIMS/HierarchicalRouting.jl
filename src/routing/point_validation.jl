
"""
    point_in_exclusion(point::Point{2,Float64}, exclusion::IGeometry{wkbPolygon})::Bool
    point_in_exclusion(point::Point{2,Float64}, exclusions::Vector{IGeometry{wkbPolygon}})::Bool

Check if a point is within an exclusion zone.

# Arguments
- `point`: Point to check.
- `exclusion::IGeometry`: One exclusion zone polygon/geometry from a DataFrame.
- `exclusions::Vector{IGeometry{wkbPolygon}}`: A vector of exclusion zone polygons/geometry.

# Returns
- `true` if point is within an exclusion zone, `false` otherwise.
"""
function point_in_exclusion(point::Point{2,Float64}, exclusion::IGeometry{wkbPolygon})::Bool
    return any(AG.contains(exclusion, AG.createpoint(point.data)))
end
function point_in_exclusion(point::Point{2,Float64}, exclusions::POLY_VEC)::Bool
    return any(AG.contains.(exclusions, Ref(AG.createpoint(point.data))))
end

"""
    containing_exclusion(point::Point{2,Float64}, exclusions::Vector{IGeometry{wkbPolygon}})::Int

Return the index of the exclusion zone that contains the point.

# Arguments
- `point`: Point to check.
- `exclusions`: A DataFrame containing exclusion zone polygons.

# Returns
- Index of the first exclusion zone that contains the point, or 0 if not found.
"""
function containing_exclusion(point::Point{2,Float64}, exclusions::POLY_VEC)::Int
    point_ag = AG.createpoint(point.data)
    exclusion_idx = findfirst(AG.contains.(exclusions, Ref(point_ag)))
    return isnothing(exclusion_idx) ? 0 : exclusion_idx
end

"""
    point_in_convexhull(
        point::Point{2,Float64},
        exclusions::Vector{IGeometry{wkbPolygon}}
    )::Int

Check if a point is within a convex hull of exclusion zones.

# Arguments
- `point`: Point to check.
- `exclusions`: Exclusion zones.

# Returns
- Index of exclusion zone if point is within a convex hull, 0 otherwise.
"""
function point_in_convexhull(point::Point{2,Float64}, exclusions::POLY_VEC)
    point_ag::AG.IGeometry{wkbPoint} = AG.createpoint(point.data)
    convex_exclusions_ag::POLY_VEC = [AG.convexhull(exclusion) for exclusion in exclusions]

    point_in_exclusion_zone = AG.contains.(convex_exclusions_ag, [point_ag])

    # If point is within an exclusion zone, add all polygon vertices to graph
    # ? Consider cases where point is within the convex hull of multiple polygons
    # findall(final_point_in_exclusion_zone)
    exclusion_index = findfirst(point_in_exclusion_zone)

    return isnothing(exclusion_index) ? 0 : exclusion_index
end
