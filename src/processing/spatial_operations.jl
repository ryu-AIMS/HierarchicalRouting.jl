
"""
    filter_and_simplify_exclusions!(
        exclusions::DataFrame;
        min_area::Float64=3E-5,
        simplify_tol::Float64=1E-3
    )::DataFrame

Remove small polygons and simplify geometry based on tolerance values.

# Arguments
- `exclusions`: The DataFrame containing exclusion zones.
- `min_area`: The minimum area for exclusion polygons.
- `simplify_tol`: The tolerance value for simplifying the exclusion polygons, \n
    i.e., larger tol = more aggressive simplification.
"""
function filter_and_simplify_exclusions!(
    exclusions::DataFrame;
    min_area::Float64=1E-5,
    simplify_tol::Float64=1E-3
)::DataFrame
    exclusions.geometry .= AG.simplify.(exclusions.geometry, simplify_tol)
    filter!(row -> AG.geomarea(row.geometry) >= min_area, exclusions)
    return exclusions
end

"""
    buffer_exclusions!(exclusions::DataFrame; buffer_dist::Float64=0.0)::DataFrame

Buffer exclusion zones by a specified distance.

# Arguments
- `exclusions`: The DataFrame containing exclusion zones.
- `buffer_dist`: The buffer distance.
"""
function buffer_exclusions!(exclusions::DataFrame; buffer_dist::Float64=0.0)::DataFrame
    exclusions.geometry .= AG.buffer.(exclusions.geometry, buffer_dist)
    return exclusions
end

"""
    linestring_to_polygon(
        linestring::IGeometry{wkbLineString}
    )::IGeometry{wkbPolygon}

Convert a LineString to a Polygon.

# Arguments
- `linestring`: The LineString to convert.

# Returns
The converted Polygon.
"""
function linestring_to_polygon(
    linestring::IGeometry{wkbLineString}
)::IGeometry{wkbPolygon}
    num_points::Int = AG.ngeom(linestring)
    points::Vector{NTuple{2,Float64}} = collect(
        zip(
            AG.getx.(Ref(linestring), 0:num_points-1),
            AG.gety.(Ref(linestring), 0:num_points-1)
        )
    )

    # Close the LineString if open
    points = points[1] != points[end] ? vcat(points, points[1]) : points

    return AG.createpolygon([points])
end

"""
    unionize_overlaps(exclusions::DataFrame)::DataFrame

Unionize overlapping exclusion zones.

# Arguments
- `exclusions`: DataFrame containing exclusion zones.
"""
function unionize_overlaps!(exclusions::DataFrame)::DataFrame
    geometries = exclusions.geometry
    n = length(geometries)

    for i in 1:n
        geom1 = geometries[i]

        if AG.ngeom(geom1) > 1
            # Unionize multi-geometries
            geom_a = linestring_to_polygon(AG.getgeom(geom1, 0))
            for j in 1:AG.ngeom(geom1)-1
                geom_b = linestring_to_polygon(AG.getgeom(geom1, j))
                if AG.intersects(geom_a, geom_b)
                    geom_a = AG.union(geom_a, geom_b)
                end
            end
            geometries[i] = geom_a
            geom1 = geom_a
        end

        for j in i+1:n
            geom2 = geometries[j]

            if AG.overlaps(geom1, geom2)
                @debug "Partial overlap: $i and $j"
                union_geom = AG.union(geom1, geom2)
                geometries[i] = union_geom
                geometries[j] = union_geom

                for k in 1:n
                    if AG.overlaps(union_geom, geometries[k])
                        geometries[k] = union_geom
                    end
                end
            end

            if AG.contains(geom1, geom2)
                @debug "Full overlap: $i contains $j"
                geometries[j] = geom1
            elseif AG.contains(geom2, geom1)
                @debug "Full overlap: $j contains $i"
                geometries[i] = geom2
            end
        end
    end

    # Remove duplicate unionized geometries
    unique_geometries = unique(geometries[.!AG.isempty.(geometries)])

    empty!(exclusions)
    append!(exclusions, [(geometry=geom,) for geom in unique_geometries])

    return exclusions
end
