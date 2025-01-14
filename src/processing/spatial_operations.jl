
"""
    simplify_exclusions!(
    exclusions::DataFrame;
    min_area::Real=50,
    simplify_tol::Real=2
    )::DataFrame

Simplify exclusions by removing small polygons, and simplifying geometry.

# Arguments
- `exclusions::DataFrame`: The DataFrame containing exclusion zones.
- `min_area::Real`: The minimum area for exclusion polygons. Default = 50
- `simplify_tol::Real`: The tolerance value for simplifying the exclusion polygons, i.e., larger tol = more aggressive simplification. Default = 2
"""
function simplify_exclusions!(
    exclusions::DataFrame;
    min_area::Real=50,
    simplify_tol::Real=2
    )
    exclusions = exclusions[AG.geomarea.(exclusions.geometry) .>= min_area, :]
    exclusions.geometry .= AG.simplify.(exclusions.geometry, simplify_tol)
    return exclusions
end

"""
    buffer_exclusions!(exclusions::DataFrame, buffer_dist=1.0)::DataFrame

Buffer exclusion zones by a specified distance.

# Arguments
- `exclusions::DataFrame`: The DataFrame containing exclusion zones.
- `buffer_dist::Real`: The buffer distance. Default = 1.0
"""
function buffer_exclusions!(exclusions::DataFrame; buffer_dist::Real=2.0)
    exclusions.geometry .= AG.buffer.(exclusions.geometry, buffer_dist)
    return exclusions
end

"""
    unionize_overlaps(exclusions::DataFrame)::DataFrame

Unionize overlapping exclusion zones.

# Arguments
- `exclusions::DataFrame`: The DataFrame containing exclusion zones.
"""
function unionize_overlaps!(exclusions::DataFrame)
    geometries = exclusions.geometry
    n = length(geometries)

    """
        linestring_to_polygon(linestring::AG.IGeometry{AG.wkbLineString})::AG.IGeometry{AG.wkbPolygon}

    Convert a LineString to a Polygon.

    # Arguments
    - `linestring::AG.IGeometry{AG.wkbLineString}`: The LineString to convert.

    # Returns
    - `polygon::AG.IGeometry{AG.wkbPolygon}`: The converted Polygon.
    """
    function linestring_to_polygon(linestring::AG.IGeometry{AG.wkbLineString})
        num_points = AG.ngeom(linestring)
        points = [(AG.getx(linestring, i), AG.gety(linestring, i)) for i in 0:num_points-1]

        # Close the LineString if open
        if points[1] != points[end]
            push!(points, points[1])
        end

        return AG.createpolygon([points])
    end

    for i in 1:n
        geom1 = geometries[i]

        if AG.ngeom(geom1) > 1
            # Unionize multi-geometries
            geom_a = linestring_to_polygon(AG.getgeom(geom1, 0))
            for j in 1:AG.ngeom(geom1) - 1
                geom_b = linestring_to_polygon(AG.getgeom(geom1, j))

                if AG.intersects(geom_a, geom_b) #|| AG.overlaps(geom_a, geom_b) || AG.contains(geom_a, geom_b) || AG.contains(geom_b, geom_a)
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
    exclusions = DataFrame(geometry = unique_geometries)

    return exclusions
end
