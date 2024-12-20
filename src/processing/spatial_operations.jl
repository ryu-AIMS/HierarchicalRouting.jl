
"""
    simplify_exclusions!(
    exclusions::DataFrame;
    min_area::Real=50,
    simplify_tol::Real=2,
    convex_flag::Bool=true
    )::DataFrame

Simplify exclusions by applying a convex hull, removing small polygons, and simplifying.

# Arguments
- `exclusions::DataFrame`: The DataFrame containing exclusion zones.
- `min_area::Real`: The minimum area for exclusion polygons. Default = 50
- `simplify_tol::Real`: The tolerance value for simplifying the exclusion polygons. Default = 2
- `convex_flag::Bool`: Apply convex hull to exclusion polygons. Default = true
"""
function simplify_exclusions!(
    exclusions::DataFrame;
    min_area::Real=50,
    simplify_tol::Real=2,
    convex_flag::Bool=true
    )
    if convex_flag
        exclusions.geometry .= AG.convexhull.(exclusions.geometry)
    end
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
function buffer_exclusions!(exclusions::DataFrame; buffer_dist::Real=1.0)
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

    for i in 1:n
        geom1 = geometries[i]

        for j in i+1:n
            geom2 = geometries[j]

            if AG.overlaps(geom1, geom2)
                @debug "Partial overlap: $i and $j"
                union_geom = AG.union(geom1, geom2)
                exclusions.geometry[i] = union_geom
                exclusions.geometry[j] = union_geom

                for k in 1:n
                    if AG.overlaps(union_geom, geometries[k])
                        exclusions.geometry[k] = union_geom
                    end
                end
            end

            if AG.contains(geom1, geom2)
                @debug "Full overlap: $i contains $j"
                exclusions.geometry[j] = geom1
            elseif AG.contains(geom2, geom1)
                @debug "Full overlap: $j contains $i"
                exclusions.geometry[i] = geom2
            end
        end
    end

    # Remove duplicate unionized geometries
    unique_geometries = unique(exclusions.geometry[.!AG.isempty.(exclusions.geometry)])
    exclusions = DataFrame(geometry = unique_geometries)

    return exclusions
end
