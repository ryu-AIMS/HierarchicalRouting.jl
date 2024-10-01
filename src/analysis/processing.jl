
"""
    to_multipolygon(raster::Raster{T, 2}) where {T<:Union{Bool,Int16}}

Convert raster to multipolygons.
Invert vertical axis.
"""
# TODO: Check polygonize coordinates output
function to_multipolygon(raster::Raster{T, 2})::GI.Wrappers.MultiPolygon where {T<:Union{Bool,Int16}}
    return GO.polygonize(.==(0), raster[:, end:-1:1])
end

"""
    to_dataframe(mp::GI.Wrappers.MultiPolygon)::DataFrame

Create a DataFrame from multipolygons
"""
function to_dataframe(mp::GI.Wrappers.MultiPolygon)::DataFrame
    return DataFrame(geometry=mp.geom)
end

"""
    simplify_exclusions!(exclusions::DataFrame; min_area::Real=200, simplify_tol::Real=10)::DataFrame

Simplify exclusions by applying a convex hull, removing small polygons, and simplifying.

# Arguments
- `exclusions::DataFrame`: The DataFrame containing exclusion zones.
- `min_area::Real`: The minimum area for exclusion polygons. Default = 200
- `simplify_tol::Real`: The tolerance value for simplifying the exclusion polygons. Default = 10
"""
function simplify_exclusions!(exclusions::DataFrame; min_area::Real=200, simplify_tol::Real=10)
    # TODO: Update default args
    exclusions.geometry .= AG.convexhull.(exclusions.geometry)
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
                # @info "Partial overlap: $i and $j"
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
                # @info "Full overlap: $i contains $j"
                exclusions.geometry[j] = geom1
            elseif AG.contains(geom2, geom1)
                # @info "Full overlap: $j contains $i"
                exclusions.geometry[i] = geom2
            end
        end
    end

    # Remove duplicate unionized geometries
    unique_geometries = unique(exclusions.geometry[.!AG.isempty.(exclusions.geometry)])
    exclusions = DataFrame(geometry = unique_geometries)

    return exclusions
end
