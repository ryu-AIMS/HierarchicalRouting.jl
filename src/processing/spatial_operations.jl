
"""
    filter_and_simplify_exclusions(
        old_geoms::POLY_VEC;
        min_area::Float64,
        simplify_tol::Float64
    )::POLY_VEC

Simplify polygon geometries based on tolerance values, whilst ensuring they remain polygons,
then remove small and empty polygons.

# Arguments
- `old_geoms`: A vector of polygon geometries to filter and simplify.
- `min_area`: The minimum area (in square degrees) for exclusion polygons.
- `simplify_tol`: The tolerance value for simplifying the exclusion polygons, \n
    i.e., larger tol = more aggressive simplification.
"""
function filter_and_simplify_exclusions(
    old_geoms::POLY_VEC;
    min_area::Float64,
    simplify_tol::Float64
)::POLY_VEC
    tmp_geoms::Vector{<:IGeometry} = AG.simplify.(old_geoms, simplify_tol)
    explode_multipolygons!(tmp_geoms)

    new_geoms::POLY_VEC = filter(
        geom -> AG.geomarea(geom) >= min_area && !AG.isempty(geom),
        tmp_geoms
    )
    return new_geoms
end

"""
    Explodes MultiPolygons in a vector of geometries to individual Polygons, ensuring the
    resulting vector is a POLY_VEC.
"""
function explode_multipolygons!(geometries::Vector{<:IGeometry})::POLY_VEC
    multi_poly_idxs = findall(==(IGeometry{wkbMultiPolygon}), typeof.(geometries))

    for idx in multi_poly_idxs
        multi_poly::IGeometry{wkbMultiPolygon} = geometries[idx]
        n::Int = AG.ngeom(multi_poly)

        tmp_vec::Vector{IGeometry{wkbPolygon}} = AG.getgeom.(Ref(multi_poly), 0:n-1)
        append!(geometries, tmp_vec)
    end

    # Delete MultiPolygons after appending individual Polygons
    deleteat!(geometries, multi_poly_idxs)

    return POLY_VEC(geometries)
end

"""
    buffer_exclusions!(
        exclusion_geometries::POLY_VEC;
        buffer_dist::Float64=0.0
    )::POLY_VEC
    buffer_exclusions(
        exclusion_geometries::POLY_VEC;
        buffer_dist::Float64=0.0
    )::POLY_VEC

Buffer exclusion zones by a specified distance.

# Arguments
- `exclusions`: The DataFrame containing exclusion zones.
- `buffer_dist`: The buffer distance.
"""
function buffer_exclusions!(
    exclusion_geometries::POLY_VEC;
    buffer_dist::Float64=0.0
)::POLY_VEC
    exclusion_geometries .= AG.buffer.(exclusion_geometries, buffer_dist)
    return exclusion_geometries
end
function buffer_exclusions(
    exclusion_geometries::POLY_VEC;
    buffer_dist::Float64=0.0
)::POLY_VEC
    copy_vec = copy(exclusion_geometries)
    buffer_exclusions!(copy_vec; buffer_dist=buffer_dist)
    return copy_vec
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
    unionize_overlaps!(geometries::POLY_VEC)::POLY_VEC
    unionize_overlaps(geometries::POLY_VEC)::POLY_VEC
    unionize_overlaps!(exclusions::DataFrame)::DataFrame

Unionize overlapping exclusion zones.

# Arguments
- `geometries`: A vector of geometries (polygons) to unionize.
- `exclusions`: DataFrame containing exclusion zones.
"""
function unionize_overlaps!(geometries::POLY_VEC)::POLY_VEC
    n = length(geometries)

    # Local helper to fix invalid geometries
    function repair(g::IGeometry)::IGeometry{wkbPolygon}
        geom_type = typeof(g)

        if geom_type == IGeometry{wkbLineString} || geom_type == IGeometry{AG.wkbLinearRing}
            pts = [AG.getpoint(g, i) for i in 0:AG.ngeom(g)-1]

            # Ensure closed loop
            ((pts[1][1] != pts[end][1]) || (pts[1][2] != pts[end][2])) && push!(pts, pts[1])

            xs::Vector{Float64} = getindex.(pts, 1)
            ys::Vector{Float64} = getindex.(pts, 2)
            g::IGeometry{wkbPolygon} = AG.createpolygon(xs, ys)
        end

        return g
    end

    # Pre-repair all geometries
    geometries .= repair.(geometries)

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
                union_geom = try
                    AG.union(geom1, geom2)
                catch
                    # Repair/buffer fallback if needed
                    AG.union(AG.buffer(repair(geom1), 0.0), AG.buffer(repair(geom2), 0.0))
                end
                union_geom = repair(union_geom)

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

    filtered = unique(geometries[.!AG.isempty.(geometries)])
    empty!(geometries)
    append!(geometries, filtered)

    return geometries
end
function unionize_overlaps(geometries::POLY_VEC)::POLY_VEC
    geoms_copy = copy(geometries)
    unionize_overlaps!(geoms_copy)
    return geoms_copy
end
function unionize_overlaps!(exclusions::DataFrame)::DataFrame
    geometries = exclusions.geometry
    unique_geometries = unionize_overlaps!(geometries)

    empty!(exclusions)
    append!(exclusions, [(geometry=geom,) for geom in unique_geometries])

    return exclusions
end
