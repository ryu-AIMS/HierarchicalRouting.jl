module HierarchicalRouting

using Base.Threads
using Statistics

using Clustering
using CoordinateTransformations
using GeometryBasics
using Rasters
using DataFrames

import ArchGDAL as AG
import GeoInterface as GI
import GeometryOps as GO
import GeoDataFrames as GDF

using FLoops
using Distances: haversine
using Graphs, SimpleWeightedGraphs

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
    process_exclusions(exclusions::DataFrame, buffer_dist=1.0, simplify_tol=0.1)::DataFrame

Process exclusions by applying buffer, convex hull, and simplification operations to exclusion zones.

# Arguments
- `exclusions::DataFrame`: The DataFrame containing exclusion zones.
- `buffer_dist::Real`: The buffer distance to apply to the exclusion polygons. Default = 1.0
- `simplify_tol::Real`: The tolerance value for simplifying the exclusion polygons. Default = 0.1
"""
function process_exclusions(exclusions::DataFrame; buffer_dist::Real=1.0, simplify_tol::Real=10)#0.7)
    exclusions.geometry .= AG.buffer.(exclusions.geometry, buffer_dist)
    exclusions.geometry .= AG.convexhull.(exclusions.geometry)
    exclusions.geometry .= AG.simplify.(exclusions.geometry, simplify_tol)
    return exclusions
end

function unionize_overlaps(exclusions::DataFrame)
    geometries = exclusions.geometry
    n = length(geometries)
    processed = falses(n)

    for i in 1:n
        if processed[i]
            continue
        end

        geom1 = geometries[i]
        for j in i+1:n
            if processed[j]
                continue
            end

            geom2 = geometries[j]
            if AG.overlaps(geom1, geom2)
                println("Overlap between $i and $j")

                union_geom = AG.union(geom1, geom2)

                exclusions.geometry[i] = union_geom
                exclusions.geometry[j] = union_geom

                for k in 1:n
                    if AG.overlaps(union_geom, geometries[k]) && !processed[k]
                        exclusions.geometry[k] = union_geom
                        processed[k] = true
                    end
                end

            end
        end
    end

    # Remove duplicate unionized geometries
    unique_geometries = unique(exclusions.geometry[.!AG.isempty.(exclusions.geometry)])
    exclusions = DataFrame(geometry = unique_geometries)

    return exclusions
end

include("analysis/clustering.jl")
include("analysis/routing_heuristics.jl")
include("analysis/feasible_paths.jl")

include("plotting/plots.jl")


export
    extract_subset,
    cluster_targets,
    centroids

export
    to_multipolygon,
    to_dataframe,
    process_exclusions,
    unionize_overlaps

export
    create_exclusion_zones,
    nearest_neighbour,
    get_waypoints,
    get_feasible_matrix

export
    plot_polygons,
    plot_mothership_route,
    plot_waypoints_and_exclusions

export EPSG


end
