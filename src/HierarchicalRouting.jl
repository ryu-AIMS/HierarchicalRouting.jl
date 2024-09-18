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
    simplify_exclusions!(exclusions::DataFrame, buffer_dist=1.0, simplify_tol=0.1)::DataFrame

Simplify exclusions by applying convex hull, and simplification operations to polygons.

# Arguments
- `exclusions::DataFrame`: The DataFrame containing exclusion zones.
- `simplify_tol::Real`: The tolerance value for simplifying the exclusion polygons. Default = 0.1
"""
function simplify_exclusions!(exclusions::DataFrame; simplify_tol::Real=10)#0.7)
    exclusions.geometry .= AG.convexhull.(exclusions.geometry)
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

include("analysis/clustering.jl")
include("analysis/routing_heuristics.jl")
include("analysis/feasible_paths.jl")

include("plotting/plots.jl")


export
    to_multipolygon,
    to_dataframe,
    simplify_exclusions!,
    buffer_exclusions!,
    unionize_overlaps!

export
    extract_subset,
    cluster_targets,
    centroids

export
    create_exclusion_zones,
    nearest_neighbour,
    get_waypoints,
    two_opt,
    return_route_distance,
    two_opt_swap

export
    get_feasible_matrix,
    shortest_feasible_path

export
    plot_polygons,
    plot_mothership_route,
    plot_route_w_exclusions,
    plot_waypoints_and_exclusions,
    plot_waypoints_and_exclusions_with_graph

export EPSG


end
