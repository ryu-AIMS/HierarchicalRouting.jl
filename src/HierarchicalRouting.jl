module HierarchicalRouting

using Statistics

import ArchGDAL as AG
import GeoInterface as GI
import GeometryOps as GO
using GeometryBasics
using CoordinateTransformations

using Rasters
using DataFrames
import GeoDataFrames as GDF

using Clustering
using Distances: haversine
using Graphs, SimpleWeightedGraphs

"""
    to_multipolygon(raster::Raster{T, 2}) where {T<:Union{Bool,Int16}}

Convert raster to multipolygons.
Invert vertical axis.
"""
function to_multipolygon(
    raster::Raster{T, 2}
)::GI.Wrappers.MultiPolygon where {T<:Union{Bool,Int16}}
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
function process_exclusions(exclusions::DataFrame; buffer_dist::Real=1.0, simplify_tol::Real=0.7)
    exclusions.geometry .= AG.buffer.(exclusions.geometry, buffer_dist)
    exclusions.geometry .= AG.convexhull.(exclusions.geometry)
    exclusions.geometry .= AG.simplify.(exclusions.geometry, simplify_tol)
    return
end


include("analysis/clustering.jl")
include("analysis/routing.jl")

include("plotting/plots.jl")


export
    extract_subset,
    cluster_targets,
    centroids

export
    to_multipolygon,
    to_dataframe,
    process_exclusions

export
    create_exclusion_zones,
    nearest_neighbour,
    get_waypoints,
    get_feasible_matrix

export
    plot_polygons,
    plot_mothership_route

export EPSG


end
