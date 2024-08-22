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

## Notes:

Write out with:

```julia
GDF.write("<path.gpkg>", df, crs=EPSG(7844))
```

Where `crs` can be any valid EPSG code.
"""
function to_dataframe(mp::GI.Wrappers.MultiPolygon)::DataFrame
    return DataFrame(geometry=mp.geom)
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
    to_dataframe

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
