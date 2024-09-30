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

include("analysis/processing.jl")
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
    target_threshold,
    cluster_targets,
    centroids

export
    create_exclusion_zones,
    nearest_neighbour,
    two_opt,
    return_route_distance,
    two_opt_swap

export
    plot_polygons,
    plot_mothership_route,
    plot_route_w_exclusions,
    plot_waypoints_and_exclusions,
    plot_waypoints_and_exclusions_with_graph

export EPSG


end
