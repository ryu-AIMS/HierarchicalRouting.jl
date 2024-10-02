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

include("problem_setup.jl")

include("analysis/processing.jl")
include("analysis/clustering.jl")
include("analysis/routing_heuristics.jl")
include("analysis/feasible_paths.jl")

include("plotting/plots.jl")


export
    load_problem


export
    create_clusters

export
    create_exclusion_zones,
    nearest_neighbour,
    two_opt,
    return_route_distance,
    two_opt_swap,
    tender_sequential_nearest_neighbour

export
    plot_polygons,
    plot_mothership_route,
    plot_route_w_exclusions,
    plot_centroids_and_exclusions,
    plot_waypoints_and_exclusions,
    plot_waypoints_and_exclusions_with_graph

export EPSG


end
