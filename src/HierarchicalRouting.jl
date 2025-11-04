module HierarchicalRouting

using Base
using Statistics

using CoordinateTransformations
using GeometryBasics
using Rasters
using DataFrames

import ArchGDAL as AG
import ArchGDAL: IGeometry, wkbPolygon, wkbLineString, wkbPoint, wkbMultiPolygon
import GeoInterface as GI
import GeometryOps as GO
import GeoDataFrames as GDF

using Distances
using Graphs, SimpleWeightedGraphs

const POLY_VEC = Vector{IGeometry{wkbPolygon}}

include("problem/problem_setup.jl")

include("processing/processing.jl")
include("processing/data_types.jl")
include("processing/reads.jl")
include("processing/spatial_operations.jl")

include("clustering/clustering.jl")

include("routing/routing_heuristics.jl")
include("routing/waypoints.jl")
include("routing/feasible_paths.jl")
include("routing/polygon_traversal.jl")
include("optimization/metric_calcs.jl")
include("optimization/solution_assessment.jl")

include("problem/solve_problem.jl")
include("problem/disturbance_events.jl")

include("plot/plot.jl")
include("show.jl")

include("processing/write_files.jl")

export Plot

export load_problem, initial_solution, improve_solution

export cluster_problem, nearest_neighbour, two_opt, tender_sequential_nearest_neighbour

export generate_cluster_df, optimize_mothership_route, disturb_clusters

end
