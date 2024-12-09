module HierarchicalRouting

using Base
using Statistics

using CoordinateTransformations
using GeometryBasics
using Rasters
using DataFrames

import ArchGDAL as AG
import GeoInterface as GI
import GeometryOps as GO
import GeoDataFrames as GDF

import Distances: euclidean # haversine
using Graphs, SimpleWeightedGraphs

include("problem/problem_setup.jl")

include("processing/processing.jl")
include("processing/data_types.jl")
include("processing/reads.jl")
include("processing/spatial_operations.jl")

include("clustering/clustering.jl")

include("routing/routing_heuristics.jl")
include("routing/feasible_paths.jl")
include("optimization/metric_calcs.jl")
include("optimization/solution_assessment.jl")

include("problem/solve_problem.jl")

include("plotting/plots.jl")

export Plot


end
