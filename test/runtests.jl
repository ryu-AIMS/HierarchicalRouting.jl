using Test
using HierarchicalRouting

include("setup.jl")

const (TEST_PROBLEM, TEST_SOLUTION) = initialize_instances()

include("problem_initialization.jl")
include("point_feasibility.jl")
include("route_consistency.jl")
