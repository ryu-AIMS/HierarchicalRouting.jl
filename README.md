# HierarchicalRouting

HierarchicalRouting is a Julia package that implements metaheuristics to generate
near‐optimal routing solutions for a fleet of coordinated vessels, including one mothership
and a given number of satellite vessels (tenders).
The package supports planning for scenarios with multiple multiple clusters of target sites
and environmental constraints (exclusion zones).

## Setup

### Initialize project

In the Julia REPL:
```julia
] instantiate
```

### Configure problem

Define a configuration file named `.config.toml` in the project root.
This file should now only contain the debug configuration.
Below is an example:

```toml
[DEBUG]
debug_mode = false          # Set to false to disable result caching for debugging purposes
```

By convention, this file is named `.config.toml` (note the leading `.`).

All other problem parameters (data paths, vessel capacities, clustering tolerance, etc.)
must be provided by the user directly as arguments to the `load_problem()` and
`solve_problem()` functions.

## Quickstart

Below is an example of how to load a problem, generate an initial solution, and improve the
solution.

```julia
using HierarchicalRouting

# Load the problem configuration
# Defaults to the first file in the target scenario directory if no argument is passed.
problem = load_problem(
    "data/targets/scenarios/output_slopes_3-10m.geojson",           # target_path
    "data/site/Moore_2024-02-14b_v060_rc1.gpkg",                    # subset_path
    "data/env_constraints/bathy/Cairns-Cooktown_bathy.tif",         # bathy_path
    "data/env_disturbances/waves/output_slope_zs_Hs_Tp.geojson",    # wave_disturbance_path
    Point{2, Float64}(146.175, -16.84),                             # depot
    -10.0,                                                          # draft_ms
    -5.0,                                                           # draft_t
    Float16(5.0),                                                   # weight_ms
    Float16(2.0),                                                   # weight_t
    Int8(3),                                                        # n_tenders
    Int16(2);                                                       # t_cap
);

# Generate an initial solution
solution_init = initial_solution(
    problem,                                                        # problem
    Int8(5),                                                        # num_clusters
    Float64(5E-5),                                                  # cluster_tolerance
    Set((2, 4))                                                     # disturbance_clusters
);

# Improve solution using simulated annealing
solution_best, z_best = improve_solution(
    solution_init,                                                  # solution_init
    HierarchicalRouting.simulated_annealing,                        # opt_function
    HierarchicalRouting.critical_path,                              # objective_function
    HierarchicalRouting.perturb_swap_solution,                      # perturb_function
    problem.tenders.exclusion;                                      # exclusion
);
```

## Visualization

To visualize the routing solution, use GeoMakie and the `Plot` module.
Below are examples to create:
- one complete plot of final clusters, exclusions, and routes; and
- a series of plots for each sequential routing plan pre-deployment, and at each(assuming 2)
disturbance event.

```julia
using GeoMakie

fig = Figure(size=(750, 880))
ax = Axis(fig[1, 1], xlabel="Longitude", ylabel="Latitude");
# Add exclusions for the mothership
HierarchicalRouting.Plot.exclusions!(ax, problem.mothership.exclusion, labels=true);
# Add exclusions for tenders
HierarchicalRouting.Plot.exclusions!(ax, problem.tenders.exclusion, labels=true);
# Add clustered points
HierarchicalRouting.Plot.clusters!(ax, clusters=solution_best.cluster_sets[end]);
# Add mothership routes
HierarchicalRouting.Plot.linestrings!(ax, solution_best.mothership_routes[end].route);
# Add tender routes
HierarchicalRouting.Plot.tenders!(ax, solution_best.tenders);
```

```julia
using GeoMakie

fig = Figure(size=(1650, 600))
ax1, ax2, ax3 =
    Axis(fig[1, 1], xlabel="Longitude", ylabel="Latitude"),
    Axis(fig[1, 2], xlabel="Longitude"),
    Axis(fig[1, 3], xlabel="Longitude");
# Add exclusions for tenders
HierarchicalRouting.Plot.exclusions!.(
    [ax1, ax2, ax3],
    [problem.tenders.exclusion],
    labels=false
);
# Add clustered points
HierarchicalRouting.Plot.clusters!(
    ax1,
    clusters=solution_best.cluster_sets[end],
    labels=true,
    centers=false,
    nodes= true,
    cluster_radius=0.025
);
HierarchicalRouting.Plot.clusters!(
    ax2,
    clusters=solution_best.cluster_sets[end],
    labels=true,
    centers=false,
    nodes= true,
    cluster_radius=0.025
);
HierarchicalRouting.Plot.clusters!(
    ax3,
    clusters=solution_best.cluster_sets[end],
    labels=true,
    centers=false,
    nodes= true,
    cluster_radius=0.025
);
# Add mothership routes
HierarchicalRouting.Plot.linestrings!.(
    [ax1, ax2, ax3],
    [solution_best.mothership_routes[1].route,
     solution_best.mothership_routes[2].route,
     solution_best.mothership_routes[3].route],
    labels=true,
    color=:black
);
# Add tender routes
HierarchicalRouting.Plot.tenders!(
    ax1,
    [solution_best.tenders[1]],
    5
);
HierarchicalRouting.Plot.tenders!(
    ax2,
    solution_best.tenders[1:3],
    5
);
HierarchicalRouting.Plot.tenders!(
    ax3,
    solution_best.tenders
);
```

## Development setup

The steps below assume you are in the project root.

### Create a sandbox

Create a sandbox environment and install development dependencies.

```bash
mkdir sandbox
cd sandbox

julia project=.
(sandbox) julia> ]add Revise Infiltrate
(sandbox) julia> ]dev ..
```

### Run a development script

- Copy the quickstart to a file (e.g., `dev_routing.jl`) and save to the sandbox directory.
- Create the `.config.toml` file and save to the sandbox directory.
- Start the Julia REPL at project root:
(Assuming VS Code is configured to default to the sandbox environment)

```julia
;cd sandbox
include("dev_routing.jl")
```

## Project structure

The project structure is as follows:

```
HierarchicalRouting/
├───outputs             # Output files from analyses
├───src                 # Source code for package
│   ├───clustering      # Clustering utilities
│   ├───optimization    # Optimization heuristics
│   ├───plotting        # Visualizations
│   ├───problem         # Problem setup and initialization
│   ├───processing      # Data processing functions
│   └───routing         # Routing algorithms
└───test
```

## Solution structure

The solution (including problem) structure is as follows:

```
Problem
│
├── depot (Point{2, Float64})
├── targets (Targets)
├── mothership (Vessel)
└── tenders (Vessel)

Targets
│
├── points (DataFrame)
├── path (String)
└── disturbance_gdf (DataFrame)

Vessel
│
├── exclusion (DataFrame)
├── capacity (Int16)
├── number (Int8)
└── weighting (Float16)

MSTSolution
│
├── cluster_sets (Vector{Vector{Cluster}})
├── mothership_routes (Vector{MothershipSolution})
└── tenders (Vector{TenderSolution})

Cluster
│
├── id (Int64)
├── centroid (Point{2, Float64})
└── nodes (Vector{Point{2, Float64}})
 
MothershipSolution
│
├── cluster_sequence (DataFrame)
└── route (Route)

Route
│
├── nodes (Vector{Point{2, Float64}})
├── dist_matrix (Matrix{Float64})
└── line_strings (Vector{LineString{2, Float64}})
 
TenderSolution
│
├── id
├── start (Point{2, Float64})
├── finish (Point{2, Float64})
├── sorties (Vector{Route})
└── dist_matrix (Matrix{Float64})

```