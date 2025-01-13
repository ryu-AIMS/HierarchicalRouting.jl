# HierarchicalRouting

HierarchicalRouting is a Julia package that implements metaheuristics to generate
near-optimal routing solutions for a fleet of coordinated vessels, including one mothership
and a given number of satellite vessels (tenders).
The package facilitates planning for scenarios involving multiple clusters of target sites
and environmental constraints.

## Setup

### Initialize project

```julia
]instantiate
```

### Configure problem

Define a `.config.toml` file to specify data paths, output directories, and problem
parameters.
Below is an example:

```toml
[data_dir]
site = "data/site"
env_constraints = "data/env_constraints"
targets = "data/targets"
target_scenarios = "data/targets/scenarios"


[output_dir]
path = "outputs"

[parameters]
EPSG_code = 7844            # Unique coordinate reference system ID for area of interest
depot_x = 50.0              # Longitudinal coordinate
depot_y = 200.0             # Latitudinal coordinate
suitable_threshold = 50.0   # Threshold value for target locations
k = 5                       # Number of clusters to target points
cluster_tolerance = 1.0     # k-means clustering change tolerance at convergence
ms_depth = -10.0            # Mothership draft
tend_depth = -5.0           # Tender draft

n_tenders = 3               # Number of tenders available
t_cap = 2                   # Max number of sites a tender can visit in each deployment sortie
```

By convention, this file is named `.config.toml` (note the leading `.`).

## Quickstart

```julia
using HierarchicalRouting

# Load the problem configuration
# Defaults to the first file in the target scenario directory if no argument is passed.
problem = HierarchicalRouting.load_problem()

# Generate an initial solution
solution_init = HierarchicalRouting.initial_solution(problem)

# Improve solution using simulated annealing
solution_best, z_best = HierarchicalRouting.improve_solution(
    solution_init,
    HierarchicalRouting.simulated_annealing,
    HierarchicalRouting.critical_path,
    HierarchicalRouting.perturb_swap_solution,
    problem.tenders.exclusion
    );
```

## Visualization

To visualize the routing solution, use the `Plot` module. Below is an example to create a
complete plot of clusters, exclusions, and routes:

```julia
using GLMakie

fig = Figure(size = (800, 600))
ax = Axis(fig[1, 1], xlabel = "Longitude", ylabel = "Latitude")

# Add exclusions for the mothership
HierarchicalRouting.Plot.exclusions!(ax, problem.mothership.exclusion, labels = false)

# Add exclusions for tenders
HierarchicalRouting.Plot.exclusions!(ax, problem.tenders.exclusion, labels = false)

# Add clustered points
HierarchicalRouting.Plot.clusters!(
    ax,
    clusters = solution_best.clusters,
    cluster_sequence = solution_best.mothership.cluster_sequence,
    labels=true,
    centers=true,
    cluster_radius = 100
)

# Add mothership routes
HierarchicalRouting.Plot.linestrings!(ax, solution_best.mothership.line_strings; labels =
    true)

# Add tender routes
HierarchicalRouting.Plot.tenders!(ax, solution_best.tenders)
```

### Components

- Clusters: Shows clustered target points and optional labels and/or cluster centers.
- Exclusions: Displays environmental constraints to avoid for the mothership or tenders.
- Routes: Display the mothership and tender routes.

## Development setup

The steps below assume you are in the project root.

### Create a sandbox

Create a sandbox environment and install development dependencies.

```bash
$ mkdir sandbox
$ cd sandbox
$ julia project=.
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

```code
HierarchicalRouting/
├───outputs             # Output files from analyses
└───src                 # Source code for package
    ├───clustering      # Clustering utilities
    ├───optimization    # Optimization heuristics
    ├───plotting        # Visualizations
    ├───problem         # Problem setup and initialization
    ├───processing      # Data processing functions
    └───routing         # Routing algorithms
```

### Data files

Place data files in the `data/` directory, organized as follows:

```
├───data
    ├───env_constraints # Environmental constraint data
    │   ├───bathy       # Bathymetry data
    │   └───waves       # Wave data
    ├───site            # Site data
    └───targets         # Target site scenarios
        └───scenarios
```

## Solution structure

The solution (including problem) structure is as follows:

```
MSTProblem
│
├── target_scenario (String)
├── depot (Point{2, Float64})
├── mothership (Vessel)
└── tenders (Vessel)

Vessel
│
├── exclusion (DataFrame)
├── capacity (Int64)
└── number (Int64)

MSTSolution
│
├── clusters (Vector{Cluster})
├── mothership (MothershipSolution)
└── tenders (Vector{TenderSolution})

Cluster
│
├── id (Int64)
├── centroid (Point{2, Float64})
└── nodes (Vector{Point{2, Float64}})
 
MothershipSolution
│
├── cluster_sequence (DataFrame)
├── route (DataFrame)
├── cost (Float64)
└── line_strings (Vector{LineString{2, Float64}})
 
TenderSolution
│
├── id
├── sorties (Vector{Sortie})
│   ├── nodes (Vector{Point{2, Float64}})
│   └── cost (Float64)
├── start (Point{2, Float64})
├── finish (Point{2, Float64})
└── line_strings (Vector{Vector{LineString{2, Float64}}})

```