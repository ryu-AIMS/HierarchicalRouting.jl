
using Glob
using TOML

# struct Threshold
#     # value::Float64
#     # type::Symbol

#     min::Float64
#     max::Float64

#     Threshold(; min::Float64=-Inf, max::Float64=Inf) = new(min, max)
# end

# struct Constraint
#     name::String
#     threshold::Threshold
#     file_path::String
# end

@kwdef struct Vessel
    exclusion::DataFrame
    capacity::Int64 = typemax(Int64)
    number::Int64 = 1
    # TODO: add attributes
    # weight::Float64
    # env_constraint::Vector{Constraint}
end

struct Problem
    target_scenario::String
    depot::Point{2, Float64}
    mothership::Vessel
    tenders::Vessel
end

"""
    load_problem(target_scenario::String="")

Load the problem data from the configuration file and return a Problem object.

# Arguments
- `target_scenario::String`: The name of the target scenario to load. Defaults to an empty string.

# Returns
- `problem::Problem`: The problem object containing the target scenario, depot, mothership, and tenders.
"""
function load_problem(target_scenario::String="")

    config = TOML.parsefile(joinpath("src",".config.toml"))

    site_dir = config["data_dir"]["site"]
    subset = GDF.read(first(glob("*.gpkg", site_dir)))

    target_scenario_dir = config["data_dir"]["target_scenarios"]
    if target_scenario == ""
        target_scenario = first(glob("*", target_scenario_dir))
    end

    depot = Point{2, Float64}(config["parameters"]["depot_x"], config["parameters"]["depot_y"])
    n_tenders = config["parameters"]["n_tenders"]
    t_cap = config["parameters"]["t_cap"]
    EPSG_code = config["parameters"]["EPSG_code"]

    env_constraints_dir = config["data_dir"]["env_constraints"]
    env_subfolders = readdir(env_constraints_dir)
    env_paths = Dict(subfolder => joinpath(env_constraints_dir, subfolder) for subfolder in env_subfolders)

    # TODO: Use these paths to process all environmental constraints
    for (subfolder, path) in env_paths
        @eval $(Symbol("env_dir_" * subfolder)) = $path
        @eval $(Symbol("rast_path_" * subfolder)) = glob("*.tif", $(Symbol("env_dir_" * subfolder)))
    end
    bathy_fullset_path = rast_path_bathy

    output_dir = config["output_dir"]["path"]

    bathy_subset_path = joinpath(output_dir, "bathy_subset.tif")

    # process exclusions
    ms_exclusion_zones_df = process_exclusions(
        bathy_fullset_path,
        config["parameters"]["ms_depth"],
        subset,
        EPSG_code,
        bathy_subset_path,
        joinpath(output_dir, "ms_exclusion.gpkg"),
        joinpath(output_dir, "ms_exclusion.tif")
    )
    ms_exclusions = ms_exclusion_zones_df |> buffer_exclusions! |> simplify_exclusions! |> unionize_overlaps! |> simplify_exclusions! |> unionize_overlaps!

    t_exclusions = process_exclusions(
        bathy_fullset_path,
        config["parameters"]["tend_depth"],
        subset,
        EPSG_code,
        bathy_subset_path,
        joinpath(output_dir, "t_exclusion.gpkg"),
        joinpath(output_dir, "t_exclusion.tif")
    )
    t_exclusions = unionize_overlaps!(
        simplify_exclusions!(
            buffer_exclusions!(
                t_exclusions, buffer_dist=0.1
            ),
            min_area=10,
            simplify_tol=1,
            convex_flag=false
        )
    )

    mothership = Vessel(exclusion = ms_exclusions)
    tenders = Vessel(exclusion = t_exclusions, capacity = t_cap, number = n_tenders)

    problem = Problem(target_scenario, depot, mothership, tenders)
    return problem
end
