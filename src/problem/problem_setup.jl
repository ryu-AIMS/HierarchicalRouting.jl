
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
    weighting::Float64 = 1.0
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
The problem object containing the target scenario, depot, mothership, and tenders.
"""
function load_problem(target_scenario::String="")

    config = TOML.parsefile(joinpath("src",".config.toml"))
    site_dir = config["data_dir"]["site"]
    depot = Point{2, Float64}(config["parameters"]["depot_x"], config["parameters"]["depot_y"])
    n_tenders = config["parameters"]["n_tenders"]
    t_cap = config["parameters"]["t_cap"]
    EPSG_code = config["parameters"]["EPSG_code"]
    target_scenario = target_scenario == "" ? first(glob("*", config["data_dir"]["target_scenarios"])) : target_scenario

    env_constraints_dir = config["data_dir"]["env_constraints"]
    env_subfolders = readdir(env_constraints_dir)
    env_paths = Dict(subfolder => joinpath(env_constraints_dir, subfolder) for subfolder in env_subfolders)
    env_data = Dict(
        subfolder => (
            env_dir = path,
            rast_file = isempty(glob("*.tif", path)) ? nothing : first(glob("*.tif", path))
        ) for (subfolder, path) in env_paths
    )
    # TODO: Process all environmental constraints to create single cumulative exclusion zone

    output_dir = config["output_dir"]["path"]

    subset = GDF.read(first(glob("*.gpkg", site_dir)))
    bathy_subset_path = joinpath(output_dir, "bathy_subset.tif")

    # Process exclusions
    if !(config["DEBUG"]["debug_mode"])
        ms_exclusion_zones_df = read_and_polygonize_exclusions(
            env_data["bathy"].rast_file,
            config["parameters"]["depth_ms"],
            subset,
            EPSG_code,
            bathy_subset_path,
            joinpath(output_dir, "ms_exclusion.gpkg"),
            joinpath(output_dir, "ms_exclusion.tif")
        )
    else
        ms_exclusion_zones_df = read_and_polygonize_exclusions(
            env_data["bathy"].rast_file,
            config["parameters"]["depth_ms"],
            subset,
            EPSG_code
        )
    end
    ms_exclusions = ms_exclusion_zones_df |> filter_and_simplify_exclusions! |> buffer_exclusions! |> unionize_overlaps! |> filter_and_simplify_exclusions! |> unionize_overlaps!

    if !(config["DEBUG"]["debug_mode"])
        t_exclusions = read_and_polygonize_exclusions(
            env_data["bathy"].rast_file,
            config["parameters"]["depth_t"],
            subset,
            EPSG_code,
            bathy_subset_path,
            joinpath(output_dir, "t_exclusion.gpkg"),
            joinpath(output_dir, "t_exclusion.tif")
        )
    else
        t_exclusions = read_and_polygonize_exclusions(
            env_data["bathy"].rast_file,
            config["parameters"]["depth_t"],
            subset,
            EPSG_code
        )
    end

    #! Functions NOT modiying in place
    t_exclusions = filter_and_simplify_exclusions!(
        t_exclusions, min_area=1E-7, simplify_tol=5E-5
    )
    t_exclusions = buffer_exclusions!(t_exclusions, buffer_dist=0.0)
    t_exclusions = unionize_overlaps!(t_exclusions)
    #? Redundant 2nd call to simplify_exclusions!()
    t_exclusions = filter_and_simplify_exclusions!(
        t_exclusions, min_area=1E-7, simplify_tol=5E-5
    )

    mothership = Vessel(exclusion = ms_exclusions, weighting = config["parameters"]["weight_ms"])
    tenders = Vessel(exclusion = t_exclusions, capacity = t_cap, number = n_tenders, weighting = config["parameters"]["weight_t"])

    problem = Problem(target_scenario, depot, mothership, tenders)
    return problem
end
