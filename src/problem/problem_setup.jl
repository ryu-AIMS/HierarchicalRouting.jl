
using Glob
using TOML

# TODO: Adapt code to use/incorporate the MSTProblem struct

struct Threshold
    # value::Float64
    # type::Symbol

    min::Float64
    max::Float64

    Threshold(; min::Float64=-Inf, max::Float64=Inf) = new(min, max)
end

struct Constraint
    name::String
    threshold::Threshold
    file_path::String
end

struct Vessel
    exclusion::DataFrame
    capacity::Int64
    number::Int64
    # TODO: add vessel weighting
    # speed::Float64
    # env_constraint::Vector{Constraint}

    # Constructor with default values for capacity and number
    function Vessel(exclusion::DataFrame, capacity::Int64 = typemax(Int64), number::Int64 = 1)
        new(exclusion, capacity, number)
    end
end

struct Problem
    depot::Point{2, Float64}
    mothership::Vessel
    tenders::Vessel
    clusters::Vector{Cluster}
end

function load_problem(target_scenario::String="")

    config = TOML.parsefile(joinpath("src",".config.toml"))

    site_dir = config["data_dir"]["site"]
    target_scenario_dir = config["data_dir"]["target_scenarios"]
    env_constraints_dir = config["data_dir"]["env_constraints"]

    output_dir = config["output_dir"]["path"]

    suitable_threshold = config["parameters"]["suitable_threshold"]
    k = config["parameters"]["k"]
    cluster_tolerance = config["parameters"]["cluster_tolerance"]
    EPSG_code = config["parameters"]["EPSG_code"]

    depot = Point{2, Float64}(config["parameters"]["depot_x"], config["parameters"]["depot_y"])
    n_tenders = config["parameters"]["n_tenders"]
    t_cap = config["parameters"]["t_cap"]

    subset = GDF.read(first(glob("*.gpkg", site_dir)))

    suitable_targets_prefix = target_scenario[1:findlast(".",target_scenario)[1]-1]

    if target_scenario == ""
        target_scenario = first(glob("*", target_scenario_dir))
    else
        suitable_targets_all_path = joinpath(target_scenario_dir, target_scenario)
    end

    #TODO: Move clustering outide of load_problem

    target_subset_path = joinpath(output_dir, "target_subset_$(suitable_targets_prefix).tif")
    target_subset_threshold_path = joinpath(output_dir, "target_subset_$(suitable_targets_prefix)_threshold=$(suitable_threshold).tif")
    clustered_targets_path = joinpath(output_dir, "clustered_$(suitable_targets_prefix)_targets_k=$(k).tif")

    clusts = cluster_targets(
        clustered_targets_path,
        target_subset_threshold_path,
        k, cluster_tolerance,
        suitable_targets_all_path,
        suitable_threshold,
        target_subset_path,
        subset,
        EPSG_code
        )

    # Dynamically discover subfolders in env_constraints_dir
    env_subfolders = readdir(env_constraints_dir)
    env_paths = Dict(subfolder => joinpath(env_constraints_dir, subfolder) for subfolder in env_subfolders)

    # TODO: Use these paths to process all environmental constraints
    for (subfolder, path) in env_paths
        @eval $(Symbol("env_dir_" * subfolder)) = $path
        @eval $(Symbol("rast_path_" * subfolder)) = glob("*.tif", $(Symbol("env_dir_" * subfolder)))
    end
    bathy_fullset_path = rast_path_bathy

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

    t_exclusion_zones_df = process_exclusions(
        bathy_fullset_path,
        config["parameters"]["tend_depth"],
        subset,
        EPSG_code,
        bathy_subset_path,
        joinpath(output_dir, "t_exclusion.gpkg"),
        joinpath(output_dir, "t_exclusion.tif")
    )
    # t_exclusions = unionize_overlaps!(simplify_exclusions!(buffer_exclusions!(t_exclusion_zones_df, 0.1); min_area=100)) #|>  |> simplify_exclusions! |> unionize_overlaps!
    # t_exclusions = t_exclusion_zones_df |> buffer_exclusions! |> simplify_exclusions! |> unionize_overlaps! |> simplify_exclusions! |> unionize_overlaps!

    mothership = Vessel(ms_exclusions)
    tenders = Vessel(t_exclusion_zones_df, t_cap, n_tenders)

    problem = Problem(depot, mothership, tenders, clusts)
    return problem
end
