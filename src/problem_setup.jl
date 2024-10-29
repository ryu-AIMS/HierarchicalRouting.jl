
using Glob
using TOML

# TODO: Adapt code to use/incorporate the MSTProblem struct

struct Threshold
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
    # TODO: add vessel weighting
    speed::Float64
    capacity::Int64
    env_constraint::Vector{Constraint}
    # Vessel(speed::Float64, capacity::Int; env_constraints::Constraints=Constraints()) = new(speed, capacity, env_constraints)
end

struct MSTProblem
    data_file_path::String
    depot::Point{2, Float64}
    ms::Vessel
    tenders::Vector{Vessel}
end

# TODO: Remove hard-coded dir paths!!
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

    subset = GDF.read(first(glob("*.gpkg", site_dir)))

    suitable_targets_prefix = target_scenario[1:findlast(".",target_scenario)[1]-1]

    if target_scenario == ""
        target_scenario = first(glob("*", target_scenario_dir))
    else
        suitable_targets_all_path = joinpath(target_scenario_dir, target_scenario)
    end

    target_subset_path = joinpath(output_dir, "target_subset_$(suitable_targets_prefix).tif")
    target_subset_threshold_path = joinpath(output_dir, "target_subset_$(suitable_targets_prefix)_threshold=$(suitable_threshold).tif")
    clustered_targets_path = joinpath(output_dir, "clustered_$(suitable_targets_prefix)_targets_k=$(k).tif")

    # process targets
    clusters_raster = process_targets(
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

    return clusters_raster, ms_exclusions, t_exclusion_zones_df
end

# TODO Generalize
function process_targets(clustered_targets_path,
    target_subset_threshold_path,
    k, cluster_tolerance,
    suitable_targets_all_path,
    suitable_threshold,
    target_subset_path,
    subset,
    EPSG_code
    )
    if isfile(clustered_targets_path)
        return Raster(clustered_targets_path, mappedcrs=EPSG(EPSG_code))
    else
        # Read deployment locations
        if isfile(target_subset_path)
            suitable_targets_subset = Raster(target_subset_path, mappedcrs=EPSG(EPSG_code))
        else
            if endswith(suitable_targets_all_path, ".geojson")
                suitable_targets_poly = GDF.read(suitable_targets_all_path)
                suitable_targets_centroids = [AG.centroid(polygon) for polygon in suitable_targets_poly.geometry]
                suitable_targets_centroids_pts = [AG.getpoint(centroid,0)[1:2] for centroid in suitable_targets_centroids]

                resolution = 0.0001
                suitable_targets_all = Rasters.rasterize(last, suitable_targets_centroids_pts; res=resolution, missingval=0, fill=1, crs=EPSG(EPSG_code))
                suitable_targets_all = reverse(suitable_targets_all; dims=Y)
            else
                suitable_targets_all = Raster(suitable_targets_all_path, mappedcrs=EPSG(EPSG_code))
                suitable_targets_all = target_threshold(suitable_targets_all, suitable_threshold)
            end
            suitable_targets_subset = crop_to_subset(suitable_targets_all, subset)
            write(target_subset_path, suitable_targets_subset; force=true)
        end
        clustered_targets = cluster_targets(suitable_targets_subset, k; tol=cluster_tolerance)
        write(clustered_targets_path, clustered_targets; force=true)
        return clustered_targets
    end
end

# TODO: Generalize for all available environmental constraints
# TODO: Generalize for ms and tender vessels
function process_exclusions(
    bathy_fullset_path,
    vessel_draft,
    subset,
    EPSG_code,
    bathy_subset_path,
    exclusion_gpkg_path,
    exclusion_tif_path
    )
    # Create exclusion zones from environmental constraints
    if isfile(exclusion_gpkg_path)
        exclusion_zones_df = GDF.read(exclusion_gpkg_path)
    else
        if isfile(exclusion_tif_path)
            exclusion_zones_int = Raster(exclusion_tif_path, mappedcrs=EPSG(EPSG_code))
            exclusion_zones_bool = Raster(exclusion_zones_int .== 0, dims(exclusion_zones_int))
        else
            # Load environmental constraints
            # Bathymetry
            if isfile(bathy_subset_path)
                bathy_subset = Raster(bathy_subset_path, mappedcrs=EPSG(EPSG_code))
            else
                bathy_dataset = Raster(bathy_fullset_path, mappedcrs=EPSG(EPSG_code))
                bathy_subset = crop_to_subset(bathy_dataset, subset)
                write(bathy_subset_path, bathy_subset; force=true)
            end
            exclusion_zones_bool = create_exclusion_zones(bathy_subset, vessel_draft)
            write(exclusion_tif_path, convert.(Int64, exclusion_zones_bool); force=true)
        end

        exclusion_zones_df = exclusion_zones_bool |> to_multipolygon |> to_dataframe
        GDF.write(
            exclusion_gpkg_path,
            exclusion_zones_df;
            crs=EPSG(EPSG_code)
        )
        exclusion_zones_df = GDF.read(exclusion_gpkg_path)
    end
    return exclusion_zones_df
end
