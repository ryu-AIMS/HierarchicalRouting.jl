
"""
    process_targets(
        clustered_targets_path::String,
        k::Int8,
        cluster_tolerance::Float64,
        suitable_targets_all_path::String,
        suitable_threshold::Float64,
        target_subset_path::String,
        subset::DataFrame,
        EPSG_code::Int16
    )::Raster{Int64}

Generate a clustered targets raster by reading in the suitable target location data,
applying thresholds and cropping to a target subset, and then clustering.

# Arguments
- `clustered_targets_path`: The path to the clustered targets raster.
- `k`: The number of clusters.
- `cluster_tolerance`: The cluster tolerance.
- `suitable_targets_all_path`: The path to the suitable targets dataset.
- `suitable_threshold`: The suitable targets threshold.
- `target_subset_path`: The path to the target subset raster.
- `subset`: The DataFrame containing the study area boundary.
- `EPSG_code`: The EPSG code for the study area.

# Returns
The clustered targets raster, classified by cluster ID number.
"""
function process_targets(
    clustered_targets_path::String,
    k::Int8,
    cluster_tolerance::Float64,
    suitable_targets_all_path::String,
    suitable_threshold::Float64,
    target_subset_path::String,
    subset::DataFrame,
    EPSG_code::Int16
)::Raster{Int64}
        if endswith(suitable_targets_all_path, ".geojson")
            resolution = 0.0001
            suitable_targets_poly::DataFrame = GDF.read(suitable_targets_all_path)
            suitable_targets_centroids::Vector{AG.IGeometry{AG.wkbPoint}} =
                AG.centroid.(suitable_targets_poly.geometry)
            suitable_targets_centroids_pts::Vector{NTuple{2, Float64}} =
                (x -> x[1:2]).(AG.getpoint.(suitable_targets_centroids, 0))
            suitable_targets_all = Rasters.rasterize(
                last, suitable_targets_centroids_pts;
                res=resolution, missingval=0, fill=1, crs=EPSG(EPSG_code)
            )
    else
            suitable_targets_all = Raster(
                suitable_targets_all_path, mappedcrs=EPSG(EPSG_code)
            )
            suitable_targets_all = target_threshold(suitable_targets_all, suitable_threshold)
        end
            suitable_targets_subset::Raster = Rasters.crop(suitable_targets_all, to=subset.geom)
    if !isfile(target_subset_path)
        write(target_subset_path, suitable_targets_subset; force=true)
    end

        clustered_targets::Raster{Int64, 2} = apply_kmeans_clustering(
            suitable_targets_subset, k; tol=cluster_tolerance
        )
    write(clustered_targets_path, clustered_targets; force=true)
    return clustered_targets
end

"""
    read_and_polygonize_exclusions(
        bathy_fullset_path::String,
        vessel_draft::Float64,
        subset::DataFrame,
        EPSG_code::Int16,
        bathy_subset_path::String,
        exclusion_gpkg_path::String,
        exclusion_tif_path::String
    )::DataFrame

Create exclusion zones from environmental constraints.

# Arguments
- `bathy_fullset_path`: The path to the full bathymetry dataset.
- `vessel_draft`: The vessel draft/depth.
- `subset`: The DataFrame containing the study area boundary.
- `EPSG_code`: The EPSG code for the study area.
- `bathy_subset_path`: The path to the subset bathymetry dataset.
- `exclusion_gpkg_path`: The path to the exclusion zones GeoPackage.
- `exclusion_tif_path`: The path to the exclusion zones raster.

# Returns
The DataFrame containing the exclusion zones.
"""
function read_and_polygonize_exclusions(
    bathy_fullset_path::String,
    vessel_draft::Float64,
    subset::DataFrame,
    EPSG_code::Int16,
    bathy_subset_path::String,
    exclusion_gpkg_path::String,
    exclusion_tif_path::String
)::DataFrame
    # TODO: Generalize for all available environmental constraints
    # TODO: Generalize for ms and tender vessels
    # Create exclusion zones from environmental constraints
    exclusion_zones_df::DataFrame
    if isfile(exclusion_gpkg_path)
        exclusion_zones_df = GDF.read(exclusion_gpkg_path)
    else
        exclusion_zones_bool::Raster{Bool}
        if isfile(exclusion_tif_path)
            exclusion_zones_int::Raster = Raster(exclusion_tif_path, mappedcrs=EPSG(EPSG_code))
            exclusion_zones_bool = exclusion_zones_int .!= 0
        else
            # Load environmental constraints
            # Bathymetry
            bathy_subset::Raster
            if isfile(bathy_subset_path)
                bathy_subset = Raster(bathy_subset_path; mappedcrs=EPSG(EPSG_code))
            else
                bathy_dataset::Raster = Raster(bathy_fullset_path; mappedcrs=EPSG(EPSG_code), lazy=true)
                bathy_subset = read(Rasters.crop(bathy_dataset; to=subset.geom))
                write(bathy_subset_path, bathy_subset; force=true)
            end

            exclusion_zones_bool = create_exclusion_zones(bathy_subset, vessel_draft)
            write(exclusion_tif_path, convert.(Int8, exclusion_zones_bool); force=true)
        end

        exclusion_zones_df = polygonize_binary(exclusion_zones_bool)
        GDF.write(
            exclusion_gpkg_path,
            exclusion_zones_df;
            crs=EPSG(EPSG_code)
        )
    end

    return exclusion_zones_df
end
function read_and_polygonize_exclusions(
    bathy_fullset_path::String,
    vessel_draft::Float64,
    subset::DataFrame,
    EPSG_code::Int16
)::DataFrame
    bathy_dataset::Raster = Raster(bathy_fullset_path; mappedcrs=EPSG(EPSG_code), lazy=true)
    bathy_subset::Raster = read(Rasters.crop(bathy_dataset; to=subset.geom))
    exclusion_zones::Raster{Bool} = create_exclusion_zones(bathy_subset, vessel_draft)

    return polygonize_binary(exclusion_zones)
end

"""
    cluster_problem(problem::Problem)::Vector{Cluster}

Read and process problem data to generate an initial solution.

# Arguments
- `problem`: The problem data.

# Returns
Vector of clustered locations.
"""
function cluster_problem(problem::Problem)::Vector{Cluster}
    config = TOML.parsefile(joinpath("src", ".config.toml"))

    suitable_threshold::Float64 = config["parameters"]["suitable_threshold"]
    k::Int8 = config["parameters"]["k"]
    cluster_tolerance::Float64 = config["parameters"]["cluster_tolerance"]
    EPSG_code::Int16 = config["parameters"]["EPSG_code"]

    site_dir::String = config["data_dir"]["site"]
    subset::DataFrame = GDF.read(first(glob("*.gpkg", site_dir)))

    output_dir::String = config["output_dir"]["path"]
    target_scenario_dir = config["data_dir"]["target_scenarios"]
    scenario_name::String = problem.scenario_name
    clustered_targets_path::String = joinpath(
        output_dir,
        "clustered_$(scenario_name)_targets_k=$(k).tif"
    )
    target_subset_path::String = joinpath(output_dir, "target_subset_$(scenario_name).tif")
    suitable_targets_all_path::String = joinpath(target_scenario_dir, scenario_name*".geojson")

    clusters::Vector{Cluster} = generate_target_clusters(
        clustered_targets_path,
        k,
        cluster_tolerance,
        suitable_targets_all_path,
        suitable_threshold,
        target_subset_path,
        subset,
        EPSG_code
    )

    return clusters
end
