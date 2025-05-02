
"""
    process_geometry_targets(
        geometries::Vector{AG.IGeometry{AG.wkbPolygon}},
        EPSG_code::Int16,
        resolution::Float64 = 0.0001
    )::Raster{Int}

Read and process target location geometries to generate a rasterized representation.

# Arguments
- `geometries`: A vector of geometries representing target locations.
- `EPSG_code`: The EPSG code for the coordinate reference system.
- `resolution`: The resolution for the rasterization process.

# Returns
- A rasterized representation of the target locations.
"""
function process_geometry_targets(
    geometries::Vector{AG.IGeometry{AG.wkbPolygon}},
    EPSG_code::Int16,
    resolution::Float64 = 0.0001
)::Raster{Int}
    # Compute centroids from the geometries
    target_centroids = AG.centroid.(geometries)
    target_centroid_pts = (
        p -> Point{2,Float64}(p[1:2])).(
        AG.getpoint.(target_centroids, 0)
    )
    targets_pts_tuple = [(t[1], t[2]) for t in target_centroid_pts]

    return Rasters.rasterize(
        last,
        targets_pts_tuple;
        res = resolution,
        missingval = 0,
        fill = 1,
        crs = EPSG(EPSG_code)
    )
end

"""
    process_raster_targets(
    targets_path::String,
    EPSG_code::Int16,
    suitable_threshold::Float64
)::Raster

Read and mask target locations from a raster file.

# Arguments
- `targets_path`: The path to the raster file containing target locations.
- `EPSG_code`: The EPSG code for the coordinate reference system.
- `suitable_threshold`: The threshold value for suitable targets.

# Returns
- A rasterized representation of the target locations.
"""
function process_raster_targets(
    targets_path::String,
    EPSG_code::Int16,
    suitable_threshold::Float64
)::Raster
    # TODO: fill in with wave data!!
    suitable_targets_all = Raster(targets_path; mappedcrs = EPSG(EPSG_code), lazy = true)
    return target_threshold(suitable_targets_all, suitable_threshold)
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
    config = TOML.parsefile(".config.toml")

    suitable_threshold::Float64 = config["parameters"]["suitable_threshold"]
    k::Int8 = config["parameters"]["k"]
    cluster_tolerance::Float64 = config["parameters"]["cluster_tolerance"]
    EPSG_code::Int16 = config["parameters"]["EPSG_code"]

    site_dir::String = config["data_dir"]["site"]
    subset::DataFrame = GDF.read(first(glob("*.gpkg", site_dir)))
    suitable_targets_filename = splitext(basename(problem.targets.path))[1]

    output_dir::String = config["output_dir"]["path"]

    clustered_targets_path::String = joinpath(
        output_dir,
        "clustered_$(suitable_targets_filename)_targets_k=$(k).tif"
    )
    target_subset_path::String = joinpath(
        output_dir,
        "target_subset_$(suitable_targets_filename).tif"
    )

    clusters::Vector{Cluster} = generate_target_clusters(
        clustered_targets_path,
        k,
        cluster_tolerance,
        problem.targets,
        suitable_threshold,
        target_subset_path,
        subset,
        EPSG_code,
        problem.depot,
        problem.tenders.exclusion,
    )

    return clusters
end

"""
    process_targets(
        target_geometries::Vector{AG.IGeometry{AG.wkbPolygon}},
        target_subset_path::String,
        subset::DataFrame,
        EPSG_code::Int16
    )::Raster{Int64}
    process_targets(
        target_path::String,
        suitable_threshold::Float64,
        target_subset_path::String,
        subset::DataFrame,
        EPSG_code::Int16
    )::Raster{Int64}

Reads target locations from target path (.geojson or raster) and returns raster of target
points, applying thresholds and cropping to a target subset.

# Arguments
- `target_geometries`: A vector of geometries representing target locations.
- `target_path`: The path to the target locations file.
- `suitable_threshold`: The threshold for suitable targets.
- `target_subset_path`: The path to the target subset raster.
- `subset`: The DataFrame containing the study area boundary.
- `EPSG_code`: The EPSG code for the study area.

# Returns
- A rasterized representation of the target locations.
"""
function process_targets(
    target_geometries::Vector{AG.IGeometry{AG.wkbPolygon}},
    target_subset_path::String,
    subset::DataFrame,
    EPSG_code::Int16,
)::Raster{Int64}
    #TODO: Remove process_geometry_targets() function call and process directly
    suitable_targets_all = process_geometry_targets(
        target_geometries,
        EPSG_code
    )
    suitable_targets_subset::Raster{Int} = Rasters.crop(suitable_targets_all, to=subset.geom)
    if !isfile(target_subset_path)
        write(target_subset_path, suitable_targets_subset; force=true)
    end
    return suitable_targets_subset
end
function process_targets(
    target_path::String,
    suitable_threshold::Float64,
    target_subset_path::String,
    subset::DataFrame,
    EPSG_code::Int16
)::Raster{Int64}
    #TODO: Remove process_raster_targets() function call and process directly
    suitable_targets_all = process_raster_targets(
        target_path,
        EPSG_code,
        suitable_threshold
    )
    suitable_targets_subset::Raster{Int} = Rasters.crop(suitable_targets_all, to=subset.geom)
    if !isfile(target_subset_path)
        write(target_subset_path, suitable_targets_subset; force=true)
    end
    return suitable_targets_subset
end
