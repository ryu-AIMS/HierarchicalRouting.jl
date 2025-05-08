
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
    raster_to_gdf(raster::Raster{Int, 2})::DataFrame

Convert a raster to a GeoDataFrame by extracting the coordinates of the raster cells where
    the value is 1.

# Arguments
- `raster`: The raster to convert.

# Returns
A GeoDataFrame containing the coordinates of the raster cells where the value is 1.
"""
function raster_to_gdf(raster::Raster{Int, 2})
    indices::Vector{CartesianIndex{2}} = findall(x -> x == 1, raster)

    lons = raster.dims[1][getindex.(indices, 1)]
    lats = raster.dims[2][getindex.(indices, 2)]

    coords = Point{2, Float64}.(zip(lons, lats))
    return DataFrame(
        ID = 1:length(coords),
        geometry = coords
    )
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
    bathy_subset_path::String="",
    exclusion_gpkg_path::String="",
    exclusion_tif_path::String=""
)::DataFrame
    # TODO: Generalize for all available environmental constraints
    # TODO: Generalize for ms and tender vessels
    # Create exclusion zones from environmental constraints
    exclusion_zones_df::DataFrame
    if !isempty(exclusion_gpkg_path) && isfile(exclusion_gpkg_path)
        exclusion_zones_df = GDF.read(exclusion_gpkg_path)
    else
        exclusion_zones_bool::Raster{Bool}
        if !isempty(exclusion_tif_path) && isfile(exclusion_tif_path)
            exclusion_zones_int::Raster = Raster(exclusion_tif_path, mappedcrs=EPSG(EPSG_code))
            exclusion_zones_bool = exclusion_zones_int .!= 0
        else
            # Load environmental constraints
            # Bathymetry
            bathy_subset::Raster
            if !isempty(bathy_subset_path) && isfile(bathy_subset_path)
                bathy_subset = Raster(bathy_subset_path; mappedcrs=EPSG(EPSG_code))
            else
                bathy_dataset::Raster = Raster(bathy_fullset_path; mappedcrs=EPSG(EPSG_code), lazy=true)
                bathy_subset = read(Rasters.crop(bathy_dataset; to=subset.geom))
                if !isempty(bathy_subset_path)
                    write(bathy_subset_path, bathy_subset; force=true)
                end
            end

            exclusion_zones_bool = create_exclusion_zones(bathy_subset, vessel_draft)
            if !isempty(exclusion_tif_path)
                write(exclusion_tif_path, convert.(Int8, exclusion_zones_bool); force=true)
            end
        end

        exclusion_zones_df = polygonize_binary(exclusion_zones_bool)
        if !isempty(exclusion_gpkg_path)
            GDF.write(
                exclusion_gpkg_path,
                exclusion_zones_df;
                crs=EPSG(EPSG_code)
            )
        end
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
    process_targets(
        target_geometries::Vector{AG.IGeometry{AG.wkbPolygon}},
        subset::DataFrame,
        EPSG_code::Int16,
        target_subset_path::String=""
    )::Raster{Int64}
    process_targets(
        target_path::String,
        suitable_threshold::Float64,
        subset::DataFrame,
        EPSG_code::Int16,
        target_subset_path::String="",
    )::Raster{Int64}

Reads target locations from target path (.geojson or raster) and returns raster of target
points, applying thresholds and cropping to a target subset.

# Arguments
- `target_geometries`: A vector of geometries representing target locations.
- `target_path`: The path to the target locations file.
- `suitable_threshold`: The threshold for suitable targets.
- `subset`: The DataFrame containing the study area boundary.
- `EPSG_code`: The EPSG code for the study area.
- `target_subset_path`: The path to the target subset raster.

# Returns
A rasterized representation of the target locations, cropped to the study area.
"""
function process_targets(
    target_geometries::Vector{AG.IGeometry{AG.wkbPolygon}},
    subset::DataFrame,
    EPSG_code::Int16,
    target_subset_path::String="",
)::Raster{Int64}
    suitable_targets_all = process_geometry_targets(
        target_geometries,
        EPSG_code
    )
    suitable_targets_subset::Raster{Int} = Rasters.crop(suitable_targets_all, to=subset.geom)
    if !isempty(target_subset_path) && !isfile(target_subset_path)
        write(target_subset_path, suitable_targets_subset; force=true)
    end
    return suitable_targets_subset
end
function process_targets(
    target_path::String,
    suitable_threshold::Float64,
    subset::DataFrame,
    EPSG_code::Int16,
    target_subset_path::String="",
)::Raster{Int64}
    suitable_targets_all = Raster(target_path; mappedcrs = EPSG(EPSG_code), lazy = true)
    suitable_targets_masked =  target_threshold(suitable_targets_all, suitable_threshold)

    suitable_targets_subset::Raster{Int} = Rasters.crop(
        suitable_targets_masked,
        to=subset.geom
    )
    if !isempty(target_subset_path) && !isfile(target_subset_path)
        write(target_subset_path, suitable_targets_subset; force=true)
    end
    return suitable_targets_subset
end
