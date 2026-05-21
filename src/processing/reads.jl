
"""
    process_geometry_targets(
        geometries::Vector{IGeometry{wkbPolygon}},
        resolution::Float64=0.0001
    )::Raster{Int}
    process_geometry_targets(
        points::Vector{Point{2,Float64}},
        resolution::Float64=0.0001
    )::Raster{Int}

Read and process target location geometries to generate a rasterized representation.

# Arguments
- `geometries`: A vector of geometries representing target locations.
- `points`: A vector of points representing target locations.
- `resolution`: The resolution for the rasterization process.

# Returns
- A rasterized representation of the target locations.
"""
function process_geometry_targets(
    geometries::POLY_VEC,
    resolution::Float64=0.0001
)::Raster{Int}
    # Compute centroids from the geometries
    target_centroids = AG.centroid.(geometries)
    target_centroid_pts = (
        p -> Point{2,Float64}(p[1:2])).(
        AG.getpoint.(target_centroids, 0)
    )
    return process_geometry_targets(target_centroid_pts, resolution)
end
function process_geometry_targets(
    points::Vector{Point{2,Float64}},
    resolution::Float64=0.0001
)::Raster{Int}
    targets_pts_tuple = [(t[1], t[2]) for t in points]

    lons = getindex.(points, 1)
    lats = getindex.(points, 2)
    extent = Rasters.Extents.Extent(
        X=(minimum(lons) - resolution, maximum(lons) + resolution),
        Y=(minimum(lats) - resolution, maximum(lats) + resolution),
    )

    return Rasters.rasterize(
        last,
        targets_pts_tuple;
        to=extent,
        res=resolution,
        missingval=0,
        fill=1,
    )
end

"""
    raster_to_gdf(raster::Raster{Int,2})::DataFrame

Convert a raster to a GeoDataFrame by extracting the coordinates of the raster cells where
    the value is 1.

# Arguments
- `raster`: The raster to convert.

# Returns
A GeoDataFrame containing the coordinates of the raster cells where the value is 1.
"""
function raster_to_gdf(raster::Raster{Int,2})::DataFrame
    indices::Vector{CartesianIndex{2}} = findall(x -> x == 1, raster)

    lons = raster.dims[1][getindex.(indices, 1)]
    lats = raster.dims[2][getindex.(indices, 2)]

    coords = Point{2,Float64}.(zip(lons, lats))
    return DataFrame(
        ID=1:length(coords),
        geometry=coords
    )
end

"""
    read_and_polygonize_exclusions(
        bathy_fullset_path::String,
        vessel_draft::Float64,
        subset::DataFrame,
        file_name::String,
        output_dir::String="";
        simplify_tol::Float64,
        min_area::Float64=3E-5,
        buffer_dist::Float64=0.0,
    )::DataFrame
    read_and_polygonize_exclusions(
        bathy_fullset_path::String,
        vessel_draft::Float64,
        subset::DataFrame;
        simplify_tol::Float64,
        min_area::Float64=3E-5,
        buffer_dist::Float64=0.0,
    )::DataFrame

Create exclusion zones from environmental constraints.

# Arguments
- `bathy_fullset_path`: The path to the full bathymetry dataset.
- `vessel_draft`: The vessel draft/depth.
- `subset`: The DataFrame containing the study area boundary.
- `file_name`: The name of the output file.
- `output_dir`: The directory to save the output files.

# Returns
The DataFrame containing the exclusion zones.
"""
function read_and_polygonize_exclusions(
    bathy_fullset_path::String,
    vessel_draft::Float64,
    subset::DataFrame,
    file_name::String,
    output_dir::String="";
    simplify_tol::Float64,
    min_area::Float64=3E-5,
    buffer_dist::Float64=0.0,
)::DataFrame
    # ensure output directory exists
    if !isempty(output_dir)
        mkpath(output_dir)
    end
    subset_name::String = split(file_name, "m_")[2]
    exclusion_gpkg_path::String = joinpath(output_dir, "$(file_name).gpkg")
    # TODO: Generalize for all available environmental constraints
    # TODO: Generalize for ms and tender vessels
    # Create exclusion zones from environmental constraints
    if isfile(exclusion_gpkg_path)
        exclusion_zones_df = try
            GDF.read(exclusion_gpkg_path)
        catch e
            @warn "Cached exclusion gpkg unreadable, regenerating: $exclusion_gpkg_path\n$e"
            rm(exclusion_gpkg_path; force=true)
            nothing
        end
        !isnothing(exclusion_zones_df) && return exclusion_zones_df
    end

    exclusion_tif_path::String = joinpath(output_dir, "$(file_name).tif")
    exclusion_zones_bool = if isfile(exclusion_tif_path)
        try
            exclusion_zones_int::Raster = Raster(exclusion_tif_path)
            exclusion_zones_int .!= 0
        catch e
            @warn "Cached exclusion tif unreadable, regenerating: $exclusion_tif_path\n$e"
            rm(exclusion_tif_path; force=true)
            nothing
        end
    end

    if isnothing(exclusion_zones_bool)
        # Load environmental constraints
        # Bathymetry
        bathy_subset_path = joinpath(output_dir, "bathy_subset_$(subset_name).tif")
        bathy_subset = if isfile(bathy_subset_path)
            try
                Raster(bathy_subset_path)
            catch e
                @warn "Cached bathy subset unreadable, regenerating: $bathy_subset_path\n$e"
                rm(bathy_subset_path; force=true)
                nothing
            end
        end

        if isnothing(bathy_subset)
            bathy_dataset = Raster(bathy_fullset_path; lazy=true)
            bathy_subset = read(Rasters.crop(bathy_dataset; to=subset.geom))
            write(bathy_subset_path, bathy_subset; force=true)
        end

        exclusion_zones_bool = create_exclusion_zones(bathy_subset, vessel_draft)
        write(exclusion_tif_path, convert.(Int8, exclusion_zones_bool); force=true)
    end

    exclusion_zones_df = polygonize_binary(exclusion_zones_bool)

    exclusion_zones_df = prepare_exclusion_geoms(
        exclusion_zones_df.geometry;
        min_area=min_area,
        buffer_dist=buffer_dist,
        simplify_tol=simplify_tol
    )
    if !isempty(output_dir) && !isfile(exclusion_gpkg_path)
        GDF.write(exclusion_gpkg_path, exclusion_zones_df; force=true)
    end
    return exclusion_zones_df
end
function read_and_polygonize_exclusions(
    bathy_fullset_path::String,
    vessel_draft::Float64,
    subset::DataFrame;
    simplify_tol::Float64,
    min_area::Float64=3E-5,
    buffer_dist::Float64=0.0,
)::DataFrame
    bathy_dataset::Raster = Raster(bathy_fullset_path; lazy=true)
    bathy_subset::Raster = read(Rasters.crop(bathy_dataset; to=subset.geom))
    exclusion_zones::Raster{Bool} = create_exclusion_zones(bathy_subset, vessel_draft)

    exclusion_zones_df::DataFrame = polygonize_binary(exclusion_zones)
    exclusion_zones_df = prepare_exclusion_geoms(
        exclusion_zones_df.geometry;
        min_area=min_area,
        buffer_dist=buffer_dist,
        simplify_tol=simplify_tol
    )
    return exclusion_zones_df
end

"""
    process_targets(
        target_geometries::Vector{IGeometry{wkbPolygon}},
        subset::DataFrame,
        target_subset_path::String=""
    )::Raster{Int64}
    process_targets(
        target_path::String,
        suitable_threshold::Float64,
        subset::DataFrame,
        target_subset_path::String="",
    )::Raster{Int64}

Reads target locations from target path (.geojson or raster) and returns raster of target
points, applying thresholds and cropping to a target subset.

# Arguments
- `target_geometries`: A vector of geometries representing target locations.
- `target_path`: The path to the target locations file.
- `suitable_threshold`: The threshold for suitable targets.
- `subset`: The DataFrame containing the study area boundary.
- `target_subset_path`: The path to the target subset raster.

# Returns
A rasterized representation of the target locations, cropped to the study area.
"""
function process_targets(
    target_geometries::POLY_VEC,
    subset::DataFrame,
    target_subset_path::String="",
)::Raster{Int64}
    suitable_targets_all = process_geometry_targets(
        target_geometries,
    )
    suitable_targets_subset::Raster{Int} =
        Rasters.crop(suitable_targets_all, to=subset.geom)
    if !isempty(target_subset_path) && !isfile(target_subset_path)
        write(target_subset_path, suitable_targets_subset; force=true)
    end
    return suitable_targets_subset
end
function process_targets(
    target_path::String,
    suitable_threshold::Float64,
    subset::DataFrame,
    target_subset_path::String="",
)::Raster{Int64}
    suitable_targets_all = Raster(target_path; lazy=true)
    suitable_targets_masked = target_threshold(suitable_targets_all, suitable_threshold)

    suitable_targets_subset::Raster{Int} = Rasters.crop(
        suitable_targets_masked,
        to=subset.geom
    )
    if !isempty(target_subset_path) && !isfile(target_subset_path)
        write(target_subset_path, suitable_targets_subset; force=true)
    end
    return suitable_targets_subset
end
