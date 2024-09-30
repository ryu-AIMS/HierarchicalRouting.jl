
function load_problem(path::String="")
    EPSG_code = 7844
    depot = Point{2, Float64}(0.0, 0.0)
    suitable_threshold = 50.0
    k = 4
    ms_depth = -10.0

    # target site area location
    subset = define_site() #GDF.read("data/site/Moore_2024-02-14b_v060_rc1.gpkg")

    suitable_targets_all_path = joinpath(path, "data/targets/Cairns-Cooktown_suitable_slopes.tif")

    target_subset_path = joinpath(path, "outputs/target_subset.tif")
    target_subset_threshold_path = joinpath(path, "outputs/target_subset_threshold=$(suitable_threshold).tif")
    clustered_targets_path = joinpath(path, "outputs/clustered_threshold=$(suitable_threshold)_targets_k=$(k).tif")

    # process targets
    clustered_targets = process_targets(
        clustered_targets_path,
        target_subset_threshold_path,
        k,
        suitable_targets_all_path,
        suitable_threshold,
        target_subset_path,
        subset,
        EPSG_code
    )

    bathy_fullset_path = joinpath(path, "data/bathy/Cairns-Cooktown/Cairns-Cooktown_bathy.tif")

    ms_exclusion_gpkg_path = joinpath(path, "outputs/ms_exclusion.gpkg")
    ms_exclusion_tif_path = joinpath(path, "outputs/ms_exclusion.tif")
    bathy_subset_path = joinpath(path, "outputs/bathy_subset.tif")

    # process exclusions
    ms_exclusion_zones_df = process_exclusions(
        ms_exclusion_gpkg_path,
        ms_exclusion_tif_path,
        bathy_subset_path,
        bathy_fullset_path,
        EPSG_code,
        ms_depth,
        subset
    )
    ms_exclusions = ms_exclusion_zones_df |> buffer_exclusions! |> simplify_exclusions! |> unionize_overlaps! |> simplify_exclusions! |> unionize_overlaps!

    return clustered_targets, ms_exclusions#, tender_exclusions
end

function define_site()
    files = readdir("data/site"; join=true)
    gpkg_file = files[endswith.(files, ".gpkg")]

    if length(gpkg_file) == 1
        return GDF.read(gpkg_file[1])
    else
        error("Expected .gpkg site area file in 'data/site', but found $(length(gpkg_file)).")
    end
end

# TODO Generalize
function process_targets(clustered_targets_path,
    target_subset_threshold_path,
    k,
    suitable_targets_all_path,
    suitable_threshold,
    target_subset_path,
    subset,
    EPSG_code
    )
    if isfile(clustered_targets_path)
        return Raster(clustered_targets_path, mappedcrs=EPSG(EPSG_code))
    else
        if isfile(target_subset_threshold_path)
            suitable_targets_subset_threshold = Raster(target_subset_threshold_path, mappedcrs=EPSG(EPSG_code))
        else
            # Read deployment locations
            if isfile(target_subset_path)
                suitable_targets_subset = Raster(target_subset_path, mappedcrs=EPSG(EPSG_code))
            else
                suitable_targets_all = Raster(suitable_targets_all_path, mappedcrs=EPSG(EPSG_code))
                suitable_targets_subset = extract_subset(suitable_targets_all, subset)
                write(target_subset_path, suitable_targets_subset; force=true)
            end
            suitable_targets_subset_threshold = target_threshold(suitable_targets_subset, suitable_threshold)
            write(target_subset_threshold_path, suitable_targets_subset_threshold; force=true)
        end

        clustered_targets = cluster_targets(suitable_targets_subset_threshold, k)
        write(clustered_targets_path, clustered_targets; force=true)
        return clustered_targets
    end
end

# TODO: Generalize for all available environmental constraints
# TODO: Generalize for ms and tender vessels
function process_exclusions(ms_exclusion_gpkg_path,
    ms_exclusion_tif_path,
    bathy_subset_path,
    bathy_fullset_path,
    EPSG_code,
    ms_depth,
    subset
    )
    # Create exclusion zones from environmental constraints
    if isfile(ms_exclusion_gpkg_path)
        ms_exclusion_zones_df = GDF.read(ms_exclusion_gpkg_path)
    else
        if isfile(ms_exclusion_tif_path)
            ms_exclusion_zones_int = Raster(ms_exclusion_tif_path, mappedcrs=EPSG(EPSG_code))
            ms_exclusion_zones_bool = Raster(ms_exclusion_zones_int .== 0, dims(ms_exclusion_zones_int))
        else
            # Load environmental constraints
            # Bathymetry
            if isfile(bathy_subset_path)
                bathy_subset = Raster(bathy_subset_path, mappedcrs=EPSG(EPSG_code))
            else
                bathy_dataset = Raster(bathy_fullset_path, mappedcrs=EPSG(EPSG_code))
                bathy_subset = extract_subset(bathy_dataset, subset)
                write(bathy_subset_path, bathy_subset; force=true)
            end
            ms_exclusion_zones_bool = create_exclusion_zones(bathy_subset, ms_depth)
            write(ms_exclusion_tif_path, convert.(Int16, ms_exclusion_zones_bool); force=true)
        end

        ms_exclusion_zones_df = ms_exclusion_zones_bool |> to_multipolygon |> to_dataframe
        GDF.write(
            ms_exclusion_gpkg_path,
            ms_exclusion_zones_df;
            crs=EPSG(EPSG_code)
        )
    end



    if isfile(ms_exclusion_gpkg_path)
        return GDF.read(ms_exclusion_gpkg_path)
    else
        if isfile(ms_exclusion_tif_path)
            ms_exclusion_zones_int = Raster(ms_exclusion_tif_path, mappedcrs=EPSG(EPSG_code))
            ms_exclusion_zones_bool = Raster(ms_exclusion_zones_int .== 0, dims(ms_exclusion_zones_int))
        else
            # Load environmental constraints
            # Bathymetry
            if isfile(bathy_subset_path)
                bathy_subset = Raster(bathy_subset_path, mappedcrs=EPSG(EPSG_code))
            else
                bathy_dataset = Raster(bathy_fullset_path, mappedcrs=EPSG(EPSG_code))
                bathy_subset = extract_subset(bathy_dataset, subset)
                write(bathy_subset_path, bathy_subset; force=true)
            end
            ms_exclusion_zones_bool = create_exclusion_zones(bathy_subset, ms_depth)
            write(ms_exclusion_tif_path, convert.(Int16, ms_exclusion_zones_bool); force=true)
        end

        ms_exclusion_zones_df = ms_exclusion_zones_bool |> to_multipolygon |> to_dataframe
        GDF.write(
            ms_exclusion_gpkg_path,
            ms_exclusion_zones_df;
            crs=EPSG(EPSG_code)
        )
        return ms_exclusion_zones_df
    end
end
