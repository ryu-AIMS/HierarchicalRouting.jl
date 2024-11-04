
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
