
"""
    process_targets(
    clustered_targets_path::String,
    k::Int,
    cluster_tolerance::Real,
    suitable_targets_all_path::String,
    suitable_threshold::Real,
    target_subset_path::String,
    subset::DataFrame,
    EPSG_code::Int
    )

Process target locations to generate clusters.

# Arguments
- `clustered_targets_path::String`: The path to the clustered targets raster.
- `k::Int`: The number of clusters.
- `cluster_tolerance::Real`: The cluster tolerance.
- `suitable_targets_all_path::String`: The path to the suitable targets dataset.
- `suitable_threshold::Real`: The suitable targets threshold.
- `target_subset_path::String`: The path to the target subset raster.
- `subset::DataFrame`: The DataFrame containing the study area boundary.
- `EPSG_code::Int`: The EPSG code for the study area.

# Returns
- `clustered_targets::Raster`: The clustered targets raster, classified by cluster ID number.
"""
function process_targets(
    clustered_targets_path,
    k,
    cluster_tolerance,
    suitable_targets_all_path,
    suitable_threshold,
    target_subset_path,
    subset,
    EPSG_code
)
    # TODO Generalize
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
        clustered_targets = cluster_raster(suitable_targets_subset, k; tol=cluster_tolerance)
        write(clustered_targets_path, clustered_targets; force=true)
        return clustered_targets
    end
end

"""
    process_exclusions(
    bathy_fullset_path,
    vessel_draft,
    subset,
    EPSG_code,
    bathy_subset_path,
    exclusion_gpkg_path,
    exclusion_tif_path
    )

Create exclusion zones from environmental constraints.

# Arguments
- `bathy_fullset_path::String`: The path to the full bathymetry dataset.
- `vessel_draft::Real`: The vessel draft/depth.
- `subset::DataFrame`: The DataFrame containing the study area boundary.
- `EPSG_code::Int`: The EPSG code for the study area.
- `bathy_subset_path::String`: The path to the subset bathymetry dataset.
- `exclusion_gpkg_path::String`: The path to the exclusion zones GeoPackage.
- `exclusion_tif_path::String`: The path to the exclusion zones raster.

# Returns
- `exclusion_zones_df::DataFrame`: The DataFrame containing the exclusion zones.
"""
function process_exclusions(
    bathy_fullset_path,
    vessel_draft,
    subset,
    EPSG_code,
    bathy_subset_path,
    exclusion_gpkg_path,
    exclusion_tif_path
)
    # TODO: Generalize for all available environmental constraints
    # TODO: Generalize for ms and tender vessels
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

"""
    process_problem(problem::Problem)

Read and process problem data to generate an initial solution.

# Arguments
- `problem::Problem`: The problem data in the form of a `HierarchicalRouting::Problem` struct.

# Returns
- `clusters::Vector{Cluster}`: The clusters of the problem.
- `cluster_centroids_df::DataFrame`: The DataFrame containing the cluster centroids.
"""
function process_problem(problem::Problem)
    config = TOML.parsefile(joinpath("src",".config.toml"))

    suitable_threshold = config["parameters"]["suitable_threshold"]
    k = config["parameters"]["k"]
    cluster_tolerance = config["parameters"]["cluster_tolerance"]
    EPSG_code = config["parameters"]["EPSG_code"]

    target_scenario = problem.target_scenario
    suitable_targets_prefix = target_scenario[1:findlast(".",target_scenario)[1]-1]

    site_dir = config["data_dir"]["site"]
    output_dir = config["output_dir"]["path"]

    clustered_targets_path = joinpath(output_dir, "clustered_$(suitable_targets_prefix)_targets_k=$(k).tif")
    target_subset_threshold_path = joinpath(output_dir, "target_subset_$(suitable_targets_prefix)_threshold=$(suitable_threshold).tif")
    target_subset_path = joinpath(output_dir, "target_subset_$(suitable_targets_prefix).tif")

    subset = GDF.read(first(glob("*.gpkg", site_dir)))

    target_scenario_dir = config["data_dir"]["target_scenarios"]
    suitable_targets_all_path = joinpath(target_scenario_dir, target_scenario)

    clusters = cluster_targets(
        clustered_targets_path,
        k,
        cluster_tolerance,
        suitable_targets_all_path,
        suitable_threshold,
        target_subset_path,
        subset,
        EPSG_code
    )

    cluster_centroids_df = DataFrame(id = Int[0], lon = Float64[problem.depot[1]], lat = Float64[problem.depot[2]])
    [push!(cluster_centroids_df, (i, clust.centroid[1], clust.centroid[2])) for (i, clust) in enumerate(clusters)]

    return clusters, cluster_centroids_df
end
