using Revise, Infiltrator

includet("MS-T_routing.jl")


crs = 7844
depot = (-16.9, 146.15)  # (-16.92, 145.78) # (lat, long)
k = 10
ms_depth = -10

subset = GDF.read("../data/site/Moore_2024-02-14b_v060_rc1.gpkg")
bathy_dataset = Raster("../data/bathy/Cairns-Cooktown/Cairns-Cooktown_bathy.tif", mappedcrs=EPSG(crs))
suitable_targets_all = Raster("../data/targets/Cairns-Cooktown_suitable_slopes.tif", mappedcrs=EPSG(crs))

# Read/mask/trim targets to deployment scale
# TODO: Separate data loading.

# Read deployment locations
if isfile("../outputs/target_subset.tif")
    suitable_targets_subset = Raster("../outputs/target_subset.tif", mappedcrs=EPSG(crs))
else
    suitable_targets_subset = extract_subset(suitable_targets_all, subset)
    write("../outputs/target_subset.tif", suitable_targets_subset; force=true)
end

# Load environmental constraints
# Bathymetry
if isfile("../outputs/bathy_subset.tif")
    target_bathy = Raster("../outputs/bathy_subset.tif", mappedcrs=EPSG(crs))
else
    target_bathy = extract_subset(bathy_dataset, subset)
    write("../outputs/bathy_subset.tif", target_bathy; force=true)
end

# Cluster targets
if isfile("../outputs/clustered_targets.tif")
    clustered_targets = Raster("../outputs/clustered_targets$(k).tif", mappedcrs=EPSG(crs))
else
    clustered_targets = cluster_targets(suitable_targets_subset, k)
    write("../outputs/clustered_targets_k=$(k).tif", clustered_targets; force=true)
end

# Create exclusion zones from environmental constraints
if isfile("../outputs/ms_exclusion.tif")
    ms_exclusion_zones_int = Raster("../outputs/ms_exclusion.tif", mappedcrs=EPSG(crs))
    ms_exclusion_zones = Raster(ms_exclusion_zones_int .!= 0, dims(ms_exclusion_zones_int))
else
    ms_exclusion_zones = create_exclusion_zones(target_bathy, ms_depth)
    write("../outputs/ms_exclusion.tif", convert.(Int16, ms_exclusion_zones); force=true)
end

########


mp = to_multipolygon(ms_exclusion_zones)
exclusion_zones = to_dataframe(mp)

# ms_exclusion_zones |> to_multipolygon |> to_dataframe |> x -> poly(x.geometry)

# Generate initial mothership route
cluster_centroids = calc_cluster_centroids(clustered_targets, depot)

cluster_sequence, mothership_dist, clust_matrix = nearest_neighbour(cluster_centroids)

waypoints = calc_waypoints(cluster_centroids, cluster_sequence)

plot_mothership_route(clustered_targets, cluster_centroids, cluster_sequence)

# # Apply 2-opt to improve the mothership route
# cluster_sequence, mothership_waypoints, mothership_dist = two_opt(cluster_centroids, depot, cluster_sequence)
