import GeoDataFrames as GDF
using Rasters
using GLMakie

crs = 7844
depot = (-16.9, 146.15)#(-16.92, 145.78) # (lat, long)
k = 3

subset = GDF.read("../data/site/Moore_2024-02-14b_v060_rc1.gpkg")
bathy_dataset = Raster("../data/bathy/Cairns-Cooktown/Cairns-Cooktown_bathy.tif", mappedcrs=EPSG(crs))
suitable_targets_all = Raster("../data/targets/Cairns-Cooktown_suitable_slopes.tif", mappedcrs=EPSG(crs))

# Read/mask/trim targets to deployment scale

# Read deployment locations
suitable_targets_subset = extract_subset(suitable_targets_all, subset)
plot(suitable_targets_subset)
write("../outputs/target_subset.tif", suitable_targets_subset; force=true)

# Environmental constraints
# Bathymetry
target_bathy = extract_subset(bathy_dataset, subset)
plot(target_bathy)
write("../outputs/bathy_subset.tif", target_bathy; force=true)

# Cluster targets
clustered_targets = cluster_targets(suitable_targets_subset, k)
write("../outputs/clustered_targets.tif", clustered_targets; force=true)
plot(clustered_targets)

cluster_centroids = calc_cluster_centroids(clustered_targets, depot)

# Generate initial mothership route
cluster_sequence, mothership_dist = nearest_neighbour(cluster_centroids)

# Plot mothership route
plot_mothership_route(clustered_targets, cluster_centroids, cluster_sequence)

# calc waypoints
waypoints = calc_waypoints(cluster_centroids, cluster_sequence)
