import GeoDataFrames as GDF
using Rasters
using GLMakie

# Read deployment locations
suitable_targets = GDF.read("../data/target/Moore_2024-02-14b_v060_rc1.gpkg")

# Environmental constraint data: Read/mask/trim to suitable targets
# Bathymetry
bathy_dataset = Raster("../data/bathy/Cairns-Cooktown/Cairns-Cooktown_bathy.tif")
target_bathy = extract_subset(bathy_dataset, suitable_targets)
plot(target_bathy)
write("../outputs/bathy_subset.tif", target_bathy; force=true)
