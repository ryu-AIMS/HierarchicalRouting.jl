import ArchGDAL as AG
using DataFrames
using Statistics

"""
    extract_subset(spatial_dataset::Raster, subset::GeoDataFrame)

Extract a subset of a raster dataset based on a GeoDataFrame.
"""
function extract_subset(spatial_dataset::Raster, subset)
    result_raster = Rasters.trim(mask(spatial_dataset; with=subset.geom))
    return result_raster
end

"""
    cluster_targets(df::DataFrame, num_clust::Int64)

Cluster the targets in a GeoDataFrame based on their geometry.
"""
function cluster_targets(df::DataFrame, num_clust::Int64)
    # Calculate centroid of geometry for each row
    centroid_shp = [AG.centroid(row.geom) for row in eachrow(df)]
    centroid_coords = [(AG.getx(centroid,0), AG.gety(centroid,0)) for centroid in centroid_shp]

    # Convert the coordinates to a format suitable for clustering (e.g., an array)
    coordinates_array = hcat([collect(c) for c in centroid_coords]...)

    # Cluster centroids using kmeans
    clustering = kmeans(coordinates_array, num_clust)

    df.cluster_id = clustering.assignments
    return df
end
