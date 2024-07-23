import ArchGDAL as AG
using Rasters
using DataFrames
using Statistics
using Distances
using Clustering

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
function cluster_targets(raster::Raster{Int16, 2}, num_clust::Int64)
    # Extract the coordinates of non-zero values
    indices = findall(x -> x != 0, raster)

    # Convert the indices to coordinates (tuples)
    coordinates = [(Tuple(index)[1], Tuple(index)[2]) for index in indices]
    # Convert the coordinates to a format suitable for clustering (e.g., an array)
    coordinates_array = hcat([collect(c) for c in coordinates]...)

    # Perform k-means clustering on the coordinates
    clustering = kmeans(coordinates_array, num_clust)

    # Create a DataFrame to store the cluster assignments
    rows = [coord[1] for coord in coordinates]
    cols = [coord[2] for coord in coordinates]
    # cluster_df = DataFrame(row = rows, col = cols, cluster_id = clustering.assignments)

    # Create a new raster to store the cluster IDs
    cluster_raster = copy(raster)
    cluster_raster .= 0  # Initialize with zeros

    # Assign the cluster IDs to the corresponding positions in the new raster
    for i in 1:length(rows)
        cluster_raster[rows[i], cols[i]] = clustering.assignments[i]
    end

    return cluster_raster#, cluster_df
end
