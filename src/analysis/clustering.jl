
"""
    extract_subset(spatial_dataset::Raster, subset::GeoDataFrame)

Extract a subset of a raster dataset based on a GeoDataFrame.
"""
function extract_subset(spatial_dataset::Raster, subset)
    return Rasters.trim(mask(spatial_dataset; with=subset.geom))
end

"""
    cluster_targets(df::DataFrame, num_clust::Int64)

Cluster the targets in a GeoDataFrame based on their geometry.
"""
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

    return cluster_raster
end

"""
    centroids(cluster_raster::Raster{Int16, 2}, depot::Tuple{Float64, Float64})

Calculate the centroids of the clusters in the raster.
Depot included as cluster 0.

# Arguments
- `clusters` : Raster containing the cluster IDs.
- `depot` : Coordinates of the depot.

# Returns
A DataFrame containing the cluster ID, latitude, and longitude of the centroids.
"""
function centroids(clusters::Raster{Int16, 2}, depot::Tuple{Float64, Float64})
    unique_clusters = unique(clusters)

    # Remove the zero cluster ID (if present)
    unique_clusters = unique_clusters[unique_clusters .!= 0]

    cluster_centroids = DataFrame(cluster_id = Int[], lat = Float64[], lon = Float64[])

    # Get the latitude and longitude coordinates
    X_dim = clusters.dims[1]
    Y_dim = clusters.dims[2]

    # Extracting coordinates using DimensionalData
    coordinates = [(x, y) for x in X_dim, y in Y_dim]

    # Add depot as cluster 0
    push!(cluster_centroids, (0, depot[1], depot[2]))
    for cluster_id in unique_clusters
        # Find the indices of the current cluster_id
        indices = findall(x -> x == cluster_id, clusters)

        # Calculate the mean coordinates for the cluster
        lat = mean([coordinates[i][2] for i in indices])
        lon = mean([coordinates[i][1] for i in indices])

        # Add the centroid to the DataFrame
        push!(cluster_centroids, (cluster_id, lat, lon))
    end

    # Sort the cluster centroids by cluster_id
    sort!(cluster_centroids, :cluster_id)
    return cluster_centroids
end
