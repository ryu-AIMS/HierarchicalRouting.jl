
"""
    extract_subset(spatial_dataset::Raster, subset)

Extract a subset from the given `spatial_dataset` raster based on the specified `subset`.

# Arguments
- `spatial_dataset::Raster`: The spatial dataset from which to extract the subset.
- `subset`: The subset criteria to apply.

# Returns
A new raster containing the subset of the original raster.
"""
function extract_subset(spatial_dataset::Raster, subset)
    return Rasters.trim(mask(spatial_dataset; with=subset.geom))
end

"""
    cluster_targets(df::DataFrame, num_clust::Int64)

Cluster the targets in a GeoDataFrame based on their geometry.

# Arguments
- `raster::Raster{Int16, 2}` : Raster containing the target geometries.
- `num_clust::Int64` : Number of clusters to create.

# Returns
A DataFrame containing the cluster ID and the target geometries.
"""
function cluster_targets(raster::Raster{Int16, 2}, num_clust::Int64)
    indices = findall(x -> x != 0, raster)
    coordinates = [(Tuple(index)[1], Tuple(index)[2]) for index in indices]
    coordinates_array = hcat([collect(c) for c in coordinates]...)

    clustering = kmeans(coordinates_array, num_clust)

    rows = [coord[1] for coord in coordinates]
    cols = [coord[2] for coord in coordinates]

    cluster_raster = copy(raster)
    cluster_raster .= 0

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

    coordinates = [(x, y) for x in clusters.dims[1], y in clusters.dims[2]]

    cluster_centroids = DataFrame(cluster_id = Int[], lat = Float64[], lon = Float64[])
    # Add depot as cluster 0
    push!(cluster_centroids, (0, depot[1], depot[2]))

    # Add each cluster centroid to the DataFrame
    for cluster_id in unique_clusters
        indices = findall(x -> x == cluster_id, clusters)

        lat = mean([coordinates[i][2] for i in indices])
        lon = mean([coordinates[i][1] for i in indices])

        push!(cluster_centroids, (cluster_id, lat, lon))
    end

    # Sort the cluster centroids by cluster_id
    sort!(cluster_centroids, :cluster_id)
    return cluster_centroids
end
