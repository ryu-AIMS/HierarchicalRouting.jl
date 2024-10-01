
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
    centroids(clusters::Raster{Int16, 2}, depot::Tuple{Float64, Float64})

Calculate the centroids of the clusters in the raster.
Depot included as cluster 0.

# Arguments
- `clusters` : Raster containing the cluster IDs.
- `depot` : Coordinates of the depot.

# Returns
A DataFrame containing the ID, latitude, and longitude of the centroids.
"""
function centroids(clusters::Raster{Int16, 2}, depot::Tuple{Float64, Float64})
    unique_clusters = unique(clusters)

    # Remove the zero cluster ID (if present)
    unique_clusters = unique_clusters[unique_clusters .!= 0]

    coordinates = [(x, y) for x in clusters.dims[1], y in clusters.dims[2]]

    cluster_centroids = DataFrame(id = Int[], lat = Float64[], lon = Float64[])
    # Add depot as cluster 0
    push!(cluster_centroids, (0, depot[1], depot[2]))

    # Add each cluster centroid to the DataFrame
    for id in unique_clusters
        indices = findall(x -> x == id, clusters)

        lat = mean([coordinates[i][2] for i in indices])
        lon = mean([coordinates[i][1] for i in indices])

        push!(cluster_centroids, (id, lat, lon))
    end

    # Sort the centroids by id
    sort!(cluster_centroids, :id)
    return cluster_centroids
end
function centroids(clusters::Raster{Int16, 2})
    unique_clusters = unique(clusters)

    # Remove the zero cluster ID (if present)
    unique_clusters = unique_clusters[unique_clusters .!= 0]

    coordinates = [(x, y) for x in clusters.dims[1], y in clusters.dims[2]]

    cluster_centroids = DataFrame(id = Int[], lat = Float64[], lon = Float64[])

    # Add each cluster centroid to the DataFrame
    for id in unique_clusters
        indices = findall(x -> x == id, clusters)

        lat = mean([coordinates[i][2] for i in indices])
        lon = mean([coordinates[i][1] for i in indices])

        push!(cluster_centroids, (id, lat, lon))
    end

    # Sort the centroids by id
    sort!(cluster_centroids, :id)
    return cluster_centroids
end
