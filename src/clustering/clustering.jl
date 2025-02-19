
using Clustering

@kwdef struct Cluster
    id::Int
    centroid::Point{2, Float64}
    nodes::Vector{Point{2, Float64}} = [centroid]
    # TODO: Add waypoints?
    # waypoints::NTuple{2, Point{2, Float64}}
end

"""
    apply_kmeans_clustering(raster::Raster{Int, 2}, k::Int64; tol::Float64=1.0)

Cluster targets sites by applying k-means to non-zero cells in a raster.

# Arguments
- `raster::Raster{Int, 2}` : Raster containing the target geometries.
- `k::Int64` : Number of clusters to create.
- `tol::Float64` : Tolerance for kmeans convergence.

# Returns
A new raster containing the cluster IDs.
"""
function apply_kmeans_clustering(raster::Raster{Int, 2}, k::Int64; tol::Float64=1.0)
    indices = findall(x -> x != 0, raster)
    coordinates = [(Tuple(index)[1], Tuple(index)[2]) for index in indices]
    coordinates_array = hcat([collect(c) for c in coordinates]...)

    clustering = kmeans(coordinates_array, k; tol=tol)

    rows = [coord[1] for coord in coordinates]
    cols = [coord[2] for coord in coordinates]

    clustered_targets = copy(raster)
    clustered_targets .= 0

    for i in 1:length(rows)
        clustered_targets[rows[i], cols[i]] = clustering.assignments[i]
    end

    return clustered_targets
end

"""
    calculate_cluster_centroids(cluster_raster::Raster{Int64, 2})

Calculate the centroids of the clusters in the raster.

# Arguments
- `clusters` : Raster containing the cluster IDs.

# Returns
A vector of `Cluster` objects.
"""
function calculate_cluster_centroids(clusters_raster::Raster{Int64, 2})
    unique_clusters = sort(unique(clusters_raster))
    unique_clusters = unique_clusters[unique_clusters .!= 0]
    clusters = Cluster[]

    x_coords = clusters_raster.dims[1]
    y_coords = clusters_raster.dims[2]
    # Push Cluster object to cluster centroid vector
    for id in unique_clusters
        nodes = [(x_coords[i[1]], y_coords[i[2]]) for i in findall(==(id), clusters_raster)]
        col_cent = mean([node[1] for node in nodes])
        row_cent = mean([node[2] for node in nodes])

        push!(clusters, Cluster(
            id = id,
            centroid = Point{2, Float64}(col_cent, row_cent),
            nodes = [Point{2, Float64}(node[1], node[2]) for node in nodes])
        )
    end

    return clusters
end

"""
    generate_target_clusters(
        clustered_targets_path::String,
        k::Int,
        cluster_tolerance::Float64,
        suitable_targets_all_path::String,
        suitable_threshold::Float64,
        target_subset_path::String,
        subset::DataFrame,
        EPSG_code::Int
    )

Generate a clustered targets raster by reading in the suitable target data,
applying thresholds and cropping to a target subset area, and then clustering.

# Arguments
- `clustered_targets_path::String` : Path to the clustered targets raster.
- `k::Int` : Number of clusters to create.
- `cluster_tolerance::Float64` : Tolerance for kmeans convergence.
- `suitable_targets_all_path::String` : Path to the suitable targets raster.
- `suitable_threshold::Float64` : Threshold for suitable targets.
- `target_subset_path::String` : Path to the target subset raster.
- `subset::DataFrame` : DataFrame containing the target geometries.
- `EPSG_code::Int` : EPSG code for the target geometries.

# Returns
A vector of Cluster objects.
"""
function generate_target_clusters(
    clustered_targets_path::String,
    k::Int,
    cluster_tolerance::Float64,
    suitable_targets_all_path::String,
    suitable_threshold::Float64,
    target_subset_path::String,
    subset::DataFrame,
    EPSG_code::Int
)
    cluster_raster = process_targets(
        clustered_targets_path,
        k,
        cluster_tolerance,
        suitable_targets_all_path,
        suitable_threshold,
        target_subset_path,
        subset,
        EPSG_code
    )
    clusters = calculate_cluster_centroids(cluster_raster)
    return clusters
end
