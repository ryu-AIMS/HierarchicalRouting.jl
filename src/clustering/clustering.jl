
using Clustering
using Random

@kwdef struct Cluster
    id::Int
    centroid::Point{2, Float64}
    nodes::Vector{Point{2, Float64}} = [centroid]
    # TODO: Add waypoints?
    # waypoints::NTuple{2, Point{2, Float64}}
end

"""
    apply_kmeans_clustering(
        raster::Raster{Int, 2}, k::Int8; tol::Float64=1.0
    )::Raster{Int64, 2}

Cluster targets sites by applying k-means to target (non-zero) cells in a raster.

# Arguments
- `raster`: Raster containing the target geometries.
- `k`: Number of clusters to create.
- `tol`: Tolerance for kmeans convergence.

# Returns
A new raster containing the cluster IDs.
"""
function apply_kmeans_clustering(
    raster::Raster{Int, 2}, k::Int8; tol::Float64=1.0
)::Raster{Int64, 2}
    indices::Vector{CartesianIndex{2}} = findall(!=(0), raster)
    n::Int = length(indices)
    coordinates_array = Matrix{Float64}(undef, 2, n)

    # Row/col indices from each CartesianIndex
    rows::Vector{Int64} = getindex.(indices, 1)
    cols::Vector{Int64} = getindex.(indices, 2)

    # Fill the coordinate matrix using the corresponding dimension arrays
    coordinates_array[1, :] .= raster.dims[1][rows]
    coordinates_array[2, :] .= raster.dims[2][cols]

    clustering = kmeans(coordinates_array, k; distance=Haversine(), tol=tol, rng=Random.seed!(1))

    clustered_targets::Raster{Int64, 2} = similar(raster, Int64)
    clustered_targets[indices] .= clustering.assignments

    return clustered_targets
end

"""
    calculate_cluster_centroids(cluster_raster::Raster{Int64, 2})::Vector{Cluster}

Calculate the centroids of the clusters in the raster.

# Arguments
- `clusters_raster`: Raster containing the cluster IDs.

# Returns
A vector of `Cluster` objects.
"""
function calculate_cluster_centroids(clusters_raster::Raster{Int64, 2})::Vector{Cluster}
    unique_clusters = sort(unique(clusters_raster[clusters_raster .!= 0]))
    clusters_vector = Vector{Cluster}(undef, length(unique_clusters))

    x_coords = clusters_raster.dims[1]
    y_coords = clusters_raster.dims[2]

    for id in unique_clusters
        nodes = [(x_coords[i[1]], y_coords[i[2]]) for i in findall(==(id), clusters_raster)]
        col_cent = mean([node[1] for node in nodes])
        row_cent = mean([node[2] for node in nodes])

        clusters_vector[id] = Cluster(
            id = id,
            centroid = Point{2, Float64}(col_cent, row_cent),
            nodes = [Point{2, Float64}(node[1], node[2]) for node in nodes]
        )
    end
    return clusters_vector
end

"""
    generate_target_clusters(
        clustered_targets_path::String,
        k::Int8,
        cluster_tolerance::Float64,
        suitable_targets_all_path::String,
        suitable_threshold::Float64,
        target_subset_path::String,
        subset::DataFrame,
        EPSG_code::Int16
    )::Vector{Cluster}

Generate a clustered targets raster by reading in the suitable target data,
applying thresholds and cropping to a target subset area, and then clustering.

# Arguments
- `clustered_targets_path`: Path to the clustered targets raster.
- `k`: Number of clusters to create.
- `cluster_tolerance`: Tolerance for kmeans convergence.
- `suitable_targets_all_path`: Path to the suitable targets raster.
- `suitable_threshold`: Threshold for suitable targets.
- `target_subset_path`: Path to the target subset raster.
- `subset`: DataFrame containing the target geometries.
- `EPSG_code`: EPSG code for the target geometries.

# Returns
A vector of Cluster objects.
"""
function generate_target_clusters(
    clustered_targets_path::String,
    k::Int8,
    cluster_tolerance::Float64,
    suitable_targets_all_path::String,
    suitable_threshold::Float64,
    target_subset_path::String,
    subset::DataFrame,
    EPSG_code::Int16
)::Vector{Cluster}
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
    return calculate_cluster_centroids(cluster_raster)
end
