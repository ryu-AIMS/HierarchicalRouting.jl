
using Clustering

@kwdef struct Cluster
    id::Int
    centroid::Point{2, Float64}
    nodes::Vector{Point{2, Float64}} = [centroid]
    # TODO: Add waypoints?
    # waypoints::NTuple{2, Point{2, Float64}}
end

"""
    cluster_targets(raster::Raster{Int, 2}, k::Int64)

Cluster targets using kmeans clustering.

# Arguments
- `raster::Raster{Int, 2}` : Raster containing the target geometries.
- `k::Int64` : Number of clusters to create.
- `tol::Float64` : Tolerance for kmeans convergence.

# Returns
A raster containing the cluster IDs.
"""
function cluster_raster(raster::Raster{Int, 2}, k::Int64; tol::Float64=1.0)
    indices = findall(x -> x != 0, raster)
    coordinates = [(Tuple(index)[1], Tuple(index)[2]) for index in indices]
    coordinates_array = hcat([collect(c) for c in coordinates]...)

    clustering = kmeans(coordinates_array, k; tol=tol)

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
    create_clusters(clusters::Raster{Int64, 2}, depot=nothing)

Calculate the centroids of the clusters in the raster.
Depot included as cluster 0.

# Arguments
- `clusters` : Raster containing the cluster IDs.
- `depot` : Optional depot point.

# Returns
A vector of Cluster objects.
"""
function create_clusters(clusters::Raster{Int64, 2}, depot=nothing)
    unique_clusters = sort(unique(clusters))

    # Remove the zero cluster ID (if present)
    unique_clusters = unique_clusters[unique_clusters .!= 0]

    cluster_vec = Cluster[]

    if depot !== nothing
        push!(cluster_vec, Cluster(id = 0, centroid = depot))
    end

    x_coords = clusters.dims[1]
    y_coords = clusters.dims[2]
    # Push Cluster object to cluster centroid vector
    for id in unique_clusters
        nodes = [(x_coords[i[1]], y_coords[i[2]]) for i in findall(==(id), clusters)]

        col_cent = mean([node[1] for node in nodes])
        row_cent = mean([node[2] for node in nodes])

        push!(cluster_vec, Cluster(id = id, centroid = Point{2, Float64}(col_cent, row_cent), nodes = [Point{2, Float64}(node[1], node[2]) for node in nodes]))
    end

    return cluster_vec
end

"""
    cluster_targets(
        clustered_targets_path::String,
        k::Int,
        cluster_tolerance::Float64,
        suitable_targets_all_path::String,
        suitable_threshold::Float64,
        target_subset_path::String,
        subset::DataFrame,
        EPSG_code::Int
    )

Cluster targets based on their geometry.

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
function cluster_targets(
    clustered_targets_path::String,
    k::Int,
    cluster_tolerance::Float64,
    suitable_targets_all_path::String,
    suitable_threshold::Float64,
    target_subset_path::String,
    subset::DataFrame,
    EPSG_code::Int
)

    clusters_raster = process_targets(
        clustered_targets_path,
        k,
        cluster_tolerance,
        suitable_targets_all_path,
        suitable_threshold,
        target_subset_path,
        subset,
        EPSG_code
    )
    clusts = create_clusters(clusters_raster)
    return clusts
end
