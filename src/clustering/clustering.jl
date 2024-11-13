
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

Cluster the targets in a GeoDataFrame based on their geometry.

# Arguments
- `raster::Raster{Int, 2}` : Raster containing the target geometries.
- `k::Int64` : Number of clusters to create.
- `tol::Float64` : Tolerance for kmeans convergence.

# Returns
A DataFrame containing the cluster ID and the target geometries.
"""
function cluster_raster(raster::Raster{Int, 2}, k::Int64; tol::Float64=1.0)
    indices = findall(x -> x != 0, raster)
    coordinates = [(Tuple(index)[1], Tuple(index)[2]) for index in indices]
    coordinates_array = hcat([collect(c) for c in coordinates]...)

    clustering = kmeans(coordinates_array, k; tol=1)

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

    # Push Cluster object to cluster centroid vector
    for id in unique_clusters
        nodes = [(i[1], i[2]) for i in findall(==(id), clusters)]

        row_cent = mean([node[1] for node in nodes])
        col_cent = mean([node[2] for node in nodes])

        push!(cluster_vec, Cluster(id = id, centroid = Point{2, Float64}(col_cent, row_cent), nodes = [Point{2, Float64}(node[2], node[1]) for node in nodes]))
    end

    return cluster_vec
end

function cluster_targets(
    clustered_targets_path::String,
    target_subset_threshold_path::String,
    k::Int,
    cluster_tolerance::Float64,
    suitable_targets_all_path::String,
    suitable_threshold::Float64,
    target_subset_path::String,
    subset::DataFrame,
    EPSG_code::Int)

    clusters_raster = process_targets(
        clustered_targets_path,
        target_subset_threshold_path,
        k, cluster_tolerance,
        suitable_targets_all_path,
        suitable_threshold,
        target_subset_path,
        subset,
        EPSG_code
    )
    clusts = create_clusters(clusters_raster)
    return clusts
end
