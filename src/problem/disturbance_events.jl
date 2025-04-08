
"""
    disturb_clusters(
    remaining_clusters::Vector{Cluster},
    disturbance_df::DataFrame;
    res::Float64 = 0.0001
    )

Disturb the clusters by randomly removing nodes from them.

# Arguments
- `remaining_clusters`: A vector of `Cluster` objects representing clusters to be disturbed.
- `disturbance_df`: A DataFrame containing disturbance data for each node.
- `res`: The resolution of the raster to be created from the disturbance data.

# Returns
A vector of `Cluster` objects with nodes removed from them.
"""
function disturb_clusters(
    remaining_clusters::Vector{Cluster},
    disturbance_df::DataFrame;
    res::Float64 = 0.0001
)::Vector{Cluster}

    num_clusters = length(remaining_clusters)
    remaining_nodes = [node for cluster in remaining_clusters for node in cluster.nodes]

    # mask df to only include nodes in the remaining clusters
    filtered_df = DataFrame()
    for cluster in remaining_clusters
        append!(filtered_df, filter(row -> row.node in cluster.nodes, disturbance_df))
    end

    disturbance_raster = Rasters.rasterize(
        last,
        [(t[1],t[2]) for t in filtered_df.node];
        res = res,
        missingval = -9999.0,
        fill = filtered_df.wave_value,
    )
    reclustered_disturbed_targets = apply_kmeans_clustering(
        disturbance_raster,
        Int8(num_clusters);
    )
    if all(x -> x == reclustered_disturbed_targets.missingval, reclustered_disturbed_targets)
        return remaining_clusters
    end
    disturbed_clusters = calculate_cluster_centroids(
        reclustered_disturbed_targets;
        cluster_ids=getfield.(remaining_clusters, :id)
    )

    return disturbed_clusters
end
