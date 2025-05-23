
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
    disturbance_df::DataFrame,
    current_location::Point{2, Float64},
    exclusions::DataFrame;
    res::Float64 = 0.0001
)::Vector{Cluster}
    num_clusters = length(remaining_clusters)
    remaining_nodes = [node for cluster in remaining_clusters for node in cluster.nodes]

    # Filter the disturbance_df so that it only includes rows with nodes in remaining_nodes
    filtered_df = filter(row -> row.node in remaining_nodes, disturbance_df)

    if nrow(filtered_df) != length(remaining_nodes)
        @warn "Warning: Expected $(length(remaining_nodes)) nodes, but filtered df contains $(nrow(filtered_df))"
        #? Why mismatch? Nodes in exclusions? Why not previously identified?
    end

    disturbance_raster = Rasters.rasterize(
        last,
        [(t[1],t[2]) for t in filtered_df.node];
        res = res,
        missingval = -9999.0,
        fill = filtered_df.disturbance_value,
    )
    disturbed_targets = apply_kmeans_clustering(
        disturbance_raster,
        Int8(num_clusters),
        current_location,
        exclusions
    )

    if all(x -> x == disturbed_targets.missingval, disturbed_targets)
        return remaining_clusters
    end

    # Update the cluster assignments based on previous numbering
    updated_disturbed_targets = update_cluster_assignments(
        disturbed_targets,
        Dict(c.id => c.centroid for c in remaining_clusters)
    )

    disturbed_clusters = calculate_cluster_centroids(
        updated_disturbed_targets;
        cluster_ids=getfield.(remaining_clusters, :id)
    )

    return disturbed_clusters
end
