
"""
    disturb_clusters(
        remaining_clusters::Vector{Cluster},
        disturbance_df::DataFrame,
        current_location::Point{2, Float64},
        exclusions::DataFrame,
        total_tender_capacity::Int;
        tol::Float64 = 0.01
    )::Vector{Cluster}

Disturb the clusters by randomly removing nodes from them.

# Arguments
- `remaining_clusters`: A vector of `Cluster` objects representing clusters to be disturbed.
- `disturbance_df`: A DataFrame containing disturbance data for each node.
- `current_location`: The current location of the mothership at the time of disturbance.
- `exclusions`: A DataFrame containing exclusion zones that should not be disturbed.
- `total_tender_capacity`: The total capacity available for the tender fleet.
- `tol`: A tolerance value for distance calculations.

# Returns
A vector of `Cluster` objects with nodes removed from them.
"""
function disturb_clusters(
    remaining_clusters::Vector{Cluster},
    disturbance_df::DataFrame,
    current_location::Point{2,Float64},
    exclusions::DataFrame,
    total_tender_capacity::Int;
    tol::Float64=0.01
)::Vector{Cluster}
    num_clusters = length(remaining_clusters)
    remaining_nodes = [node for cluster in remaining_clusters for node in cluster.nodes]

    # Filter the disturbance_df so that it only includes rows with nodes in remaining_nodes
    filtered_df = filter(row -> row.node in remaining_nodes, disturbance_df)

    if nrow(filtered_df) != length(remaining_nodes)
        @warn "Warning: Expected $(length(remaining_nodes)) nodes, but filtered df contains $(nrow(filtered_df))"
    end

    disturbed_targets = disturb_remaining_clusters(
        filtered_df,
        num_clusters,
        current_location,
        exclusions,
        total_tender_capacity;
        tol
    )

    # Update the cluster assignments based on previous numbering
    updated_disturbed_targets = update_cluster_assignments(
        disturbed_targets,
        Dict(c.id => c.centroid for c in remaining_clusters)
    )

    disturbed_clusters = calculate_cluster_centroids(
        updated_disturbed_targets;
    )

    return disturbed_clusters
end
