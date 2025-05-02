
using Clustering
using Random
using Hungarian

@kwdef struct Cluster
    id::Int
    centroid::Point{2, Float64}
    nodes::Vector{Point{2, Float64}} = [centroid]
    # TODO: Add waypoints?
    # waypoints::NTuple{2, Point{2, Float64}}
end

"""
    generate_cluster_df(
        clusters::Vector{Cluster},
        depot::Point{2, Float64}
    )::DataFrame

Generate a DataFrame representing clusters vector containing the cluster centroids and their
IDs.

# Arguments
- `clusters`: A vector of `Cluster` objects.
- `depot`: A point representing the start/end of the route.

# Returns
- A DataFrame with the clusters: ID's and lon/lat of centroids.
"""
function generate_cluster_df(
    clusters::Vector{Cluster},
    depot::Point{2, Float64}
)::DataFrame
    cluster_centroids_df::DataFrame = DataFrame(
        id  = [0; 1:length(clusters)],
        lon = [depot[1]; [clust.centroid[1] for clust in clusters]],
        lat = [depot[2]; [clust.centroid[2] for clust in clusters]]
    )
    return cluster_centroids_df
end

"""
    apply_kmeans_clustering(
        raster::Raster{Int, 2},
        k::Int8,
        current_location::Point{2, Float64},
        exclusions::DataFrame;
        tol::Float64=1.0,
        dist_weighting::Float64=2E-5
    )::Raster{Int64, 2}
    apply_kmeans_clustering(
        raster::Raster{Float64, 2},
        k::Int8,
        current_location::Point{2, Float64},
        exclusions::DataFrame;
        tol::Float64=1.0,
        dist_weighting::Float64=2E-5
    )::Raster{Int64, 2}

Cluster targets sites by applying k-means to target (non-zero) cells in a raster.
- Float64 raster is assumed to contain disturbance values, addressed buy 3d clustering
- Int raster is assumed to contain cluster ID values, addressed by 2d spatial clustering.
Clustering considers feasible distances from the current location as a 3rd dimenseion, and
    excludes points in exclusion zones.

# Arguments
- `raster`: Raster containing the target geometries.
- `k`: Number of clusters to create.
- `current_location`: Current location of the mothership.
- `exclusions`: DataFrame containing the exclusion zones.
- `tol`: Tolerance for kmeans convergence.
- `dist_weighting`: Weighting factor for the distances (in kms) to be stored in 3d array, compared
    against lat/lons.

# Returns
A raster containing new clusters, with cluster IDs assigned to each target site.
"""
function apply_kmeans_clustering(
    raster::Raster{Int, 2},
    k::Int8,
    current_location::Point{2, Float64},
    exclusions::DataFrame;
    tol::Float64=1.0,
    dist_weighting::Float64=2E-5
)::Raster{Int64, 2}
    indices::Vector{CartesianIndex{2}} = findall(x -> x != raster.missingval, raster)
    rows::Vector{Int64} = getindex.(indices, 1)
    cols::Vector{Int64} = getindex.(indices, 2)

    points = Point{2,Float64}.(zip(raster.dims[1][rows], raster.dims[2][cols]))

    dist_vector = dist_weighting .* get_feasible_distances(
        current_location,
        points,
        exclusions
    )

    feasible_idxs = findall(x -> x != Inf, dist_vector)
    feasible_points = points[feasible_idxs]

    # 3D coordinate matrix of feasible points for clustering
    coordinates_array = Matrix{Float64}(undef, 3, length(feasible_points))
    coordinates_array[1, :] .= getindex.(feasible_points, 1)
    coordinates_array[2, :] .= getindex.(feasible_points, 2)
    coordinates_array[3, :] = dist_vector[feasible_idxs]'

    clustering = kmeans(coordinates_array, k; tol=tol, rng=Random.seed!(1))

    clustered_targets::Raster{Int64, 2} = similar(raster, Int64)
    clustered_targets .= clustered_targets.missingval
    clustered_targets[indices[feasible_idxs]] .= clustering.assignments

    return clustered_targets
end
function apply_kmeans_clustering(
    raster::Raster{Float64, 2},
    k::Int8,
    current_location::Point{2, Float64},
    exclusions::DataFrame;
    tol::Float64=1.0,
    dist_weighting::Float64=2E-5
)::Raster{Int64, 2}
    # TODO: Split this function into two separate functions, it does more than original fn
    #! 1st: 3D clustering to generate disturbance clusters,
    #! 2nd: 3D clustering to generate target clusters (as above)
    indices::Vector{CartesianIndex{2}} = findall(!=(raster.missingval), raster)
    n::Int = length(indices) # number of target sites remaining

    if n <= k
        @warn "No disturbance, as (deployment targets <= clusters required)"
        empty_raster = similar(raster, Int64, missingval=0)
        empty_raster .= empty_raster.missingval
        return empty_raster
    end

    # 3D coordinate matrix for clustering
    coordinates_array_3d = Matrix{Float64}(undef, 3, n)
    coordinates_array_3d[1, :] .= raster.dims[1][getindex.(indices, 1)]
    coordinates_array_3d[2, :] .= raster.dims[2][getindex.(indices, 2)]
    coordinates_array_3d[3, :] .= raster[indices]

    # Create k_d clusters to create disturbance on subset
    k_d_lower = min(n, k+1)
    k_d_upper = min(max(k+1, n, k^2), n)
    k_d = rand(k_d_lower:k_d_upper)

    disturbance_clusters = kmeans(
        coordinates_array_3d,
        k_d;
        tol=tol,
        rng=Random.seed!(1)
    )

    # Create a score based on the disturbance values for each cluster
    disturbance_scores = Vector{Float64}(undef, n)
    # Calculate the mean disturbance value for each cluster with stochastic perturbation
    w = 1.0 # weight for the environmental disturbance value
    t = 1.0 # perturbation weighting factor
    cluster_disturbance_vals = w*[
        mean(coordinates_array_3d[3, disturbance_clusters.assignments .== i])
        for i in 1:k_d
    ] .+ t*rand(-1.0:0.01:1.0, k_d)
    # Assign the disturbance value to every node in the cluster
    disturbance_scores .= cluster_disturbance_vals[disturbance_clusters.assignments]

    # remove nodes in the cluster with the highest disturbance score
    max_disturbance_score = maximum(disturbance_scores)
    surviving_mask = disturbance_scores .!= max_disturbance_score

    coordinates_array_2d_disturbed = coordinates_array_3d[1:2, surviving_mask]
    indices = indices[surviving_mask]
    n = length(indices) # update number of target sites remaining

    if k > n
        #! Too many nodes/clusters removed! Change threshold,
        #! or use a different method e.g. remove cluster with highest scores
        error(
            "$k clusters required from $(n) remaining node/s.\nToo many nodes removed!"
        )
    end

    remaining_pts = Point{2,Float64}.(
        coordinates_array_2d_disturbed[1, :],
        coordinates_array_2d_disturbed[2, :]
    )

    # Fill vector with feasible distance from depot to each target site
    dist_vector = get_feasible_distances(
        current_location,
        remaining_pts,
        exclusions
    )

    # Filter out infeasible points using infeasible_point_indxs
    feasible_idxs = findall(x -> x != Inf, dist_vector)
    coordinates_array_2d_disturbed = coordinates_array_2d_disturbed[:, feasible_idxs]
    filtered_dists = dist_weighting .* dist_vector[feasible_idxs]

    coordinates_array_3d_disturbed = [coordinates_array_2d_disturbed; filtered_dists']

    #re-cluster the remaining nodes into k clusters
    clustering = kmeans(
        coordinates_array_3d_disturbed,
        k;
        tol=tol,
        rng=Random.seed!(2)
    )

    clustered_targets = similar(raster, Int64, missingval=0)
    clustered_targets .= clustered_targets.missingval
    clustered_targets[indices[feasible_idxs]] .= clustering.assignments

    return clustered_targets
end

"""
    update_cluster_assignments(
        cluster_raster::Raster{Int64, 2},
        prev_centroids::Dict{Int64, Point{2,Float64}}
    )::Raster{Int64, 2}

Update cluster assignments in the raster to match previous cluster numbering, based on
matching closest clster centroids.

# Arguments
- `cluster_raster`: Raster containing the new cluster IDs.
- `prev_centroids`: Dictionary mapping previous cluster IDs to their centroids.

# Returns
- A new raster with updated cluster assignments to match the previous numbering.
"""
function update_cluster_assignments(
    cluster_raster::Raster{Int64, 2},
    prev_centroids::Dict{Int64, Point{2,Float64}}
)::Raster{Int64, 2}
    unique_new = Set(cluster_raster[cluster_raster .!= cluster_raster.missingval])

    # Compute new centroids from the cluster_raster
    new_centroids = Dict{Int64, Point{2,Float64}}()
    for new_id in unique_new
        indices = findall(x -> x == new_id, cluster_raster)
        mean_lon = mean(cluster_raster.dims[1][getindex.(indices, 1)])
        mean_lat = mean(cluster_raster.dims[2][getindex.(indices, 2)])
        new_centroids[new_id] = (mean_lon, mean_lat)
    end

    # Map cluster IDs from new to previous clusters
    cluster_mapping = one_to_one_mapping_hungarian(new_centroids, prev_centroids)

    #? default args for similar?
    updated_raster = similar(cluster_raster, Int64, missingval=cluster_raster.missingval)
    updated_raster .= updated_raster.missingval

    for new_id in unique_new
        updated_raster[cluster_raster .== new_id] .= cluster_mapping[new_id]
    end

    return updated_raster
end

function one_to_one_mapping_hungarian(
    new_centroids::Dict{Int64,Point{2,Float64}},
    prev_centroids::Dict{Int64,Point{2,Float64}}
)::Dict{Int64,Int64}
    new_ids  = collect(keys(new_centroids))
    prev_ids = collect(keys(prev_centroids))
    n, m = length(new_ids), length(prev_ids)

    costs = Array{Float64}(undef, n, m)
    for (i, nid) in enumerate(new_ids), (j, pid) in enumerate(prev_ids)
        costs[i, j] = GO.distance(new_centroids[nid], prev_centroids[pid])
    end

    assignment, _ = hungarian(costs)

    return Dict{Int64,Int64}(
        new_ids[i] => prev_ids[assignment[i]] for i in 1:length(new_ids)
    )
end

"""
    calculate_cluster_centroids(
    cluster_raster::Raster{Int64, 2};
    cluster_ids=[]
)::Vector{Cluster}

Calculate the centroids of the clusters in the raster.

# Arguments
- `clusters_raster`: Raster containing the cluster IDs.
- `cluster_ids`: Optional list of cluster IDs to assign to the clusters.

# Returns
A vector of `Cluster` objects.
"""
function calculate_cluster_centroids(
    clusters_raster::Raster{Int64, 2};
    cluster_ids=[]
)::Vector{Cluster}
    unique_clusters = Vector{Int64}(undef, 0)
    valid_clusters = unique(clusters_raster[clusters_raster .!= clusters_raster.missingval])

    unique_clusters = isempty(cluster_ids) ? sort(valid_clusters) : cluster_ids

    if !isempty(cluster_ids) && length(cluster_ids) != length(valid_clusters)
        error("Length of cluster IDs given do not match number of clusters in raster.")
    end

    clusters_vector = Vector{Cluster}(undef, length(unique_clusters))

    x_coords = clusters_raster.dims[1]
    y_coords = clusters_raster.dims[2]

    for (id, ex_id) in enumerate(unique_clusters)
        node_indices = findall(==(ex_id), clusters_raster)
        nodes = [(x_coords[i[1]], y_coords[i[2]]) for i in node_indices]
        col_cent = mean([node[1] for node in nodes])
        row_cent = mean([node[2] for node in nodes])

        clusters_vector[id] = Cluster(
            id = ex_id,
            centroid = Point{2, Float64}(col_cent, row_cent),
            nodes = Point{2, Float64}.(nodes)
        )
    end
    return clusters_vector
end

"""
    generate_target_clusters(
        target_points::Raster{Int, 2},
        clustered_targets_path::String,
        k::Int8,
        cluster_tolerance::Float64,
        current_location::Point{2, Float64},
        exclusions::DataFrame;
    )::Vector{Cluster}

Generate a vector of clustered targets from a raster of target locations, using kmeans.

# Arguments
- `target_points`: Raster containing the target points.
- `clustered_targets_path`: Path to the clustered targets raster.
- `k`: Number of clusters to create.
- `cluster_tolerance`: Tolerance for kmeans convergence.
- `current_location`: Current location of the mothership.
- `exclusions`: DataFrame containing the exclusion zones.

# Returns
A vector of Cluster objects.
"""
function generate_target_clusters(
    target_points::Raster{Int, 2},
    clustered_targets_path::String,
    k::Int8,
    cluster_tolerance::Float64,
    current_location::Point{2, Float64},
    exclusions::DataFrame;
)::Vector{Cluster}
    clustered_targets::Raster{Int64, 2} = apply_kmeans_clustering(
        target_points, k, current_location, exclusions; tol=cluster_tolerance
    )
    write(clustered_targets_path, clustered_targets; force = true)

    return calculate_cluster_centroids(clustered_targets)
end
