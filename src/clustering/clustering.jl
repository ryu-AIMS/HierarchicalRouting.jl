
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
        raster::Raster{Int, 2}, k::Int8; tol::Float64=1.0
    )::Raster{Int64, 2}
    apply_kmeans_clustering(
        raster::Raster{Float64, 2}, k::Int8; tol::Float64=1.0
    )::Raster{Int, 2}

Cluster targets sites by applying k-means to target (non-zero) cells in a raster.
- Float64 raster is assumed to contain disturbance values, addressed buy 3d clustering
- Int raster is assumed to contain cluster ID values, addressed by 2d spatial clustering.

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
    indices::Vector{CartesianIndex{2}} = findall(x -> x != raster.missingval, raster)
    n::Int = length(indices)
    coordinates_array = Matrix{Float64}(undef, 2, n)

    # Row/col indices from each CartesianIndex
    rows::Vector{Int64} = getindex.(indices, 1)
    cols::Vector{Int64} = getindex.(indices, 2)

    # Fill the coordinate matrix using the corresponding dimension arrays
    coordinates_array[1, :] .= raster.dims[1][rows]
    coordinates_array[2, :] .= raster.dims[2][cols]

    clustering = kmeans(coordinates_array, k; tol=tol, rng=Random.seed!(1))

    clustered_targets::Raster{Int64, 2} = similar(raster, Int64)
    clustered_targets[indices] .= clustering.assignments

    return clustered_targets
end
function apply_kmeans_clustering(
    raster::Raster{Float64, 2}, k::Int8; tol::Float64=1.0
)::Raster{Int64, 2}
    # TODO: Split this function into two separate functions, it does more than original fn
    #! 1st: 3D clustering to generate disturbance clusters,
    #! 2nd: 2D clustering to generate target clusters (as above)
    indices::Vector{CartesianIndex{2}} = findall(!=(raster.missingval), raster)
    n::Int = length(indices) # number of target sites remaining

    if n <= k
        @warn "No disturbance, as (deployment targets <= clusters required)"
        empty_raster = similar(raster, Int64, missingval=0)
        return empty_raster
    end

    # 2D coordinate matrix for clustering
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

    # remove nodes with the highest disturbance score - i.e. one cluster
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

    #re-cluster the remaining nodes into k clusters
    clustering = kmeans(
        coordinates_array_2d_disturbed,
        k;
        tol=tol,
        rng=Random.seed!(1)
    )

    clustered_targets = similar(raster, Int64, missingval=0)
    clustered_targets[indices] .= clustering.assignments

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
    cluster_mapping = Dict{Int64, Int64}()
    for new_id in unique_new
        new_centroid = new_centroids[new_id]
        prev_closest = nothing
        best_distance = Inf
        # Iterate over the previous clusters.
        for (prev_id, prev_centroid) in prev_centroids
            dist = GO.distance(new_centroid, prev_centroid)
            if dist < best_distance
                best_distance = dist
                prev_closest = prev_id
            end
        end
        cluster_mapping[new_id] = prev_closest
    end

    #? default args for similar?
    updated_raster = similar(cluster_raster, Int64, missingval=cluster_raster.missingval)
    for new_id in unique_new
        updated_raster[cluster_raster .== new_id] .= cluster_mapping[new_id]
    end

    return updated_raster
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
    targets::Targets,
    suitable_threshold::Float64,
    target_subset_path::String,
    subset::DataFrame,
    EPSG_code::Int16
)::Vector{Cluster}
    cluster_raster = process_targets(
        targets,
        k,
        cluster_tolerance,
        suitable_threshold,
        target_subset_path,
        subset,
        EPSG_code
    )

    write(clustered_targets_path, cluster_raster; force = true)

    return calculate_cluster_centroids(cluster_raster)
end

"""
    process_targets(
        targets::Targets,
        k::Int8,
        cluster_tolerance::Float64,
        suitable_targets_all_path::String,
        suitable_threshold::Float64,
        target_subset_path::String,
        subset::DataFrame,
        EPSG_code::Int16
    )::Raster{Int64}

Generate a clustered targets raster by reading in the suitable target location data,
applying thresholds and cropping to a target subset, and then clustering.

# Arguments
- `targets`: The targets object containing the target geometries.
- `k`: The number of clusters.
- `cluster_tolerance`: The cluster tolerance.
- `suitable_threshold`: The suitable targets threshold.
- `target_subset_path`: The path to the target subset raster.
- `subset`: The DataFrame containing the study area boundary.
- `EPSG_code`: The EPSG code for the study area.

# Returns
The clustered targets raster, classified by cluster ID number.
"""
function process_targets(
    targets::Targets,
    k::Int8,
    cluster_tolerance::Float64,
    suitable_threshold::Float64,
    target_subset_path::String,
    subset::DataFrame,
    EPSG_code::Int16
)::Raster{Int64}
    if endswith(targets.path, ".geojson")
        suitable_targets_all = process_geometry_targets(
            targets.gdf.geometry,
            EPSG_code
        )
    else
        suitable_targets_all = process_raster_targets(
            targets,
            EPSG_code,
            suitable_threshold
        )
    end

    suitable_targets_subset::Raster{Int} = Rasters.crop(suitable_targets_all, to=subset.geom)
    if !isfile(target_subset_path)
        write(target_subset_path, suitable_targets_subset; force=true)
    end

    clustered_targets::Raster{Int64, 2} = apply_kmeans_clustering(
        suitable_targets_subset, k; tol=cluster_tolerance
    )
    return clustered_targets
end
