
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
    cluster_problem(
        problem::Problem,
        dist_weighting::Float64=5E-6
    )::Vector{Cluster}

Cluster the problem data into groups based on the target locations and the depot location.
The clustering is done using k-means clustering, and the centroids of the clusters are
    calculated.

# Arguments
- `problem`: The problem data.
- `dist_weighting`: Weighting factor for the distances in 3D clustering, used in combination
    with lat/lons at first 2 dimensions. Higher values will give more weight to distance
    from current location (depot). Default = 5E-6.

# Returns
Vector of clustered locations.
"""
function cluster_problem(
    problem::Problem,
    dist_weighting::Float64=5E-6
)::Vector{Cluster}
    points::Vector{Point{2, Float64}} = problem.targets.points.geometry
    current_location::Point{2, Float64} = problem.depot
    exclusions::DataFrame = problem.tenders.exclusion

    dist_vector = dist_weighting .* get_feasible_distances(
        current_location,
        points,
        exclusions
    )

    feasible_idxs = findall(.!isinf.(dist_vector))
    feasible_points = points[feasible_idxs]

    # 3D coordinate matrix of feasible points for clustering
    coordinates_array = Matrix{Float64}(undef, 2, length(feasible_points))
    coordinates_array[1, :] .= getindex.(feasible_points, 1)
    coordinates_array[2, :] .= getindex.(feasible_points, 2)
    # coordinates_array[3, :] = dist_vector[feasible_idxs]'

    points_df = DataFrame(
        LON=getindex.(feasible_points, 1),
        LAT=getindex.(feasible_points, 2),
        geometry = feasible_points
    )
    clustering_assignments = capacitated_kmeans(
        points_df;
        max_reef_number = 6,
        max_iter = 1000,
        n_restarts = 50,
    )

    clustered_targets_df::DataFrame = DataFrame(
        id = clustering_assignments,
        geometry = feasible_points
    )

    clustered_targets::Vector{Cluster} = calculate_cluster_centroids(clustered_targets_df)

    return clustered_targets
end

"""
    apply_kmeans_clustering(
        raster::Raster{Float64, 2},
        k::Int8,
        current_location::Point{2, Float64},
        exclusions::DataFrame;
        tol::Float64=1.0,
        dist_weighting::Float64=2E-5
    )::Raster{Int64, 2}

Cluster targets locations by applying k-means to target (non-zero) cells in a Float64 raster
    containing disturbance values.
Clustering considers feasible distances from the current location as a 3rd dimension, and
    excludes points in exclusion zones.

# Arguments
- `raster`: Raster containing the target geometries.
- `k`: Number of clusters to create.
- `current_location`: Current location of the mothership.
- `exclusions`: DataFrame containing the exclusion zones.
- `tol`: Tolerance for kmeans convergence.
- `dist_weighting`: Weighting factor for the distances (in kms) to be stored in 3d array,
    compared against lat/lons.

# Returns
A raster containing new clusters, with cluster IDs assigned to each target locations.
"""
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
    n_sites::Int = length(indices) # number of target sites remaining

    if n_sites <= k
        @warn "No disturbance, as (deployment targets <= clusters required)"
        empty_raster = similar(raster, Int64, missingval=0)
        empty_raster .= empty_raster.missingval
        return empty_raster
    end

    # 3D coordinate matrix for clustering
    coordinates_array_3d = Matrix{Float64}(undef, 3, n_sites)
    coordinates_array_3d[1, :] .= raster.dims[1][getindex.(indices, 1)]
    coordinates_array_3d[2, :] .= raster.dims[2][getindex.(indices, 2)]
    coordinates_array_3d[3, :] .= raster[indices]

    # Create k_d clusters to create disturbance on subset
    k_d_lower = min(n_sites, k+1)
    k_d_upper = min(max(k+1, n_sites, k^2), n_sites)
    k_d = rand(k_d_lower:k_d_upper)

    disturbance_clusters = kmeans(
        coordinates_array_3d,
        k_d;
        tol=tol,
        rng=Random.seed!(1)
    )

    # Create a score based on the disturbance values for each cluster
    disturbance_scores = Vector{Float64}(undef, n_sites)
    # Calculate the mean disturbance value for each cluster with stochastic perturbation
    w = 1.0 # weight for the environmental disturbance value
    t = 1.0 # perturbation weighting factor
    cluster_disturbance_vals =
        w * [
            mean(coordinates_array_3d[3, disturbance_clusters.assignments .== i])
            for i in 1:k_d
        ] .+
        t * rand(-1.0:0.01:1.0, k_d)
    # Assign the disturbance value to every node in the cluster
    disturbance_scores .= cluster_disturbance_vals[disturbance_clusters.assignments]

    # remove nodes in the cluster with the highest disturbance score
    max_disturbance_score = maximum(disturbance_scores)
    surviving_mask = disturbance_scores .!= max_disturbance_score

    coordinates_array_2d_disturbed = coordinates_array_3d[1:2, surviving_mask]
    indices = indices[surviving_mask]
    n_sites = length(indices) # update number of target sites remaining

    if k > n_sites
        #! Too many nodes/clusters removed! Change threshold,
        #! or use a different method e.g. remove cluster with highest scores
        error(
            "Too many nodes removed!\n$(n_sites) remaining node/s, $k clusters required."
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
    capacitated_kmeans(
        reef_data;
        max_reef_number::Int = 6,
        max_split_distance::Float64 = 12.0,
        max_k::Int = 6,
        max_iter::Int = 1000,
        n_restarts::Int = 5
    )

Cluster locations, ensuring that no cluster has more than `max_reef_number`, and all points
are assigned to a cluster.

# Arguments
- `reef_data`: DataFrame with `.LAT` and `.LON` columns.
- `max_reef_number`: The maximum number of reefs per cluster.
- `max_split_distance`: The maximum distance between clusters to allow for splitting.
- `max_iter`: The maximum number of iterations to run the k-means algorithm.
- `n_restarts`: The number of times to run k-means with different initial centroids.

# Returns
A vector of cluster assignments for each reef.
"""
function capacitated_kmeans(
    reef_data;
    max_reef_number::Int = 6,
    max_split_distance::Float64 = 12.0,
    max_k::Int = 6,
    max_iter::Int = 1000,
    n_restarts::Int = 5,
)
    n_reefs = length(reef_data.LAT)
    k = Ref(ceil(Int, n_reefs/max_reef_number))
    coordinates_array = hcat(reef_data.LON, reef_data.LAT)' # 2×n for kmeans

    function quick_distance(i::Int, j::Int)
        if i == j
            return 0.0
        elseif i > j
            i, j = j, i
        end
        R = 6371.0
        lat1, lon1 = deg2rad(reef_data.LAT[i]), deg2rad(reef_data.LON[i])
        lat2, lon2 = deg2rad(reef_data.LAT[j]), deg2rad(reef_data.LON[j])
        dlat, dlon = (lat2 - lat1), (lon2 - lon1)
        a = sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2
        c = 2 * atan(sqrt(a), sqrt(1 - a))
        return R * c
    end
    function quick_distance(i::Int, (lon2, lat2)::Tuple{Float64, Float64})
        R = 6371.0
        lat1, lon1 = deg2rad(reef_data.LAT[i]), deg2rad(reef_data.LON[i])
        lat2, lon2 = deg2rad(lat2), deg2rad(lon2)
        dlat, dlon = (lat2 - lat1), (lon2 - lon1)
        a = sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2
        c = 2 * atan(sqrt(a), sqrt(1 - a))
        return R * c
    end

    function calc_centroid(cluster_indices)
        lon_sum = 0.0
        lat_sum = 0.0
        for i in cluster_indices
            lon_sum += reef_data.LON[i]
            lat_sum += reef_data.LAT[i]
        end
        return (lon_sum / length(cluster_indices), lat_sum / length(cluster_indices))
    end

    function single_run()
        clustering = kmeans(coordinates_array, k[]; maxiter=max_iter)
        clustering_assignment = copy(clustering.assignments)

        for _ in 1:max_iter
            # build clusters & centroids
            clusters = findall.(.==(1:k[]), Ref(clustering_assignment))
            centroids = calc_centroid.(clusters)

            # enforce max cluster size
            # for each over-capacity cluster, reassign its furthest points
            updated = false
            for c in 1:k[]
                point_idxs = clusters[c]
                while length(point_idxs) > max_reef_number
                    # find furthest point from centroid
                    dists = quick_distance.(point_idxs, Ref((centroids[c])))
                    idx = point_idxs[argmax(dists)]

                    # find under-capacity clusters within max_split_distance
                    available_clusters = findall(length.(clusters) .< max_reef_number)
                    # Main.@infiltrate
                    dists = quick_distance.(Ref(idx), centroids[available_clusters])
                    close_clusters = available_clusters[dists .≤ max_split_distance]

                    if isempty(close_clusters)
                        # no available AND close clusters --> create new cluster
                        k[] += 1
                        return single_run()
                    end

                    # pick the closest among them
                    eligible_centroids = centroids[close_clusters]
                    eligible_distances = quick_distance.(Ref(idx), eligible_centroids)
                    target_cluster = close_clusters[argmin(eligible_distances)]

                    # reassign point
                    clustering_assignment[idx] = target_cluster
                    deleteat!(point_idxs, findfirst(==(idx), point_idxs))
                    push!(clusters[target_cluster], idx)
                    updated = true
                end
                if updated
                    break
                end
            end

            updated || break
        end
        return clustering_assignment
    end

    # Run k-means multiple times to find best result
    best_clustering_assignment = zeros(Int, n_reefs)
    best_score = Inf
    for _ in 1:n_restarts
        # Reset k every time
        clustering_assignment = single_run()
        k[] = maximum(clustering_assignment) <= max_k ? maximum(clustering_assignment) : max_k
        clusters = findall.(.==(1:k[]), Ref(clustering_assignment))
        centroids = calc_centroid.(clusters)
        cluster_score = sum(
            [sum(quick_distance.(clusters[i], Ref(centroids[i]))) for i in 1:k[]]
        )
        if cluster_score < best_score
            best_score, best_clustering_assignment = cluster_score, clustering_assignment
        end
    end

    return best_clustering_assignment
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
    function calculate_cluster_centroids(
        clusters::DataFrame
    )::Vector{Cluster}

Calculate the centroids of the clusters in the raster.

# Arguments
- `clusters_raster`: Raster containing the cluster IDs.
- `clusters`: DataFrame containing the cluster IDs and their geometries.
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
function calculate_cluster_centroids(
    clusters::DataFrame
)::Vector{Cluster}
    unique_clusters = sort(unique(clusters.id))

    clusters_vector = Vector{Cluster}(undef, length(unique_clusters))

    for id in unique_clusters
        cluster_mask = clusters.id .== id
        clustered_points = clusters.geometry[cluster_mask]

        lon_cent = mean(getindex.(clustered_points, 1))
        lat_cent = mean(getindex.(clustered_points, 2))

        clusters_vector[id] = Cluster(
            id = id,
            centroid = Point{2, Float64}(lon_cent, lat_cent),
            nodes = clustered_points
        )
    end
    return clusters_vector
end
