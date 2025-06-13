
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
    total_tender_capacity::Int = problem.tenders.capacity * problem.tenders.number

    dist_vector = get_feasible_distances(
        current_location,
        points,
        exclusions
    )

    feasible_idxs = findall(.!isinf.(dist_vector))
    feasible_points = points[feasible_idxs]

    # 3D coordinate matrix of feasible points for clustering
    coordinates_array = Matrix{Float64}(undef, 3, length(feasible_points))
    coordinates_array[1, :] .= getindex.(feasible_points, 1)
    coordinates_array[2, :] .= getindex.(feasible_points, 2)
    coordinates_array[3, :] .= dist_weighting .* dist_vector[feasible_idxs]

    clustering_assignments = capacity_constrained_kmeans(
        coordinates_array;
        max_cluster_size=total_tender_capacity,
    )

    clustered_targets_df::DataFrame = DataFrame(
        id = clustering_assignments,
        geometry = feasible_points
    )

    clustered_targets::Vector{Cluster} = calculate_cluster_centroids(clustered_targets_df)

    return clustered_targets
end

"""
    disturb_remaining_clusters(
        raster::Raster{Float64, 2},
        k::Int8,
        current_location::Point{2, Float64},
        exclusions::DataFrame;
        tol::Float64=1.0,
        dist_weighting::Float64=2E-5
    )::Raster{Int64, 2}
    disturb_remaining_clusters(
        unvisited_pts::DataFrame,
        k::Int8,
        current_location::Point{2, Float64},
        exclusions::DataFrame,
        total_tender_capacity::Int;
        tol::Float64=1.0,
        dist_weighting::Float64=2E-5
    )::DataFrame

- Disturb remaining clusters by simulating a disturbance event to remove nodes.
- Re-cluster the remaining nodes into `k` clusters.

# Arguments
- `raster`: Raster containing the target geometries.
- `unvisited_pts`: DataFrame containing the disturbance values for each node.
- `k`: Number of clusters to create.
- `current_location`: Current location of the mothership.
- `exclusions`: DataFrame containing the exclusion zones.
- `total_tender_capacity`: Total capacity of the tender fleet.
- `tol`: Tolerance for kmeans convergence.
- `dist_weighting`: Weighting factor for the distances (in kms) to be stored in 3d array,
    compared against lat/lons.

# Returns
A raster/DataFrame containing new, disturbed clusters. Cluster ID is assigned to each target
    location.
"""
function disturb_remaining_clusters(
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

    # 3D coordinate matrix for disturbance clustering
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
function disturb_remaining_clusters(
    unvisited_pts::DataFrame,
    k::Int,
    current_location::Point{2, Float64},
    exclusions::DataFrame,
    total_tender_capacity::Int;
    tol::Float64=1.0,
    dist_weighting::Float64=2E-5
)::DataFrame
    n_sites::Int = size(unvisited_pts, 1) # number of target sites remaining

    if n_sites <= k
        @warn "No disturbance, as (deployment targets <= clusters required)"
        return DataFrame()
    end

    # 3D coordinate matrix for disturbance clustering
    coordinates_3d = Matrix{Float64}(undef, 3, n_sites)
    coordinates_3d[1, :] .= getindex.(unvisited_pts.node, 1)
    coordinates_3d[2, :] .= getindex.(unvisited_pts.node, 2)
    coordinates_3d[3, :] .= unvisited_pts.disturbance_value

    # Create k_d clusters to create disturbance on subset
    k_d_lower = k+1
    k_d_upper = n_sites
    k_d = rand(k_d_lower:k_d_upper)

    disturbance_clusters = kmeans(
        coordinates_3d,
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
            mean(coordinates_3d[3, disturbance_clusters.assignments .== i])
            for i in 1:k_d
        ] .+
        t * rand(-1.0:0.01:1.0, k_d)

    # Assign the disturbance value to every node in the cluster
    disturbance_scores .= cluster_disturbance_vals[disturbance_clusters.assignments]

    # remove nodes in the cluster with the highest disturbance score
    max_disturbance_score = maximum(disturbance_scores)
    surviving_mask = disturbance_scores .!= max_disturbance_score

    disturbed_coordinates_2d = coordinates_3d[1:2, surviving_mask]
    n_sites = sum(surviving_mask)

    if k > n_sites
        #! Too many nodes/clusters removed! Change threshold,
        #! or use a different method e.g. remove cluster with highest scores
        error(
            "Too many nodes removed!\n$(n_sites) remaining node/s, $k clusters required."
        )
    end

    remaining_pts = Point{2,Float64}.(
        disturbed_coordinates_2d[1, :],
        disturbed_coordinates_2d[2, :]
    )

    # Fill vector with feasible distance from depot to each target site
    dist_vector = get_feasible_distances(current_location, remaining_pts, exclusions)

    # Mask out infeasible points
    feasible_idxs = findall(x -> x != Inf, dist_vector)
    feasible_pts = remaining_pts[feasible_idxs]

    disturbed_coordinates_3d = Matrix{Float64}(undef, 3, length(feasible_pts))
    disturbed_coordinates_3d[1, :] .= getindex.(feasible_pts, 1)
    disturbed_coordinates_3d[2, :] .= getindex.(feasible_pts, 2)
    disturbed_coordinates_3d[3, :] .= dist_weighting .* dist_vector[feasible_idxs]

    #re-cluster the remaining nodes into k clusters
    clustering_assignments = capacity_constrained_kmeans(
        disturbed_coordinates_3d;
        max_cluster_size = total_tender_capacity,
        k_spec = k,
    )

    return DataFrame(
        id = clustering_assignments,
        geometry = feasible_pts
    )
end

"""
    capacity_constrained_kmeans(
        coordinates::Matrix{Float64};
        max_cluster_size::Int64,
        max_split_distance::Int64 = 12000,
        k_spec::Int = 0,
        max_iter::Int64 = 1000,
        n_restarts::Int64 = 20
    )::Vector{Int64}

Cluster locations, ensuring that no cluster has more than `max_cluster_size`, and all points
are assigned to a cluster.

# Arguments
- `coordinates`: A matrix of coordinates where each column represents a reef location, and
    each row is a coordinate, i.e. longitude, latitude, (optionally distance).
- `max_cluster_size`: The maximum number of reefs per cluster.
- `max_split_distance`: The maximum distance (m) between clusters to allow for splitting.
- `k_spec`: The specified number of clusters to create. If 0, it will be calculated based on
    the number of reefs and `max_cluster_size`, allowing more clusters to be spawned.
- `max_iter`: The maximum number of iterations to run the k-means algorithm.
- `n_restarts`: The number of times to run k-means with different initial centroids.

# Returns
A vector of cluster assignments for each reef.
"""
function capacity_constrained_kmeans(
    coordinates::Matrix{Float64};
    max_cluster_size::Int64,
    max_split_distance::Int64 = 12000,
    k_spec::Int = 0,
    max_iter::Int64 = 1000,
    n_restarts::Int64 = 20,
)::Vector{Int64}
    n_reefs::Int64 = size(coordinates, 2)
    k = k_spec == 0 ? Ref(ceil(Int, n_reefs/max_cluster_size)) : Ref(k_spec)

    # Run k-means multiple times to find best result
    best_clustering_assignment::Vector{Int} = zeros(Int, n_reefs)
    best_score::Float64 = Inf
    for _ in 1:n_restarts
        # Reset k every time
        clustering_assignment::Vector{Int64} = _constrained_kmeans_single_iteration(
            coordinates,
            k,
            k_spec,
            max_cluster_size,
            max_split_distance,
            max_iter
        )
        k[] = k_spec == 0 ? maximum(clustering_assignment) : k_spec
        clusters_list::Vector{Vector{Int64}} = findall.(
            .==(1:k[]), Ref(clustering_assignment)
        )
        centroids = Tuple(
            Tuple(mean(coordinates[1:2,c]; dims=2))
            for c in clusters_list
        )
        cluster_score::Float64 = sum(
            haversine(coordinates[1:2, p], centroids[i])
            for i in 1:k[] for p in clusters_list[i]
        )
        if cluster_score < best_score
            best_score, best_clustering_assignment = cluster_score, clustering_assignment
        end
    end

    return best_clustering_assignment
end

"""
    _constrained_kmeans_single_iteration(
        coordinates::Matrix{Float64},
        k::Ref{Int},
        k_spec::Int = 0,
        max_cluster_size::Int64 = 6,
        max_split_distance::Int64 = 12000,
        max_iter::Int64 = 1000,
    )::Vector{Int64}

Run a single iteration of k-means clustering with constraints on cluster size and distance
    between clusters.

# Arguments
- `coordinates`: A matrix of coordinates where each column represents a reef's
    longitude and latitude (and optionally a third dimension for distance).
- `k`: A reference to the number of clusters to create.
- `k_spec`: The specified number of clusters to create. If 0, it will be calculated based on
    the number of reefs and `max_cluster_size`, allowing more clusters to be spawned.
- `max_cluster_size`: The maximum number of reefs per cluster.
- `max_split_distance`: The maximum distance (m) between clusters to allow for splitting.
- `max_iter`: The maximum number of iterations to run the k-means algorithm.

# Returns
A vector of cluster assignments for each reef, ensuring that no cluster exceeds the
    maximum size and that points are reassigned to the closest cluster within the distance
    constraints.
"""
function _constrained_kmeans_single_iteration(
    coordinates::Matrix{Float64},
    k::Ref{Int},
    k_spec::Int = 0,
    max_cluster_size::Int64 = 6,
    max_split_distance::Int64 = 12000,
    max_iter::Int64 = 1000,
)::Vector{Int64}
    clustering = kmeans(coordinates, k[]; maxiter=max_iter)
    clustering_assignment::Vector{Int64} = copy(clustering.assignments)

    # Build clusters & centroids
    clusters_list::Vector{Vector{Int64}} = findall.(
        .==(1:k[]),
        Ref(clustering_assignment)
    )
    centroids = Tuple(
        Tuple(mean(coordinates[1:2,c]; dims=2))
        for c in clusters_list
    )

    available_clusters = Vector{Int64}(undef, k[])
    dists_pt_to_centroids = Vector{Float64}(undef, k[])

    # Enforce max cluster size by reassigning furthest points for over-capacity clusters
    for c in 1:k[]
        point_idxs::Vector{Int64} = clusters_list[c]
        while length(point_idxs) > max_cluster_size
            # Find furthest point from centroid
            dists_centroid_to_pts::Vector{Float64} = [
                haversine(coordinates[1:2,p], centroids[c])
                for p in point_idxs
            ]
            idx::Int64 = point_idxs[argmax(dists_centroid_to_pts)]

            # Find under-capacity clusters within max_split_distance
            available_clusters = findall(length.(clusters_list) .< max_cluster_size)
            dists_pt_to_centroids = [
                haversine(coordinates[1:2,idx], centroids[c])
                for c in available_clusters
            ]
            close_clusters::Vector{Int64} = available_clusters[
                dists_pt_to_centroids .≤ [max_split_distance]
            ]


            if isempty(close_clusters)
                if iszero(k_spec)
                    # no available AND close clusters --> create new cluster
                    k[] += 1
                    return _constrained_kmeans_single_iteration(
                        coordinates,
                        k,
                        max_cluster_size,
                        max_split_distance,
                        max_iter
                    )
                else
                    close_clusters = available_clusters
                end
            end

            # Pick the closest among them
            eligible_centroids = centroids[close_clusters]
            eligible_distances::Vector{Float64} = [
                haversine(coordinates[1:2,idx], c)
                for c in eligible_centroids
            ]
            target_cluster::Int64 = close_clusters[argmin(eligible_distances)]

            # Reassign point
            clustering_assignment[idx] = target_cluster
            deleteat!(point_idxs, findfirst(==(idx), point_idxs))
            push!(clusters_list[target_cluster], idx)
        end
    end
    return clustering_assignment
end

"""
    update_cluster_assignments(
        cluster_raster::Raster{Int64, 2},
        prev_centroids::Dict{Int64, Point{2,Float64}}
    )::Raster{Int64, 2}
    update_cluster_assignments(
        cluster_df::DataFrame,
        prev_centroids::Dict{Int64, Point{2,Float64}}
    )::DataFrame

Update cluster assignments in the raster to match previous cluster numbering, based on
matching closest clster centroids.

# Arguments
- `cluster_raster`: Raster containing the new cluster IDs.
- `cluster_df`: DataFrame containing the new cluster IDs and their geometries.
- `prev_centroids`: Dictionary mapping previous cluster IDs to their centroids.

# Returns
A new raster or DataFrame with updated cluster assignments to match the previous numbering.
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
function update_cluster_assignments(
    cluster_df::DataFrame,
    prev_centroids::Dict{Int64, Point{2,Float64}}
)::DataFrame
    unique_new = Set(cluster_df.id)

    # Compute new centroids from the cluster_raster
    new_centroids = Dict{Int64, Point{2,Float64}}()
    for new_id in unique_new
        indices = findall(==(new_id), cluster_df.id)
        mean_lon = mean(getindex.(cluster_df.geometry[indices], 1))
        mean_lat = mean(getindex.(cluster_df.geometry[indices], 2))
        new_centroids[new_id] = (mean_lon, mean_lat)
    end

    # Map cluster IDs from new to previous clusters
    cluster_mapping = one_to_one_mapping_hungarian(new_centroids, prev_centroids)
    new_ids = Vector{Int64}(undef, length(cluster_df.id))

    for (i, new_id) in enumerate(cluster_df.id)
        new_ids[i] = cluster_mapping[new_id]
    end
    updated_df = DataFrame(
        id = new_ids,
        geometry = cluster_df.geometry
    )
    return updated_df
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

Calculate the centroids of the clusters in the raster/DataFrame, and return a vector of
    `Cluster` objects.

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

    for (idx, ex_id) in enumerate(unique_clusters)
        node_indices = findall(==(ex_id), clusters_raster)
        nodes = [(x_coords[i[1]], y_coords[i[2]]) for i in node_indices]
        col_cent = mean([node[1] for node in nodes])
        row_cent = mean([node[2] for node in nodes])

        clusters_vector[idx] = Cluster(
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

    for (idx, ex_id) in enumerate(unique_clusters)
        cluster_mask = clusters.id .== ex_id
        clustered_points = clusters.geometry[cluster_mask]

        lon_cent = mean(getindex.(clustered_points, 1))
        lat_cent = mean(getindex.(clustered_points, 2))

        clusters_vector[idx] = Cluster(
            id = ex_id,
            centroid = Point{2, Float64}(lon_cent, lat_cent),
            nodes = clustered_points
        )
    end
    return clusters_vector
end
