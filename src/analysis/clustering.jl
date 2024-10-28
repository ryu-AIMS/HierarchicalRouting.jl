
mutable struct ClusterFields
    centroid::Point{2, Float64}
    nodes::Vector{Point{2, Float64}}
    # TODO: Add waypoints
    # waypoints::NTuple{2, Point{2, Float64}}

    function ClusterFields(centroid::Point{2, Float64})
        new(centroid, [centroid])  # nodes defaults to [centroid]
    end

    function ClusterFields(centroid::Point{2, Float64}, nodes::Vector{Point{2, Float64}})
        new(centroid, nodes)
    end
end

struct Cluster
    id::Int
    attributes::ClusterFields
end

"""
    cluster_targets(raster::Raster{Int, 2}, num_clust::Int64)

Cluster the targets in a GeoDataFrame based on their geometry.

# Arguments
- `raster::Raster{Int, 2}` : Raster containing the target geometries.
- `num_clust::Int64` : Number of clusters to create.

# Returns
A DataFrame containing the cluster ID and the target geometries.
"""
function cluster_targets(raster::Raster{Int, 2}, num_clust::Int64)
    indices = findall(x -> x != 0, raster)
    coordinates = [(Tuple(index)[1], Tuple(index)[2]) for index in indices]
    coordinates_array = hcat([collect(c) for c in coordinates]...)

    clustering = kmeans(coordinates_array, num_clust)

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

    coords = [(x, y) for x in clusters.dims[1], y in clusters.dims[2]]

    cluster_vec = Cluster[]

    if depot !== nothing
        push!(cluster_vec, Cluster(0, ClusterFields(depot, [depot])))
    end

    # Push Cluster object to cluster centroid vector
    for id in unique_clusters
        nodes = [(i[1], i[2]) for i in findall(==(id), clusters)]

        # lon = mean([coords[i[1], i[2]][1] for i in nodes])
        # lat = mean([coords[i[1], i[2]][2] for i in nodes])
        row_cent = mean([node[1] for node in nodes])
        col_cent = mean([node[2] for node in nodes])

        # push!(cluster_vec, Cluster(id, ClusterFields(Point{2, Float64}(lon, lat), [Point{2, Float64}(coords[i[1], i[2]][1], coords[i[1], i[2]][2]) for i in nodes])))
        push!(cluster_vec, Cluster(id, ClusterFields(Point{2, Float64}(row_cent, col_cent), [Point{2, Float64}(node[1], node[2]) for node in nodes])))
    end

    return cluster_vec
end
