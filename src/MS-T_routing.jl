import ArchGDAL as AG
using Rasters
using DataFrames
using Statistics
using Distances
using Clustering

"""
    extract_subset(spatial_dataset::Raster, subset::GeoDataFrame)

Extract a subset of a raster dataset based on a GeoDataFrame.
"""
function extract_subset(spatial_dataset::Raster, subset)
    result_raster = Rasters.trim(mask(spatial_dataset; with=subset.geom))
    return result_raster
end

"""
    cluster_targets(df::DataFrame, num_clust::Int64)

Cluster the targets in a GeoDataFrame based on their geometry.
"""
function cluster_targets(df::DataFrame, num_clust::Int64)
    # Calculate centroid of geometry for each row
    centroid_shp = [AG.centroid(row.geom) for row in eachrow(df)]
    centroid_coords = [(AG.getx(centroid,0), AG.gety(centroid,0)) for centroid in centroid_shp]

    # Convert the coordinates to a format suitable for clustering (e.g., an array)
    coordinates_array = hcat([collect(c) for c in centroid_coords]...)

    # Cluster centroids using kmeans
    clustering = kmeans(coordinates_array, num_clust)

    df.cluster_id = clustering.assignments
    return df
end
function cluster_targets(raster::Raster{Int16, 2}, num_clust::Int64)
    # Extract the coordinates of non-zero values
    indices = findall(x -> x != 0, raster)

    # Convert the indices to coordinates (tuples)
    coordinates = [(Tuple(index)[1], Tuple(index)[2]) for index in indices]
    # Convert the coordinates to a format suitable for clustering (e.g., an array)
    coordinates_array = hcat([collect(c) for c in coordinates]...)

    # Perform k-means clustering on the coordinates
    clustering = kmeans(coordinates_array, num_clust)

    # Create a DataFrame to store the cluster assignments
    rows = [coord[1] for coord in coordinates]
    cols = [coord[2] for coord in coordinates]
    # cluster_df = DataFrame(row = rows, col = cols, cluster_id = clustering.assignments)

    # Create a new raster to store the cluster IDs
    cluster_raster = copy(raster)
    cluster_raster .= 0  # Initialize with zeros

    # Assign the cluster IDs to the corresponding positions in the new raster
    for i in 1:length(rows)
        cluster_raster[rows[i], cols[i]] = clustering.assignments[i]
    end

    return cluster_raster#, cluster_df
end

function calc_cluster_centroids(cluster_raster::Raster{Int16, 2})
    # Extract the unique cluster IDs
    unique_clusters = unique(cluster_raster)

    # Remove the zero cluster ID (if present)
    unique_clusters = unique_clusters[unique_clusters .!= 0]

    # Initialize the DataFrame for cluster centroids
    cluster_centroids = DataFrame(cluster_id = Int[], lat = Float64[], lon = Float64[])

    # Get the latitude and longitude coordinates
    X_dim = cluster_raster.dims[1]
    Y_dim = cluster_raster.dims[2]

    # Extracting coordinates using DimensionalData
    coordinates = [(x, y) for x in X_dim, y in Y_dim]

    for cluster_id in unique_clusters
        # Find the indices of the current cluster_id
        indices = findall(x -> x == cluster_id, cluster_raster)

        # Calculate the mean coordinates for the cluster
        lat = mean([coordinates[i][2] for i in indices])
        lon = mean([coordinates[i][1] for i in indices])

        # Add the centroid to the DataFrame
        push!(cluster_centroids, (cluster_id, lat, lon))
    end

    return cluster_centroids
end

function distance_matrix(cluster_centroids::DataFrame, depot::Tuple{Float64, Float64})
    # Number of centroids
    num_centroids = nrow(cluster_centroids)

    # Initialize the distance matrix
    dist_matrix = Matrix{Float64}(undef, num_centroids + 1, num_centroids + 1)

    # Get the coordinates of the centroids
    centroid_coords = [(row.lat, row.lon) for row in eachrow(cluster_centroids)]

    # Compute distances from the depot to each centroid
    dist_matrix[1, 1] = 0.0  # depot to depot distance
    for i in 1:num_centroids
        dist_matrix[1, i + 1] = dist_matrix[i + 1, 1] = haversine(depot, centroid_coords[i])
        # haversine(depot, centroid_coords[i])
    end

    # Compute distances between centroids
    for i in 1:num_centroids
        for j in 1:num_centroids
            dist_matrix[i + 1, j + 1] = haversine(centroid_coords[i], centroid_coords[j])
        end
    end

    return dist_matrix
end

"""
    nearest_neighbour(dist_matrix::Matrix{Float64})

Apply the nearest neighbor algorithm starting from the depot (1st row/col) and returning to the depot.
"""
function nearest_neighbour(cluster_centroids::DataFrame, depot::Tuple{Float64, Float64})

    dist_matrix = distance_matrix(cluster_centroids, depot)

    num_clusters = size(dist_matrix, 1) - 1  # excludes the depot
    visited = falses(num_clusters + 1)
    tour = Int[]
    total_distance = 0.0

    # Start at the depot
    current_location = 1
    push!(tour, current_location)
    visited[current_location] = true

    while length(tour) <= num_clusters
        # Find the nearest unvisited neighbor
        distances = dist_matrix[current_location, :]
        distances[visited] .= Inf  # visited distances = inf
        nearest_idx = argmin(distances)

        # Update tour and total distance
        push!(tour, nearest_idx)
        total_distance += dist_matrix[current_location, nearest_idx]
        visited[nearest_idx] = true
        current_location = nearest_idx
    end

    # Return to the depot
    push!(tour, 1)
    total_distance += dist_matrix[current_location, 1]

    # Adjust cluster_sequence to zero-based indexing
    cluster_sequence = tour .- 1

    # Generate mothership waypoints
    cluster_centroids = vcat(DataFrame(cluster_id = [0], lat = [depot[1]], lon = [depot[2]]), sort!(cluster_centroids, :cluster_id))
    ordered_centroids = cluster_centroids[[findfirst(==(id), cluster_centroids.cluster_id) for id in cluster_sequence], :]
    mothership_waypoints = [(row.lat, row.lon) for row in eachrow(ordered_centroids)]

    mothership_dist = 0.0
    for i in 1:length(mothership_waypoints) - 1
        mothership_dist += haversine(mothership_waypoints[i], mothership_waypoints[i + 1])
    end

    return cluster_sequence, total_distance, mothership_waypoints, mothership_dist
end

function plot_mothership_route(clustered_targets::Raster{Int16, 2}, cluster_centroids::DataFrame, mothership_waypoints::Vector{Tuple{Float64, Float64}}, depot::Tuple{Float64, Float64})
    fig = Figure()
    ax = Axis(fig[1, 1], title = "Mothership Route", xlabel = "Longitude", ylabel = "Latitude")

    # Plot the clustered targets
    image!(ax, clustered_targets, colormap=:viridis)

    # Plot the cluster centroids
    scatter!(ax, cluster_centroids.lon, cluster_centroids.lat, color=:red, markersize=25, label="Cluster Centroids")

    # Plot the depot
    scatter!(ax, [depot[2]], [depot[1]], color=:black, markersize=25, label="Depot")

    # Extract latitude and longitude from mothership waypoints
    route_lats = [wp[1] for wp in mothership_waypoints]
    route_lons = [wp[2] for wp in mothership_waypoints]

    # Plot the mothership route
    lines!(ax, route_lons, route_lats, color=:black, linewidth=2, label="Mothership Route")

    axislegend(ax)
    display(fig)
end
