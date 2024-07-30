using Statistics

import ArchGDAL as AG
import GeoInterface as GI
import GeometryOps as GO
using GeometryBasics
using CoordinateTransformations

using Rasters
using DataFrames
import GeoDataFrames as GDF

using Clustering, Distances

using GLMakie, GeoMakie


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

    return cluster_raster
end

function create_exclusion_zones(target_bathy::Raster, ms_depth)
    # Create exclusion zones based on the bathymetry data
    exclusion_zones = target_bathy .<= ms_depth
    return exclusion_zones
end

########

function plot_polygons(multipolygon::GI.Wrappers.MultiPolygon)
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], title = "Polygonized Raster Data")

    for polygon in GI.coordinates(multipolygon)
        for ring in polygon
            xs = [point[1] for point in ring]
            ys = [point[2] for point in ring]
            lines!(ax, xs, ys, color = :blue)
            # Optionally close the polygon by connecting the last point to the first
            lines!(ax, [xs..., xs[1]], [ys..., ys[1]], color = :blue)
        end
    end

    display(fig)
end

"""
    multipolygon_to_dataframe(multipolygon::GeoInterface.Wrappers.MultiPolygon)

Convert multipolygons to a GeoDataFrame.

# Arguments
- `multipolygon` :
"""
function multipolygon_to_dataframe(multipolygon::GeoInterface.Wrappers.MultiPolygon)
    data = Vector(undef, length(GeoInterface.coordinates(multipolygon)))

    # Extract coordinates from the MultiPolygon
    for (polygon_id, polygon) in enumerate(GeoInterface.coordinates(multipolygon))
        # Convert coordinates to a Polygon object
        rings = [Point(coord...) for coord in polygon[1]]
        poly = Polygon(rings)
        data[polygon_id] = (polygon_id, poly)
    end

    return DataFrame(data, [:polygon_id, :geometry])
end

function convert_raster_to_polygon(raster::Raster{Bool,2})
    multipolygon = polygonize(raster)
    multipolygon = GeoInterface.MultiPolygon(GeoInterface.coordinates(multipolygon))

    # Convert the MultiPolygon to a DataFrame
    multipolygon_df = multipolygon_to_dataframe(multipolygon)
    # plot_polygons(multipolygon)
    return multipolygon_df
end

########

function calc_cluster_centroids(cluster_raster::Raster{Int16, 2}, depot::Tuple{Float64, Float64})
    unique_clusters = unique(cluster_raster)

    # Remove the zero cluster ID (if present)
    unique_clusters = unique_clusters[unique_clusters .!= 0]

    cluster_centroids = DataFrame(cluster_id = Int[], lat = Float64[], lon = Float64[])

    # Get the latitude and longitude coordinates
    X_dim = cluster_raster.dims[1]
    Y_dim = cluster_raster.dims[2]

    # Extracting coordinates using DimensionalData
    coordinates = [(x, y) for x in X_dim, y in Y_dim]

    # Add depot as cluster 0
    push!(cluster_centroids, (0, depot[1], depot[2]))
    for cluster_id in unique_clusters
        # Find the indices of the current cluster_id
        indices = findall(x -> x == cluster_id, cluster_raster)

        # Calculate the mean coordinates for the cluster
        lat = mean([coordinates[i][2] for i in indices])
        lon = mean([coordinates[i][1] for i in indices])

        # Add the centroid to the DataFrame
        push!(cluster_centroids, (cluster_id, lat, lon))
    end
    # sort the cluster centroids by cluster_id
    sort!(cluster_centroids, :cluster_id)
    return cluster_centroids
end

function distance_matrix(cluster_centroids::DataFrame)
    # Number of centroids
    num_centroids = nrow(cluster_centroids)

    # Initialize the distance matrix
    dist_matrix = Matrix{Float64}(undef, num_centroids, num_centroids)

    # Get the coordinates of the centroids
    centroid_coords = [(row.lat, row.lon) for row in eachrow(cluster_centroids)]

    # Compute distances between centroids
    for i in 1:num_centroids
        for j in 1:num_centroids
            dist_matrix[i, j] = haversine(centroid_coords[i], centroid_coords[j])
        end
    end

    return dist_matrix
end

function is_feasible_path(start_pt::Tuple{Float64, Float64}, end_pt::Tuple{Float64, Float64}, env_constraint::DataFrame)
    # Create a line from start_pt to end_pt
    line = GI.LineString([start_pt, end_pt])

    # Check if the line intersects with any polygon in the env_constraint
    for row in eachrow(env_constraint)
        polygon = row.geometry
        if GI.intersects(line, polygon)
            return false
        end
    end
    return true
end
function apply_constraints(waypoints::Vector{Tuple{Float64, Float64}}, target_bathy::Raster{Int16, 2})
    # Initialize the total distance
    total_distance = 0.0

    # Iterate over the waypoints and calculate the distance between them
    for i in 1:length(waypoints)-1
        # Check if the path between the waypoints is feasible
        if !is_feasible_path(waypoints[i], waypoints[i+1], target_bathy)
            # find shortest path around constraint
        else
            # Calculate the distance between the waypoints
            distance = haversine(waypoints[i], waypoints[i+1])
        end
        total_distance += distance
    end

    return total_distance
end

########

"""
    nearest_neighbour(dist_matrix::Matrix{Float64})

Apply the nearest neighbor algorithm starting from the depot (1st row/col) and returning to the depot.
"""
function nearest_neighbour(cluster_centroids::DataFrame)

    dist_matrix = distance_matrix(cluster_centroids)

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
    # cluster_centroids = vcat(DataFrame(cluster_id = [0], lat = [depot[1]], lon = [depot[2]]), sort!(cluster_centroids, :cluster_id))
    ordered_centroids = cluster_centroids[[findfirst(==(id), cluster_centroids.cluster_id) for id in cluster_sequence], :]
    # mothership_waypoints = [(row.lat, row.lon) for row in eachrow(ordered_centroids)]

    return ordered_centroids, total_distance, dist_matrix
end

function plot_mothership_route(clustered_targets::Raster{Int16, 2}, cluster_centroids::DataFrame, cluster_sequence::DataFrame)
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], title = "Mothership Route", xlabel = "Longitude", ylabel = "Latitude")

    # Plot the clustered targets
    image!(ax, clustered_targets, colormap = :viridis)

    # Plot the cluster centroids
    scatter!(ax, cluster_centroids.lon, cluster_centroids.lat, markersize = 5, color = :red, label = "Cluster Centroids")

    # Annotate the cluster centroids with their cluster_id
    for i in 1:nrow(cluster_centroids)
        text!(ax, cluster_centroids.lon[i]+0.001, cluster_centroids.lat[i], text = string(cluster_centroids.cluster_id[i]), align = (:center, :center), color = :black)
    end

    # Generate the mothership route from the cluster sequence
    mothership_route = [(cluster_centroids.lon[findfirst(==(id), cluster_centroids.cluster_id)], cluster_centroids.lat[findfirst(==(id), cluster_centroids.cluster_id)]) for id in cluster_sequence.cluster_id]

    # Extract latitude and longitude from mothership route
    route_lats = [wp[2] for wp in mothership_route]
    route_lons = [wp[1] for wp in mothership_route]

    # Plot the mothership route
    lines!(ax, route_lons, route_lats, color = :blue, linewidth = 2, label = "Mothership Route")

    axislegend(ax)
    display(fig)
end
function plot_mothership_route(clustered_targets::Raster{Int16, 2}, waypoints::Vector{Tuple{Float64, Float64}}, cluster_sequence::DataFrame)
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], title = "Mothership Route", xlabel = "Longitude", ylabel = "Latitude")

    # Plot the clustered targets
    image!(ax, clustered_targets, colormap = :viridis)

    # Extract latitude and longitude from waypoints
    waypoint_lons = [wp[2] for wp in waypoints]
    waypoint_lats = [wp[1] for wp in waypoints]

    # Plot the cluster centroids
    scatter!(ax, waypoint_lons, waypoint_lats, markersize = 5, color = :red, label = "Cluster Centroids")

    # Annotate the cluster centroids with their cluster_id
    for i in 1:length(waypoints)
        text!(ax, waypoint_lons[i] + 0.001, waypoint_lats[i], text = string(i), align = (:center, :center), color = :black)
    end

    # Generate the mothership route including the depot at the start and end
    mothership_route = [(depot[2], depot[1])]  # Start with depot
    append!(mothership_route, [(wp[2], wp[1]) for wp in waypoints])  # Add waypoints
    push!(mothership_route, (depot[2], depot[1]))  # End with depot

    # Extract latitude and longitude from mothership route
    route_lats = [wp[2] for wp in mothership_route]
    route_lons = [wp[1] for wp in mothership_route]

    # Plot the mothership route
    lines!(ax, route_lons, route_lats, color = :blue, linewidth = 2, label = "Mothership Route")

    axislegend(ax)
    display(fig)
end

function calc_waypoints(cluster_centroids::DataFrame, cluster_sequence::DataFrame)::Vector{Tuple{Float64, Float64}}
    # Ensure cluster_sequence is a vector of cluster IDs
    cluster_ids = cluster_sequence.cluster_id

    # Generate waypoints for each cluster in the sequence
    waypoints = Vector{Tuple{Float64, Float64}}()
    # # Add the depot as the first waypoint
    # push!(waypoints, (cluster_centroids.lat[1], cluster_centroids.lon[1]))
    for i in 1:length(cluster_ids)-1
        current_centroid = cluster_centroids[findfirst(x -> x == cluster_ids[i], cluster_centroids.cluster_id), :]
        next_centroid = cluster_centroids[findfirst(x -> x == cluster_ids[i + 1], cluster_centroids.cluster_id), :]
        # Calculate the centroid between the current and next cluster
        centroid = ((current_centroid.lat + next_centroid.lat) / 2, (current_centroid.lon + next_centroid.lon) / 2)
        push!(waypoints, centroid)
    end
    # # Add the depot as the last waypoint
    # push!(waypoints, (cluster_centroids.lat[1], cluster_centroids.lon[1]))
    return waypoints
end
