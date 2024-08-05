using GLMakie, GeoMakie


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

function plot_mothership_route(clustered_targets::Raster{Int16, 2}, cluster_centroids::DataFrame, cluster_sequence::DataFrame)
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], title = "Mothership Route", xlabel = "Longitude", ylabel = "Latitude")

    image!(ax, clustered_targets, colormap = :viridis)

    scatter!(ax, cluster_centroids.lon, cluster_centroids.lat, markersize = 5, color = :red, label = "Cluster Centroids")

    # Annotate the cluster centroids with their cluster_id
    for i in 1:nrow(cluster_centroids)
        text!(ax, cluster_centroids.lon[i]+0.001, cluster_centroids.lat[i], text = string(cluster_centroids.cluster_id[i]), align = (:center, :center), color = :black)
    end

    # Generate the mothership route from the cluster sequence
    mothership_route = [(cluster_centroids.lon[findfirst(==(id), cluster_centroids.cluster_id)], cluster_centroids.lat[findfirst(==(id), cluster_centroids.cluster_id)]) for id in cluster_sequence.cluster_id]

    # Extract latitude and longitude from mothership route
    route_lats = last.(mothership_route)  # [wp[2] for wp in mothership_route]
    route_lons = first.(mothership_route)  # [wp[1] for wp in mothership_route]

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
