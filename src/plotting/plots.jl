using GLMakie, GeoMakie


function plot_polygons(multipolygon::GI.Wrappers.MultiPolygon)
    fig = Figure(size = (800, 600))
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

function plot_mothership_route(clustered_targets::Raster{Int16, 2}, waypoints::Vector{Tuple{Float64, Float64}})
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], title = "Mothership Route", xlabel = "Longitude", ylabel = "Latitude")

    image!(ax, clustered_targets, colormap = :viridis)

    waypoint_lons = [wp[2] for wp in waypoints]
    waypoint_lats = [wp[1] for wp in waypoints]

    scatter!(ax, waypoint_lons, waypoint_lats, markersize = 5, color = :red, label = "Waypoints")

    for i in 1:length(waypoints)
        text!(ax, waypoint_lons[i] + 0.001, waypoint_lats[i], text = string(i), align = (:center, :center), color = :black)
    end

    # # Generate the mothership route including the depot at the start and end
    # mothership_route = [(depot[2], depot[1])]  # Start with depot
    # append!(mothership_route, [(wp[2], wp[1]) for wp in waypoints])  # Add waypoints
    # push!(mothership_route, (depot[2], depot[1]))  # End with depot

    # # Extract latitude and longitude from mothership route
    # route_lats = [wp[2] for wp in mothership_route]
    # route_lons = [wp[1] for wp in mothership_route]

    lines!(ax, waypoint_lons, waypoint_lats, color = :blue, linewidth = 2, label = "Mothership Route")

    axislegend(ax)
    display(fig)
end

function plot_route_w_exclusions(clustered_targets::Raster{Int16, 2}, waypoints::Vector{Point{2, Float64}}, ms_exclusion_zones::DataFrame)
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], title = "Mothership Route with Exclusions", xlabel = "Longitude", ylabel = "Latitude")

    image!(ax, clustered_targets, colormap = :viridis)

    waypoint_lons = [wp[2] for wp in waypoints]
    waypoint_lats = [wp[1] for wp in waypoints]

    scatter!(ax, waypoint_lons, waypoint_lats, markersize = 5, color = :red, label = "Waypoints")

    for i in 1:length(waypoints)
        text!(ax, waypoint_lons[i] + 0.001, waypoint_lats[i], text = string(i), align = (:center, :center), color = :black)
    end

    lines!(ax, waypoint_lons, waypoint_lats, color = :blue, linewidth = 2, label = "Mothership Route")

    for i in 1:nrow(ms_exclusion_zones)
        exclusion = ms_exclusion_zones[i, :]
        exclusion_polygon = exclusion.geometry
        for ring in GI.coordinates(exclusion_polygon)
            xs = [point[1] for point in ring]
            ys = [point[2] for point in ring]
            lines!(ax, xs, ys, color = :red)
            # Optionally close the polygon by connecting the last point to the first
            lines!(ax, [xs..., xs[1]], [ys..., ys[1]], color = :red)
        end
    end

    axislegend(ax)
    display(fig)

end
function plot_route_w_exclusions(clustered_targets::Raster{Int16, 2}, waypoints::Vector{Point{2, Float64}}, ms_exclusions::Raster{Bool, 2})
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], title = "Mothership Route with Exclusions", xlabel = "Longitude", ylabel = "Latitude")

    # cmap 0 to white
    custom_cmap = cgrad(:viridis, 256)
    custom_cmap = [RGBAf(1.0, 1.0, 1.0, 0.0); custom_cmap.colors[2:end]]  # Map 0 values to white/transparent

    heatmap!(ax, clustered_targets, colormap = custom_cmap)

    waypoint_lons = [wp[2] for wp in waypoints]
    waypoint_lats = [wp[1] for wp in waypoints]

    scatter!(ax, waypoint_lons, waypoint_lats, markersize = 5, color = :red, label = "Waypoints")

    for i in 1:length(waypoints)
        text!(ax, waypoint_lons[i] + 0.001, waypoint_lats[i], text = string(i), align = (:center, :center), color = :black)
    end

    lines!(ax, waypoint_lons, waypoint_lats, color = :blue, linewidth = 2, label = "Mothership Route")

    image!(ax, Raster(ms_exclusions .== 0), alpha = 0.5)

    axislegend(ax)
    display(fig)

end

function plot_waypoints_and_exclusions(waypoints::Vector{Point{2, Float64}}, exclusions::DataFrame)
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], title = "Mothership Route", xlabel = "Longitude", ylabel = "Latitude")

    waypoint_lons = [wp[1] for wp in waypoints]
    waypoint_lats = [wp[2] for wp in waypoints]

    # Plot cluster centroids
    scatter!(ax, waypoint_lons, waypoint_lats, markersize = 5, color = :red, label = "Cluster Centroids")

    # Annotate centroids with cluster_ids
    for i in 1:length(waypoints)
        text!(ax, waypoint_lons[i], waypoint_lats[i] + 0.01, text = string(i), align = (:center, :center), color = :black)
    end

    # MS route lat and longs
    route_lats = [wp[2] for wp in waypoints]
    route_lons = [wp[1] for wp in waypoints]

    # Exclusion zones
    for (i, zone) in enumerate(eachrow(exclusions))
        polygon = zone[:geometry]
        for ring in GeoInterface.coordinates(polygon)
            xs = [coord[1] for coord in ring]
            ys = [coord[2] for coord in ring]
            poly!(ax, xs, ys, color = (:gray, 0.5), strokecolor = :black, label = "Exclusion Zone")

            centroid_x, centroid_y = mean(xs), mean(ys)
            text!(ax, centroid_x, centroid_y, text = string(i), align = (:center, :center), color = :blue)
        end
    end

    # axislegend(ax)
    display(fig)
end

function plot_waypoints_and_exclusions_with_graph(g::SimpleWeightedGraph{Int64, Float64},
    idx_to_point::Dict{Int64, Point{2, Float64}},
    waypoints::Vector{Point{2, Float64}},
    exclusions::DataFrame)
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], title = "Mothership Route with Graph", xlabel = "Longitude", ylabel = "Latitude")

    # Waypoints
    waypoint_lons = [wp[1] for wp in waypoints]
    waypoint_lats = [wp[2] for wp in waypoints]
    scatter!(ax, waypoint_lons, waypoint_lats, markersize = 5, color = :red, label = "Waypoints")
    for i in 1:length(waypoints)
        text!(ax, waypoint_lons[i], waypoint_lats[i] + 0.01, text = string(i), align = (:center, :center), color = :black)
    end

    # Exclusion zones (polygons)
    for (i, zone) in enumerate(eachrow(exclusions))
        polygon = zone[:geometry]
        for ring in GeoInterface.coordinates(polygon)
            xs = [coord[1] for coord in ring]
            ys = [coord[2] for coord in ring]
            poly!(ax, xs, ys, color = (:gray, 0.5), strokecolor = :black, label = "Exclusion Zone")

            centroid_x, centroid_y = mean(xs), mean(ys)
            text!(ax, centroid_x, centroid_y, text = string(i), align = (:center, :center), color = :blue)
        end
    end

    # Weighted graph
    x_coords = [idx_to_point[i][1] for i in 1:nv(g)]
    y_coords = [idx_to_point[i][2] for i in 1:nv(g)]

    scatter!(ax, x_coords, y_coords, markersize = 6, color = :blue, label = "Graph Nodes")
    [text!(ax, x_coords[i], y_coords[i], text = string(i), align = (:center, :center), color = :red) for i in 1:nv(g)]

    for e in edges(g)
        x1, y1 = idx_to_point[e.src][1], idx_to_point[e.src][2]
        x2, y2 = idx_to_point[e.dst][1], idx_to_point[e.dst][2]
        lines!(ax, [x1, x2], [y1, y2], color = :green)
    end

    edge_weights = [g.weights[e.src, e.dst] for e in edges(g)]
    for (e, weight) in zip(edges(g), edge_weights)
        x1, y1 = idx_to_point[e.src][1], idx_to_point[e.src][2]
        x2, y2 = idx_to_point[e.dst][1], idx_to_point[e.dst][2]

        lines!(ax, [x1, x2], [y1, y2], color = :green)

        # # Add edge weights
        # mid_x, mid_y = (x1 + x2) / 2, (y1 + y2) / 2
        # text!(ax, mid_x, mid_y, text = string(round(weight, digits = 2)), align = (:center, :center), color = :purple)
    end

    display(fig)
end
