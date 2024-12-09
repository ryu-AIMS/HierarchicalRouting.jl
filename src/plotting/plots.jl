module Plot

import ..HierarchicalRouting

using Rasters
using DataFrames
using GeometryBasics

using GLMakie, GeoMakie


function points(clusters::Vector{HierarchicalRouting.Cluster})
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], xlabel = "Longitude", ylabel = "Latitude")

    points!(ax, clusters)

    return fig, ax
end
function points!(
    ax::Axis,
    clusters::Vector{HierarchicalRouting.Cluster};
)
    colormap = distinguishable_colors(length(clusters) + 1)[2:end]

    for (idx, cluster) in enumerate(clusters)
        color = colormap[idx]
        for point in cluster.nodes
            scatter!(ax, point, color = color, markersize = 5)
        end
    end
    return ax
end

function clusters(
    cluster_sequence::DataFrame;
    cluster_radius = 100,
    centers = false
)
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], xlabel = "Longitude", ylabel = "Latitude")

    clusters!(ax, cluster_sequence, cluster_radius = cluster_radius, centers = centers)

    return fig, ax
end
function clusters!(
    ax::Axis,
    cluster_sequence::DataFrame;
    cluster_radius = 100,
    centers = false
)
    sequence_id = [row.id for row in eachrow(cluster_sequence)[2:end-1]]
    centroids = hcat(cluster_sequence.lon, cluster_sequence.lat)
    # series(centroids)

    # cluster circles
    colormap = distinguishable_colors(length(sequence_id)+1)[2:end]

    # Precompute circle offsets
    circle_offset_lon = cluster_radius .* cos.(range(0, 2π, length=100))
    circle_offset_lat = cluster_radius .* sin.(range(0, 2π, length=100))

    for (idx, seq) in enumerate(sequence_id)
        color = colormap[seq]

        center_lon = centroids[:, 1][idx + 1]
        center_lat = centroids[:, 2][idx + 1]

        circle_lons = center_lon .+ circle_offset_lon
        circle_lats = center_lat .+ circle_offset_lat

        circle = hcat(circle_lons, circle_lats)

        poly!(ax, circle, color = (color, 0.2), strokecolor = color, label = "Cluster Centroids")

        if centers
            scatter!(ax, [center_lon], [center_lat], markersize = 10, color = (color, 0.2), strokewidth = 0)
            # series((center_lon, center_lat), color = (color, 0.2), markersize = 10)
        end
    end
    return ax
end

function exclusions(
    exclusions::DataFrame;
    labels::Bool = false
)
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], xlabel = "Longitude", ylabel = "Latitude")

    exclusions!(ax, exclusions, labels = labels)

    return fig, ax
end
function exclusions!(
    ax::Axis,
    exclusions::DataFrame;
    labels::Bool = false
)
    for (i, zone) in enumerate(eachrow(exclusions))
        polygon = zone[:geometry]
        for ring in GeoInterface.coordinates(polygon)
            xs, ys = [coord[1] for coord in ring], [coord[2] for coord in ring]
            poly!(ax, xs, ys, color = (:gray, 0.5), strokecolor = :black)#, label = "Exclusion Zone")

            if labels
                centroid_x, centroid_y = mean(xs), mean(ys)
                text!(ax, centroid_x, centroid_y, text = string(i), align = (:center, :center), color = :blue)
            end
        end
    end
    return ax
end

function linestrings(
    line_strings::Vector{LineString{2, Float64}};
    labels::Bool = false
)
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], xlabel = "Longitude", ylabel = "Latitude")

    linestrings!(ax, line_strings, labels = labels)

    return fig, ax
end
    function linestrings!(
    ax::Axis,
    line_strings::Vector{LineString{2, Float64}};
    labels::Bool = false
)
    n_graphs = length(line_strings)
    color_palette = [cgrad(:rainbow, n_graphs)[i] for i in 1:n_graphs]

    waypoints = [line[1][1] for line in line_strings]

    waypoint_lons, waypoint_lats = [wp[1] for wp in waypoints], [wp[2] for wp in waypoints]

    # Mark waypoints with 'x'
    scatter!(ax, waypoint_lons, waypoint_lats, marker = 'x', markersize = 10, color = :black)#, label = "Cluster Centroids")

    # Plot LineStrings
    for (idx, line_string) in enumerate(line_strings)
        color = color_palette[idx]
        points = [Point(p[1], p[2]) for l in line_string for p in l.points]
        lines!(ax, points, color = color, linewidth = 2)
    end

    if labels
        # Annotate waypoints by sequence
        text!(ax, waypoint_lons, waypoint_lats .+ 30, text = string.((0:length(waypoint_lons)-1)), align = (:center, :center), color = :black)
    end
    return ax
end

function tenders(
    tender_soln::Vector{HierarchicalRouting.TenderSolution}
)
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], xlabel = "Longitude", ylabel = "Latitude")

    tenders!(ax, tender_soln)

    return fig, ax
end
function tenders!(
    ax::Axis,
    tender_soln::Vector{HierarchicalRouting.TenderSolution}
)
    colormap = distinguishable_colors(length(tender_soln) + 1)[2:end]

    for t_soln in tender_soln
        color = colormap[t_soln.id]

        for route in t_soln.sorties
            nodes = [t_soln.start]
            append!(nodes, [node for node in route.nodes])
            append!(nodes, [t_soln.finish])

            node_lons, node_lats = [wp[1] for wp in nodes], [wp[2] for wp in nodes]

            lines!(ax, node_lons, node_lats, color = color, linewidth = 1)
        end
    end
    return ax
end


end
