module Plot

import ..HierarchicalRouting

using DataFrames
using Rasters
using Statistics

using GeometryBasics

using GLMakie, GeoMakie

"""
    clusters(
        clusters::Union{Vector{HierarchicalRouting.Cluster}, Nothing} = nothing,
        cluster_sequence::Union{DataFrame, Nothing} = nothing,
        cluster_radius::Real = 0,
        centers = false,
        labels = false
    )

Create a plot of nodes by cluster.

# Arguments
- `clusters::Union{Vector{HierarchicalRouting.Cluster}, Nothing}`: Clusters.
- `cluster_sequence::Union{DataFrame, Nothing}`: Cluster by sequence visited.
- `cluster_radius::Real`: Radius of circle to represent clusters.
- `nodes::Bool`: Plot nodes flag.
- `centers::Bool`: Plot cluster centers flag.
- `labels::Bool`: Plot cluster labels flag.

# Returns
- `fig, ax`: Figure and Axis objects.
"""
function clusters(
    ;
    clusters::Union{Vector{HierarchicalRouting.Cluster}, Nothing} = nothing,
    cluster_sequence::Union{DataFrame, Nothing} = nothing,
    cluster_radius::Real = 0,
    nodes::Bool = true,
    centers::Bool = false,
    labels::Bool = false
)
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], xlabel = "Longitude", ylabel = "Latitude")

    clusters!(
        ax,
        clusters = clusters,
        cluster_sequence = cluster_sequence,
        cluster_radius = cluster_radius,
        nodes = nodes,
        centers = centers,
        labels = labels
    )

    return fig, ax
end
"""
    clusters!(
        ax::Axis;
        clusters::Union{Vector{HierarchicalRouting.Cluster}, Nothing} = nothing,
        cluster_sequence::Union{DataFrame, Nothing} = nothing,
        cluster_radius::Real = 0,
        nodes::Bool = true,
        centers::Bool = false,
        labels::Bool = false
    )

Plot nodes by cluster.

# Arguments
- `ax::Axis`: Axis object.
- `clusters::Union{Vector{HierarchicalRouting.Cluster}, Nothing}`: Clusters.
- `cluster_sequence::Union{DataFrame, Nothing}`: Cluster by sequence visited.
- `cluster_radius::Real`: Radius of circle to represent clusters.
- `centers::Bool`: Plot cluster centers flag.
- `labels::Bool`: Plot cluster labels flag.

# Returns
- `ax`: Axis object.
"""
function clusters!(
    ax::Axis;
    clusters::Union{Vector{HierarchicalRouting.Cluster}, Nothing} = nothing,
    cluster_sequence::Union{DataFrame, Nothing} = nothing,
    cluster_radius::Real = 0,
    nodes::Bool = true,
    centers::Bool = false,
    labels::Bool = false
)
    # Validate inputs
    if isnothing(clusters) && isnothing(cluster_sequence)
        error("At least one of `clusters` or `cluster_sequence` must be provided.")
    end

    sequence_id, colormap, centroids = nothing, nothing, nothing
    if !isnothing(cluster_sequence)
        sequence_id = [row.id for row in eachrow(cluster_sequence)[2:end-1]]
        colormap = distinguishable_colors(length(sequence_id) + 2)[3:end]
        centroids = hcat(cluster_sequence.lon, cluster_sequence.lat)[2:end-1,:]
    elseif !isnothing(clusters)
        sequence_id = 1:length(clusters)
        colormap = distinguishable_colors(length(clusters) + 2)[3:end]
        centroids = hcat([cluster.centroid[1] for cluster in clusters], [cluster.centroid[2] for cluster in clusters])
    end

    circle_offsets = cluster_radius > 0 ? (
        cluster_radius .* cos.(range(0, 2π, length=100)),
        cluster_radius .* sin.(range(0, 2π, length=100))
    ) : nothing

    for (idx, seq) in enumerate(sequence_id)
        color = colormap[seq]

        # plot nodes
        if nodes && !isnothing(clusters) && !isempty(clusters[seq].nodes)
            scatter!(ax, clusters[seq].nodes, color = color, markersize = 10, marker = :x)
        end

        center_lon, center_lat = !isnothing(cluster_sequence) ?
            centroids[idx, :] :
            centroids[seq, :]

        if cluster_radius > 0
            circle_lons = center_lon .+ circle_offsets[1]
            circle_lats = center_lat .+ circle_offsets[2]

            poly!(ax, hcat(circle_lons, circle_lats), color = (color, 0.2), strokecolor = color, label = "Cluster Centroids")
        end

        if centers
            scatter!(ax, [center_lon], [center_lat], markersize = 10, color = (color, 0.2), strokewidth = 0)
        end
        if labels
            text!(ax, center_lon, center_lat, text = string(seq), align = (:center, :center), color = color)
        end
    end
    return ax
end

"""
    exclusions(
        exclusions::DataFrame;
        labels::Bool = false
    )

Create a plot of exclusion zones.

# Arguments
- `exclusions::DataFrame`: Exclusion zone polygons.
- `labels::Bool`: Plot exclusion zones flag.

# Returns
- `fig, ax`: Figure and Axis objects.
"""
function exclusions(
    exclusions::DataFrame;
    labels::Bool = false
)
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], xlabel = "Longitude", ylabel = "Latitude")

    exclusions!(ax, exclusions, labels = labels)

    return fig, ax
end
"""
    exclusions!(
        ax::Axis,
        exclusions::DataFrame;
        labels::Bool = false
    )

Plot exclusion zones.

# Arguments
- `ax::Axis`: Axis object.
- `exclusions::DataFrame`: Exclusion zone polygons.
- `labels::Bool`: Plot exclusion zones flag.

# Returns
- `ax`: Axis object.
"""
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
                text!(ax, centroid_x, centroid_y, text = string(i), align = (:center, :center), color = :grey)
            end
        end
    end
    return ax
end

"""
    linestrings(
        route::HierarchicalRouting.Route;
        markers::Bool = false,
        labels::Bool = false,
        color = nothing
    )

Create a plot of LineStrings for mothership route.

# Arguments
- `route::HierarchicalRouting.Route`: Route including nodes and LineStrings.
- `markers::Bool`: Plot waypoints flag.
- `labels::Bool`: Plot LineString labels flag.
- `color`: LineString color.

# Returns
- `fig, ax`: Figure and Axis objects.
"""
function linestrings(
    route::HierarchicalRouting.Route;
    markers::Bool = false,
    labels::Bool = false,
    color = nothing
)
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], xlabel = "Longitude", ylabel = "Latitude")

    linestrings!(ax, route, markers = markers, labels = labels, color = color)

    return fig, ax
end
"""
    linestrings!(
        ax::Axis,
        route::HierarchicalRouting.Route;
        markers::Bool = false,
        labels::Bool = false,
        color = nothing
    )

Plot LineStrings for mothership route.

# Arguments
- `ax::Axis`: Axis object.
- `route::HierarchicalRouting.Route`: Route including nodes and LineStrings.
- `markers::Bool`: Plot waypoints flag.
- `labels::Bool`: Plot LineString labels flag.
- `color`: LineString color.

# Returns
- `ax`: Axis object.
"""
function linestrings!(
    ax::Axis,
    route::HierarchicalRouting.Route;
    markers::Bool = false,
    labels::Bool = false,
    color = nothing
)
    line_strings = route.line_strings
    waypoints = route.nodes[1:end-1]

    n_graphs = length(line_strings)
    color_palette = isnothing(color) ? cgrad(:rainbow, n_graphs) : nothing
    color = isnothing(color) ? color_palette : fill(color, n_graphs)

    waypoint_matrix = hcat([wp[1] for wp in waypoints], [wp[2] for wp in waypoints])

    # Mark waypoints with 'x'
    if markers
        scatter!(waypoint_matrix, marker = 'x', markersize = 10, color = :black)#, label = "Waypoints")
        # series(waypoint_matrix, marker = 'x', markersize = 10, color = :black, label = "Waypoints")
    end

    # Plot LineStrings
    for (idx, line_string) in enumerate(line_strings)
        line_color = color[idx]

        # If line_string is a single LineString, iterate over its points directly.
        points = hasproperty(line_string, :points) ?
                 [Point(p[1], p[2]) for p in line_string.points] :
                 [Point(p[1], p[2]) for l in line_string for p in l.points]
        line_width = line_color == :black ? 3 : 2
        lines!(ax, points, color = line_color, linewidth = line_width)
    end

    if labels
        # Annotate waypoints by sequence
        text!(ax, waypoint_matrix[:,1], waypoint_matrix[:,2] .+ 0.002, text = string.(0:size(waypoint_matrix,1)-1), align = (:center, :center), color = :black)
    end
    return ax
end

"""
    tenders(
        tender_soln::Vector{HierarchicalRouting.TenderSolution}
    )

Create a plot of tender routes within each cluster.

# Arguments
- `tender_soln::Vector{HierarchicalRouting.TenderSolution}`: Tender solutions.

# Returns
- `fig, ax`: Figure and Axis objects.
"""
function tenders(
    tender_soln::Vector{HierarchicalRouting.TenderSolution}
)
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], xlabel = "Longitude", ylabel = "Latitude")

    tenders!(ax, tender_soln)

    return fig, ax
end
"""
    tenders!(
        ax::Axis,
        tender_soln::Vector{HierarchicalRouting.TenderSolution}
    )
    function tenders!(
        ax::Axis,
        tender_soln::Vector{HierarchicalRouting.TenderSolution},
        num_clusters::Int64
    )

Plot tender routes within each cluster, colored by cluster, sequentially shaded by sortie.

# Arguments
- `ax`: Axis object.
- `tender_soln`: Tender solutions.
- `num_clusters`: Number of clusters to color/plot.

# Returns
- `ax`: Axis object.
"""
function tenders!(
    ax::Axis,
    tender_soln::Vector{HierarchicalRouting.TenderSolution}
)
    colormap = distinguishable_colors(length(tender_soln) + 2)[3:end]

    # TODO: Plot critical path (longest) thicker than other paths
    for t_soln in tender_soln
        base_hue = convert_rgb_to_hue(colormap[t_soln.id])
        s = length(t_soln.sorties)
        palette = sequential_palette(base_hue, s+3)[3:end]

        for (sortie, color) in zip(t_soln.sorties, palette[1:s])
            linestrings!(ax, sortie, color = color)
        end
    end
    return ax
end
function tenders!(
    ax::Axis,
    tender_soln::Vector{HierarchicalRouting.TenderSolution},
    num_clusters::Int64
)
    colormap = distinguishable_colors(num_clusters + 2)[3:end]

    # TODO: Plot critical path (longest) thicker than other paths
    for t_soln in tender_soln
        base_hue = convert_rgb_to_hue(colormap[t_soln.id])
        s = length(t_soln.sorties)
        palette = sequential_palette(base_hue, s+3)[3:end]

        for (sortie, color) in zip(t_soln.sorties, palette[1:s])
            linestrings!(ax, sortie, color = color)
        end
    end
    return ax
end

function convert_rgb_to_hue(base_color::RGB{Colors.FixedPointNumbers.N0f8})
    base_color_float64 = ColorTypes.RGB{Float64}(
        float(base_color.r),
        float(base_color.g),
        float(base_color.b)
    )
    return ColorTypes.HSV(base_color_float64).h
end

export clusters, clusters!, exclusions, exclusions!, linestrings, linestrings!, tenders, tenders!
end
