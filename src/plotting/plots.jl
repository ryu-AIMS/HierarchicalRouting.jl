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
    centers = false,
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
        if !isnothing(clusters) && !isempty(clusters[seq].nodes)
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
            text!(ax, center_lon, center_lat, text = string(seq), align = (:center, :center), color = :black)
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
                text!(ax, centroid_x, centroid_y, text = string(i), align = (:center, :center), color = :blue)
            end
        end
    end
    return ax
end

"""
    linestrings(
        line_strings::Vector{LineString{2, Float64}};
        labels::Bool = false
    )

Create a plot of LineStrings for mothership route.

# Arguments
- `line_strings::Vector{LineString{2, Float64}}`: LineStrings for mothership route.
- `labels::Bool`: Plot LineString labels flag.

# Returns
- `fig, ax`: Figure and Axis objects.
"""
function linestrings(
    line_strings::Vector{LineString{2, Float64}};
    labels::Bool = false
)
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], xlabel = "Longitude", ylabel = "Latitude")

    linestrings!(ax, line_strings, labels = labels)

    return fig, ax
end
"""
    linestrings!(
        ax::Axis,
        line_strings::Vector{LineString{2, Float64}};
        labels::Bool = false
    )

Plot LineStrings for mothership route.

# Arguments
- `ax::Axis`: Axis object.
- `line_strings::Vector{LineString{2, Float64}}`: LineStrings for mothership route.
- `labels::Bool`: Plot LineString labels flag.

# Returns
- `ax`: Axis object.
"""
function linestrings!(
    ax::Axis,
    line_strings::Vector{LineString{2, Float64}};
    markers::Bool = false,
    labels::Bool = false,
    color = nothing
)
    n_graphs = length(line_strings)
    color_palette = isnothing(color) ? cgrad(:rainbow, n_graphs) : nothing
    color = isnothing(color) ? color_palette : fill(color, n_graphs)

    waypoints = [line[1][1] for line in line_strings]
    waypoint_matrix = hcat([wp[1] for wp in waypoints], [wp[2] for wp in waypoints])

    # Mark waypoints with 'x'
    if markers
        scatter!(waypoint_matrix, marker = 'x', markersize = 10, color = :black)#, label = "Waypoints")
        # series(waypoint_matrix, marker = 'x', markersize = 10, color = :black, label = "Waypoints")
    end

    # Plot LineStrings
    for (idx, line_string) in enumerate(line_strings)
        line_color = color[idx]
        points = [Point(p[1], p[2]) for l in line_string for p in l.points]
        lines!(ax, points, color = line_color, linewidth = 2)
    end

    if labels
        # Annotate waypoints by sequence
        text!(ax, waypoint_matrix[:,1], waypoint_matrix[:,2] .+ 30, text = string.(0:size(waypoint_matrix,1)-1), align = (:center, :center), color = :black)
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

Plot tender routes within each cluster.

# Arguments
- `ax::Axis`: Axis object.
- `tender_soln::Vector{HierarchicalRouting.TenderSolution}`: Tender solutions.

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
        color = colormap[t_soln.id]
        [linestrings!(ax, s.line_strings, color = color, labels = false) for s in t_soln.sorties]
    end
    return ax
end


end
