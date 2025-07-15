module Plot

import ..HierarchicalRouting

using DataFrames
using Rasters
using Statistics

using GeometryBasics

using GLMakie, GeoMakie

"""
    clusters(
        ;
        clusters::Union{Vector{HierarchicalRouting.Cluster},Nothing}=nothing,
        cluster_sequence::Union{DataFrame,Nothing}=nothing,
        cluster_radius::Real=0,
        nodes::Bool=true,
        centers::Bool=false,
        labels::Bool=false
    )::Tuple{Figure, Axis}

Create a plot of nodes by cluster.

# Arguments
- `clusters`: Clusters to plot.
- `cluster_sequence`: Cluster by sequence visited.
- `cluster_radius`: Radius of circle to represent clusters.
- `nodes`: Plot nodes flag.
- `centers`: Plot cluster centers flag.
- `labels`: Plot cluster labels flag.

# Returns
- `fig, ax`: Figure and Axis objects.
"""
function clusters(
    ;
    clusters::Union{Vector{HierarchicalRouting.Cluster},Nothing}=nothing,
    cluster_sequence::Union{DataFrame,Nothing}=nothing,
    cluster_radius::Real=0,
    nodes::Bool=true,
    centers::Bool=false,
    labels::Bool=false
)::Tuple{Figure,Axis}
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], xlabel="Longitude", ylabel="Latitude")

    clusters!(
        ax,
        clusters=clusters,
        cluster_sequence=cluster_sequence,
        cluster_radius=cluster_radius,
        nodes=nodes,
        centers=centers,
        labels=labels
    )

    return fig, ax
end
"""
    clusters!(
        ax::Axis;
        clusters::Union{Vector{HierarchicalRouting.Cluster},Nothing}=nothing,
        cluster_sequence::Union{DataFrame,Nothing}=nothing,
        cluster_radius::Real=0,
        nodes::Bool=true,
        centers::Bool=false,
        labels::Bool=false
    )::Axis

Plot nodes by cluster.

# Arguments
- `ax`: Axis object.
- `clusters`: Clusters to plot.
- `cluster_sequence`: Cluster by sequence visited.
- `cluster_radius`: Radius of circle to represent clusters.
- `nodes`: Plot nodes flag.
- `centers`: Plot cluster centers flag.
- `labels`: Plot cluster labels flag.

# Returns
- `ax`: Axis object.
"""
function clusters!(
    ax::Axis;
    clusters::Union{Vector{HierarchicalRouting.Cluster},Nothing}=nothing,
    cluster_sequence::Union{DataFrame,Nothing}=nothing,
    cluster_radius::Real=0,
    nodes::Bool=true,
    centers::Bool=false,
    labels::Bool=false
)::Axis
    # Validate inputs
    if isnothing(clusters) && isnothing(cluster_sequence)
        error("At least one of `clusters` or `cluster_sequence` must be provided.")
    end

    sequence_id, colormap, centroids = nothing, nothing, nothing
    if !isnothing(cluster_sequence)
        sequence_id = [row.id for row in eachrow(cluster_sequence)[2:end-1]]
        colormap = distinguishable_colors(length(sequence_id) + 2)[3:end]
        centroids = hcat(cluster_sequence.lon, cluster_sequence.lat)[2:end-1, :]
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
            scatter!(ax, clusters[seq].nodes, color=color, markersize=10, marker=:x)
        end

        center_lon, center_lat = !isnothing(cluster_sequence) ?
                                 centroids[idx, :] :
                                 centroids[seq, :]

        if cluster_radius > 0
            circle_lons = center_lon .+ circle_offsets[1]
            circle_lats = center_lat .+ circle_offsets[2]

            poly!(ax, hcat(circle_lons, circle_lats), color=(color, 0.2), strokecolor=color, label="Cluster Centroids")
        end

        if centers
            scatter!(ax, [center_lon], [center_lat], markersize=10, color=(color, 0.2), strokewidth=0)
        end
        if labels
            text!(
                ax,
                center_lon,
                center_lat,
                text=HierarchicalRouting.generate_letter_id(seq - 1),
                font="bold",
                fontsize=24,
                align=(:center, :center),
                color=color
            )
        end
    end
    return ax
end

"""
    exclusions(
        exclusions::DataFrame;
        labels::Bool=false
    )::Tuple{Figure, Axis}

Create a plot of exclusion zones.

# Arguments
- `exclusions`: Exclusion zone polygons.
- `labels`: Plot exclusion zones flag.

# Returns
- `fig, ax`: Figure and Axis objects.
"""
function exclusions(
    exclusions::DataFrame;
    labels::Bool=false
)::Tuple{Figure,Axis}
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], xlabel="Longitude", ylabel="Latitude")

    exclusions!(ax, exclusions, labels=labels)

    return fig, ax
end
"""
    exclusions!(
        ax::Axis,
        exclusions::DataFrame;
        labels::Bool=false
    )::Axis

Plot exclusion zones.

# Arguments
- `ax`: Axis object.
- `exclusions`: Exclusion zone polygons.
- `labels`: Plot exclusion zones flag.

# Returns
- `ax`: Axis object.
"""
function exclusions!(
    ax::Axis,
    exclusions::DataFrame;
    labels::Bool=false
)::Axis
    for (i, zone) in enumerate(eachrow(exclusions))
        polygon = zone[:geometry]
        for ring in GeoInterface.coordinates(polygon)
            xs, ys = [coord[1] for coord in ring], [coord[2] for coord in ring]
            poly!(ax, xs, ys, color=(:gray, 0.5), strokecolor=:black)#, label = "Exclusion Zone")

            if labels
                centroid_x, centroid_y = mean(xs), mean(ys)
                text!(ax, centroid_x, centroid_y, text=string(i), align=(:center, :center), color=:grey)
            end
        end
    end
    return ax
end

"""
    linestrings(
        route::HierarchicalRouting.Route;
        markers::Bool=false,
        labels::Bool=false,
        color=nothing
    )::Tuple{Figure, Axis}

Create a plot of LineStrings for mothership route.

# Arguments
- `route`: Route including nodes and LineStrings.
- `markers`: Plot waypoints flag.
- `labels`: Plot LineString labels flag.
- `color`: LineString color.

# Returns
- `fig, ax`: Figure and Axis objects.
"""
function linestrings(
    route::HierarchicalRouting.Route;
    markers::Bool=false,
    labels::Bool=false,
    color=nothing
)::Tuple{Figure,Axis}
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], xlabel="Longitude", ylabel="Latitude")

    linestrings!(ax, route, markers=markers, labels=labels, color=color)

    return fig, ax
end
"""
    linestrings!(
        ax::Axis,
        route::HierarchicalRouting.Route;
        markers::Bool=false,
        labels::Bool=false,
        color=nothing
    )::Axis

Plot LineStrings for mothership route.

# Arguments
- `ax`: Axis object.
- `route`: Route including nodes and LineStrings.
- `markers`: Plot waypoints flag.
- `labels`: Plot LineString labels flag.
- `color`: LineString color.

# Returns
- `ax`: Axis object.
"""
function linestrings!(
    ax::Axis,
    route::HierarchicalRouting.Route;
    markers::Bool=false,
    labels::Bool=false,
    color=nothing
)::Axis
    line_strings = route.line_strings
    waypoints = route.nodes[1:end-1]

    n_graphs = length(line_strings)
    color_palette = isnothing(color) ? cgrad(:rainbow, n_graphs) : nothing
    color = isnothing(color) ? color_palette : fill(color, n_graphs)

    waypoint_matrix = hcat([wp[1] for wp in waypoints], [wp[2] for wp in waypoints])

    # Mark waypoints with 'x'
    if markers
        scatter!(waypoint_matrix, marker='x', markersize=10, color=:black)#, label = "Waypoints")
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
        lines!(ax, points, color=line_color, linewidth=line_width)
    end

    if labels
        # Annotate waypoints by sequence
        text!(ax, waypoint_matrix[:, 1], waypoint_matrix[:, 2] .+ 0.002, text=string.(0:size(waypoint_matrix, 1)-1), align=(:center, :center), color=:black)
    end
    return ax
end

function route!(
    ax::Axis,
    ms::HierarchicalRouting.MothershipSolution;
    markers::Bool=false,
    labels::Bool=false,
    color=nothing
)
    return linestrings!(ax, ms.route; markers, labels, color)
end

function route!(ax::Axis, tender_soln::Vector{HierarchicalRouting.TenderSolution})
    return tenders!(ax, tender_soln)
end

"""
    tenders(
        tender_soln::Vector{HierarchicalRouting.TenderSolution}
    )::Tuple{Figure, Axis}

Create a plot of tender routes within each cluster.

# Arguments
- `tender_soln::Vector{HierarchicalRouting.TenderSolution}`: Tender solutions.

# Returns
- `fig, ax`: Figure and Axis objects.
"""
function tenders(
    tender_soln::Vector{HierarchicalRouting.TenderSolution}
)::Tuple{Figure,Axis}
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], xlabel="Longitude", ylabel="Latitude")

    tenders!(ax, tender_soln)

    return fig, ax
end

"""
    tenders!(
        ax::Axis,
        tender_soln::Vector{HierarchicalRouting.TenderSolution}
    )::Axis
    function tenders!(
        ax::Axis,
        tender_soln::Vector{HierarchicalRouting.TenderSolution},
        num_clusters::Int64
    )::Axis

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
)::Axis
    # Create custom colormap, skipping the first two colors (yellow and black)
    colormap = distinguishable_colors(length(tender_soln) + 2)[3:end]

    # TODO: Plot critical path (longest) thicker than other paths
    for (t_n, t_soln) in enumerate(tender_soln)
        base_hue = convert_rgb_to_hue(colormap[t_n])
        s = length(t_soln.sorties)
        palette = sequential_palette(base_hue, s + 3)[3:end]

        for (sortie, color) in zip(t_soln.sorties, palette[1:s])
            linestrings!(ax, sortie, color=color)
        end
    end
    return ax
end
function tenders!(
    ax::Axis,
    tender_soln::Vector{HierarchicalRouting.TenderSolution},
    num_clusters::Int64
)::Axis
    colormap = distinguishable_colors(num_clusters + 2)[3:end]

    # TODO: Plot critical path (longest) thicker than other paths
    for t_soln in tender_soln
        base_hue = convert_rgb_to_hue(colormap[t_soln.id])
        s = length(t_soln.sorties)
        palette = sequential_palette(base_hue, s + 3)[3:end]

        for (sortie, color) in zip(t_soln.sorties, palette[1:s])
            linestrings!(ax, sortie, color=color)
        end
    end
    return ax
end

"""
    solution(
        problem::HierarchicalRouting.Problem,
        soln::HierarchicalRouting.MSTSolution;
        cluster_radius::Float64=0.0,
        show_mothership_exclusions::Bool=false,
        show_tenders_exclusions::Bool=true,
        show_mothership::Bool=true,
        show_tenders::Bool=true,
        fig_size=(750, 880)
    )::Figure

Create a plot of the full routing solution, including:
- exclusion zones for the **mothership** and **tenders**,
- mothership route,
- tender sorties (coloured by cluster), and
- clustered target points (coloured by cluster).

# Arguments
- `problem`: The hierarchical routing problem instance.
- `soln`: The full solution to the problem.
- `cluster_radius`: Radius of the cluster circles to display around cluster centroids.
- `show_mothership_exclusions`: Whether to show **mothership** exclusion zones.
- `show_tenders_exclusions`: Whether to show **tender** exclusion zones.
- `show_mothership`: Whether to show the **mothership** route.
- `show_tenders`: Whether to show **tender** routes.
- `fig_size`: Size of the figure.

# Returns
- `fig`: The created Figure object containing the plot.
"""
function solution(
    problem::HierarchicalRouting.Problem,
    soln::HierarchicalRouting.MSTSolution;
    cluster_radius::Float64=0.0,
    show_mothership_exclusions::Bool=false,
    show_tenders_exclusions::Bool=true,
    show_mothership::Bool=true,
    show_tenders::Bool=true,
)::Figure
    fig = Figure(size=(750, 880))
    ax = Axis(fig[1, 1], xlabel="Longitude", ylabel="Latitude")

    # Exclusions
    show_mothership_exclusions && exclusions!(ax, problem.mothership.exclusion; labels=false)
    show_tenders_exclusions && exclusions!(ax, problem.tenders.exclusion; labels=false)

    # Clusters
    clusters!(
        ax,
        clusters=soln.cluster_sets[end],
        nodes=true,
        centers=false,
        labels=true,
        cluster_radius=cluster_radius
    )

    # Mothership route
    if show_mothership
        route!(ax, soln.mothership_routes[end]; markers=true, labels=true, color=:black)
    end

    # Tender sorties/routes
    show_tenders && route!(ax, soln.tenders[end])

    return fig
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
