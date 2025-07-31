module Plot

using ..HierarchicalRouting:
    Problem,
    Cluster,
    MothershipSolution,
    TenderSolution,
    MSTSolution,
    Route,
    generate_letter_id,
    critical_path,
    tender_clust_dist,
    mothership_dist_within_clusts

using DataFrames
using Rasters
using Statistics

using GeometryBasics
using GLMakie, GeoMakie

"""
    clusters(
        clusters::Vector{Cluster};
        cluster_radius::Union{Float64, Int64}=0,
        nodes::Bool=true,
        centers::Bool=false,
        labels::Bool=false
    )::Tuple{Figure,Axis}
    clusters(
        cluster_sequence::DataFrame;
        cluster_radius::Union{Float64, Int64}=0,
        centers::Bool=false,
        labels::Bool=false
    )::Tuple{Figure,Axis}

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
    clusters::Vector{Cluster};
    cluster_radius::Union{Float64,Int64}=0,
    nodes::Bool=true,
    centers::Bool=false,
    labels::Bool=false
)::Tuple{Figure,Axis}
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], xlabel="Longitude", ylabel="Latitude")

    clusters!(
        ax,
        clusters,
        cluster_radius,
        nodes,
        centers,
        labels
    )

    return fig, ax
end
function clusters(
    cluster_sequence::DataFrame;
    cluster_radius::Union{Float64,Int64}=0,
    centers::Bool=false,
    labels::Bool=false
)::Tuple{Figure,Axis}
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], xlabel="Longitude", ylabel="Latitude")

    clusters!(
        ax,
        cluster_sequence;
        cluster_radius,
        centers,
        labels
    )

    return fig, ax
end

function _calc_radius_offset(radius::Float64)::Union{Tuple,Nothing}
    offsets = if radius > 0
        (
            radius .* cos.(range(0, 2π, length=100)),
            radius .* sin.(range(0, 2π, length=100))
        )
    else
        nothing
    end

    return offsets
end

"""
    clusters!(
        ax::Axis;
        clusters::Vector{Cluster};
        cluster_radius::Float64=0.0,
        nodes::Bool=true,
        centers::Bool=false,
        labels::Bool=false
    )::Axis
    clusters!(
        ax::Axis;
        cluster_sequence::DataFrame;
        cluster_radius::Float64=0.0,
        centers::Bool=false,
        labels::Bool=false
    )::Axis
    clusters!(
        ax::Axis,
        cluster_radius::Float64,
        sequence_ids::Vector{Int},
        centroids::Vector{Point{2,Float64}},
        centers::Bool=false,
        labels::Bool=false;
        nodes=false,
        clusters=nothing
    )::Axis

Plot nodes by cluster.

# Arguments
- `ax`: Axis object.
- `clusters`: Clusters to plot.
- `cluster_sequence`: Cluster by sequence visited.
- `cluster_radius`: Radius of circle to represent clusters.
- `sequence_ids`: Sequence IDs for clusters.
- `centroids`: Cluster centroids.
- `nodes`: Plot nodes flag.
- `centers`: Plot cluster centers flag.
- `labels`: Plot cluster labels flag.

# Returns
- `ax`: Axis object.
"""
function clusters!(
    ax::Axis,
    clusters::Vector{Cluster};
    cluster_radius::Float64=0.0,
    nodes::Bool=true,
    centers::Bool=false,
    labels::Bool=false
)::Axis
    sequence_ids = getfield.(clusters, :id)
    centroids = getfield.(clusters, :centroid)

    return clusters!(ax, cluster_radius, sequence_ids, centroids, centers, labels; nodes, clusters)
end
function clusters!(
    ax::Axis,
    cluster_sequence::DataFrame;
    cluster_radius::Float64=0.0,
    centers::Bool=false,
    labels::Bool=false
)::Axis
    sequence_ids = cluster_sequence.id[2:end-1]

    centroids = collect(zip(cluster_sequence.lon, cluster_sequence.lat))[2:end-1]
    ordered_centroids = centroids[sortperm(sequence_ids)]
    ordered_ids = sort(sequence_ids)

    return clusters!(ax, cluster_radius, ordered_ids, ordered_centroids, centers, labels)
end
function clusters!(
    ax::Axis,
    cluster_radius::Float64,
    sequence_ids::Vector{Int},
    centroids::Vector{Point{2,Float64}},
    centers::Bool=false,
    labels::Bool=false;
    nodes=false,
    clusters=nothing
)::Axis
    colormap = create_colormap(sequence_ids)
    circle_offsets = _calc_radius_offset(cluster_radius)

    for seq in sequence_ids
        color = colormap[seq]

        # Plot nodes
        if !isnothing(clusters) && nodes
            scatter!(ax, clusters[seq].nodes, color=color, markersize=10, marker=:x)
        end

        center_lon, center_lat = centroids[seq]
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
                text=generate_letter_id(seq - 1),
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
    route(
        route::Route;
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
function route(
    route::Route;
    markers::Bool=false,
    labels::Bool=false,
    color=nothing
)::Tuple{Figure,Axis}
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], xlabel="Longitude", ylabel="Latitude")

    route!(ax, route, markers=markers, labels=labels, color=color)

    return fig, ax
end
"""
    route!(
        ax::Axis,
        route::Route;
        markers::Bool=false,
        labels::Bool=false,
        color=nothing
    )::Axis
    route!(
        ax::Axis,
        ms::MothershipSolution;
        markers::Bool=false,
        labels::Bool=false,
        color=nothing
    )::Axis
    route!(ax::Axis, tender_soln::Vector{TenderSolution})::Axis

Plot LineStrings for mothership route.

# Arguments
- `ax`: Axis object.
- `route`: Route including nodes and LineStrings.
- `ms`: Solution instance containing mothership route.
- `tender_soln`: Solution instance containing tender routes.
- `markers`: Plot waypoints flag.
- `labels`: Plot LineString labels flag.
- `color`: LineString color.

# Returns
- `ax`: Axis object.
"""
function route!(
    ax::Axis,
    route::Route;
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
    ms::MothershipSolution;
    markers::Bool=false,
    labels::Bool=false,
    color=nothing
)::Axis
    return route!(ax, ms.route; markers, labels, color)
end
function route!(ax::Axis, tender_soln::Vector{TenderSolution})::Axis
    return tenders!(ax, tender_soln)
end
function route!(
    ax::Axis,
    tender_soln::Vector{TenderSolution},
    colormap::Vector{RGB{Colors.FixedPointNumbers.N0f8}}
)::Axis
    for t_soln in tender_soln
        base_hue = convert_rgb_to_hue(colormap[t_soln.id])
        s = length(t_soln.sorties)
        palette = sequential_palette(base_hue, s + 3)[3:end]

        for (sortie, color) in zip(t_soln.sorties, palette[1:s])
            route!(ax, sortie, color=color)
        end
    end
    return ax
end

"""
    tenders(
        tender_soln::Vector{TenderSolution}
    )::Tuple{Figure, Axis}

Create a plot of tender routes within each cluster.

# Arguments
- `tender_soln::Vector{TenderSolution}`: Tender solutions.

# Returns
- `fig, ax`: Figure and Axis objects.
"""
function tenders(
    tender_soln::Vector{TenderSolution}
)::Tuple{Figure,Axis}
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], xlabel="Longitude", ylabel="Latitude")

    tenders!(ax, tender_soln)

    return fig, ax
end

"""
    tenders!(
        ax::Axis,
        tender_soln::Vector{TenderSolution}
    )::Axis
    function tenders!(
        ax::Axis,
        tender_soln::Vector{TenderSolution},
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
    tender_soln::Vector{TenderSolution}
)::Axis
    colormap = create_colormap(getfield.(tender_soln, :id))
    return route!(ax, tender_soln, colormap)
end
function tenders!(
    ax::Axis,
    tender_soln::Vector{TenderSolution},
    num_clusters::Int64
)::Axis
    colormap = create_colormap(1:num_clusters)
    return route!(ax, tender_soln, colormap)
end

"""
    solution(
        problem::Problem,
        soln::MSTSolution;
        cluster_radius::Float64=0.0,
        show_mothership_exclusions::Bool=true,
        show_tenders_exclusions::Bool=true,
        show_mothership::Bool=true,
        show_tenders::Bool=true,
    )::Figure
    solution(
        problem::Problem,
        soln_a::MSTSolution,
        soln_b::MSTSolution;
        cluster_radius::Float64=0.0,
        show_mothership_exclusions::Bool=true,
        show_tenders_exclusions::Bool=true,
        show_mothership::Bool=true,
        show_tenders::Bool=true,
    )::Figure

Create a plot of the full routing solution, including:
- exclusion zones for the **mothership** and **tenders**,
- mothership route,
- tender sorties (coloured by cluster), and
- clustered target points (coloured by cluster).

# Arguments
- `problem`: The hierarchical routing problem instance.
- `soln`: The full solution to the problem.
- `soln_a`: The first solution to compare.
- `soln_b`: The second solution to compare.
- `cluster_radius`: Radius of the cluster circles to display around cluster centroids.
- `show_mothership_exclusions`: Whether to show **mothership** exclusion zones.
- `show_tenders_exclusions`: Whether to show **tender** exclusion zones.
- `show_mothership`: Whether to show the **mothership** route.
- `show_tenders`: Whether to show **tender** routes.

# Returns
- `fig`: The created Figure object containing the plot.
"""
function solution(
    problem::Problem,
    soln::MSTSolution;
    cluster_radius::Float64=0.0,
    show_mothership_exclusions::Bool=true,
    show_tenders_exclusions::Bool=true,
    show_mothership::Bool=true,
    show_tenders::Bool=true,
)::Figure
    fig = Figure(size=(750, 880))
    ax = Axis(fig[1, 1], xlabel="Longitude", ylabel="Latitude")

    # Exclusions
    show_mothership_exclusions && exclusions!(ax, problem.mothership.exclusion; labels=false)
    show_tenders_exclusions && exclusions!(ax, problem.tenders.exclusion; labels=false)

    # Tender sorties/routes
    show_tenders && route!(ax, soln.tenders[end])

    # Mothership route
    if show_mothership
        route!(ax, soln.mothership_routes[end]; markers=true, labels=true, color=:black)
    end

    # Clusters
    clusters!(
        ax,
        soln.cluster_sets[end];
        labels=true,
        cluster_radius
    )

    # Annotate critical path cost
    vessel_weightings = (problem.mothership.weighting, problem.tenders.weighting)
    critical_path_dist = critical_path(soln, vessel_weightings)
    total_dist = critical_distance_path(soln, vessel_weightings)

    annotate_cost!(
        ax,
        critical_path_dist;
        position=(0.95, 0.07),
        fontsize=14,
        color=:black
    )
    annotate_cost!(
        ax,
        total_dist;
        position=(0.95, 0.01),
        fontsize=14,
        color=:black,
        metric="critical_distance_path()\ntotal dist"
    )

    highlight_critical_path!(ax, soln, vessel_weightings)

    return fig
end
function solution(
    problem::Problem,
    soln_a::MSTSolution,
    soln_b::MSTSolution;
    cluster_radius::Float64=0.0,
    show_mothership_exclusions::Bool=true,
    show_tenders_exclusions::Bool=true,
    show_mothership::Bool=true,
    show_tenders::Bool=true,
)::Figure
    fig = Figure(size=(1350, 750))  ## 2 fig plot
    ax1, ax2 =
        Axis(fig[1, 1], xlabel="Longitude", ylabel="Latitude"),
        Axis(fig[1, 2], xlabel="Longitude")

    # Exclusions
    if show_mothership_exclusions
        exclusions!.([ax1, ax2], Ref(problem.mothership.exclusion); labels=false)
    end
    if show_tenders_exclusions
        exclusions!.([ax1, ax2], Ref(problem.tenders.exclusion); labels=false)
    end

    # Clusters
    clusters!.(
        [ax1, ax2],
        [soln_a.cluster_sets[end], soln_b.cluster_sets[end]];
        cluster_radius=cluster_radius,
        nodes=true,
        centers=false,
        labels=true,
    )

    # Tender sorties/routes
    if show_tenders
        tenders!.(
            [ax1, ax2],
            [soln_a.tenders[end], soln_b.tenders[end]]
        )
    end

    # Mothership route
    if show_mothership
        route!.(
            [ax1, ax2],
            [soln_a.mothership_routes[end], soln_b.mothership_routes[end]];
            markers=true,
            labels=true,
            color=:black
        )
    end

    # Annotate critical path costs
    vessel_weightings = (problem.mothership.weighting, problem.tenders.weighting)
    critical_path_dist_a = critical_path(soln_a, vessel_weightings)
    critical_path_dist_b = critical_path(soln_b, vessel_weightings)
    total_dist_a = critical_distance_path(soln_a, vessel_weightings)
    total_dist_b = critical_distance_path(soln_b, vessel_weightings)

    annotate_cost!.(
        [ax1, ax2],
        [critical_path_dist_a, critical_path_dist_b];
        position=(0.95, 0.07),
        fontsize=14,
        color=:black
    )

    annotate_cost!.(
        [ax1, ax2],
        [total_dist_a, total_dist_b];
        position=(0.95, 0.01),
        fontsize=14,
        color=:black,
        metric="critical_distance_path()\ntotal dist"
    )

    return fig
end

"""
    function solution_disturbances(
        problem::Problem,
        solution_disturbed::MSTSolution,
        disturbance_clusters::Set{Int64};
        cluster_radius::Float64=0.0,
        show_mothership_exclusions::Bool=true,
        show_tenders_exclusions::Bool=true,
        show_mothership::Bool=true,
        show_tenders::Bool=true,
    )::Figure

Create a plot of the solution at each progressive disturbance event, including:
- exclusion zones for the **mothership** and **tenders**,
- mothership route,
- tender sorties (coloured by cluster), and
- clustered points remaining remaining after disturbance event (coloured by cluster).

#! NOTE: This function assumes 2 disturbance events.
#TODO: Generalize to any number of disturbance events.

# Arguments
- `problem`: The hierarchical routing problem instance.
- `solution_disturbed`: The solution with disturbances.
- `disturbance_clusters`: Set of disturbance cluster IDs.
- `cluster_radius`: Radius of the cluster circles to display around cluster centroids.
- `show_mothership_exclusions`: Whether to show **mothership** exclusion zones.
- `show_tenders_exclusions`: Whether to show **tender** exclusion zones.
- `show_mothership`: Whether to show the **mothership** route.
- `show_tenders`: Whether to show **tender** routes.

# Returns
- `fig`: The created Figure object containing the plot.
"""
function solution_disturbances(
    problem::Problem,
    solution_disturbed::MSTSolution,
    disturbance_clusters::Set{Int64};
    cluster_radius::Float64=0.0,
    show_mothership_exclusions::Bool=false,
    show_tenders_exclusions::Bool=true,
    show_mothership::Bool=true,
    show_tenders::Bool=true,
)::Figure
    #! NOTE: This function assumes 2 disturbance events.
    #TODO: Generalize to any number of disturbance events.
    fig = Figure(size=(1650, 600))  ## 3 fig plot
    ax1, ax2, ax3 =
        Axis(fig[1, 1], xlabel="Longitude", ylabel="Latitude"),
        Axis(fig[1, 2], xlabel="Longitude"),
        Axis(fig[1, 3], xlabel="Longitude")

    # Exclusions
    if show_mothership_exclusions
        exclusions!.([ax1, ax2, ax3], Ref(problem.mothership.exclusion); labels=false)
    end
    if show_tenders_exclusions
        exclusions!.([ax1, ax2, ax3], Ref(problem.tenders.exclusion); labels=false)
    end

    # Clusters
    clusters!.(
        [ax1, ax2, ax3],
        [solution_disturbed.cluster_sets[end], solution_disturbed.cluster_sets[end], solution_disturbed.cluster_sets[end]];
        cluster_radius=cluster_radius,
        nodes=true,
        centers=false,
        labels=true,
    )

    # Tender sorties/routes
    if show_tenders
        ordered_disturbances = sort(unique(disturbance_clusters))
        route!.(
            [ax1, ax2, ax3],
            [
                solution_disturbed.tenders[1][1:ordered_disturbances[1]-1],
                solution_disturbed.tenders[2][1:ordered_disturbances[2]-1],
                solution_disturbed.tenders[3]
            ],
            length.(solution_disturbed.tenders[1:3])
        )
    end

    # Mothership route
    if show_mothership
        route!.(
            [ax1, ax2, ax3],
            [
                solution_disturbed.mothership_routes[1].route,
                solution_disturbed.mothership_routes[2].route,
                solution_disturbed.mothership_routes[3].route
            ];
            labels=true,
            color=:black
        )
    end

    return fig
end

function annotate_cost!(
    ax::Axis,
    cost::Float64;
    position::Tuple{Float64,Float64}=(0.95, 0.02),
    fontsize::Int=14,
    color::Symbol=:black,
    metric::String="Critical path"
)::Axis
    # Annotate the cost of the critical path on the plot
    cost_km = cost / 1000  # Convert cost to km
    text!(
        ax,
        position...,
        text="$metric: $(round(cost_km, digits=2)) km",
        align=(:right, :bottom),
        space=:relative,
        fontsize=fontsize,
        color=color
    )
    return ax
end

function highlight_critical_path!(
    ax::Axis,
    soln::MSTSolution,
    vessel_weightings::NTuple{2,AbstractFloat}=(1.0, 1.0);
    color=:red,
    linewidth=4
)
    # Unpack
    tenders = soln.tenders[end]
    ms_route = soln.mothership_routes[end].route
    num_clusters = length(tenders)

    # Within clusters
    clust_sorties = tender_clust_dist.(tenders)
    clust_sorties = map(x -> isempty(x) ? [0.0] : x, clust_sorties)
    longest_sortie_cost = maximum.(clust_sorties) .* vessel_weightings[2]
    ms_within = mothership_dist_within_clusts(ms_route)[1:num_clusters]
    mothership_sub_cost = vessel_weightings[1] .* ms_within

    # For each cluster, identify critical path
    for j in 1:num_clusters
        if longest_sortie_cost[j] ≥ mothership_sub_cost[j]
            # Draw the longest tender sortie
            routes = tenders[j].sorties
            idx = argmax(tender_clust_dist(tenders[j]))
            route!(ax, routes[idx]; color=color)
        else
            # Draw the mothership route
            start_point = ms_route.nodes[2j]
            end_point = ms_route.nodes[2j+1]
            start_segment = findfirst(
                ==(start_point),
                getindex.(getfield.(ms_route.line_strings, :points), 1)
            )
            end_segment = findfirst(
                ==(end_point),
                getindex.(getfield.(ms_route.line_strings, :points), 2)
            )

            lines!.(
                Ref(ax),
                ms_route.line_strings[start_segment:end_segment];
                color=color,
                linewidth=linewidth
            )
        end
    end

    # Highlight mothership segments between clusters
    start_segment_points::Vector{Point{2,Float64}} = ms_route.nodes[1:2:end-1]
    end_segment_points::Vector{Point{2,Float64}} = ms_route.nodes[2:2:end]

    start_segment = findfirst.(
        .==(start_segment_points),
        Ref(getindex.(getfield.(ms_route.line_strings, :points), 1))
    )
    end_segment = findfirst.(
        .==(end_segment_points),
        Ref(getindex.(getfield.(ms_route.line_strings, :points), 2))
    )
    segments = [ms_route.line_strings[i:j] for (i, j) in zip(start_segment, end_segment)]

    lines!.(
        Ref(ax),
        segments;
        color=color,
        linewidth=linewidth
    )

    return nothing
end

function create_colormap(ids::Vector{Int})::Vector{RGB{Colors.FixedPointNumbers.N0f8}}
    # Create custom colormap, skipping the first two colors (yellow and black)
    max_id = maximum(ids)
    colormap = distinguishable_colors(max_id + 2)[3:end]
    return colormap
end
function create_colormap(ids::UnitRange{Int64})
    return create_colormap(collect(ids))
end

function convert_rgb_to_hue(base_color::RGB{Colors.FixedPointNumbers.N0f8})
    base_color_float64 = ColorTypes.RGB{Float64}(
        float(base_color.r),
        float(base_color.g),
        float(base_color.b)
    )
    return ColorTypes.HSV(base_color_float64).h
end

export clusters, clusters!, exclusions, exclusions!, route, route!, tenders, tenders!
end
