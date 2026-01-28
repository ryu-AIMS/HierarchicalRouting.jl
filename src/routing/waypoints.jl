
"""
    adjust_waypoint(
        waypoint::Point{2,Float64},
        exclusions::Vector{IGeometry{wkbPolygon}},
    )::Point{2,Float64}

Adjust waypoint if inside exclusion zone to closest boundary point outside exclusion zone.

# Returns
Adjusted waypoint if inside exclusion zone, else original waypoint.
"""
function adjust_waypoint(
    waypoint::Point{2,Float64},
    exclusions::POLY_VEC,
)::Point{2,Float64}
    waypoint_geom = AG.createpoint(waypoint.data)
    convex_hulls::POLY_VEC = [AG.convexhull(exclusion) for exclusion in exclusions]
    containing_polygons::POLY_VEC = filter(
        hull -> AG.contains(hull, waypoint_geom),
        convex_hulls
    )

    if isempty(containing_polygons)
        return waypoint
    end

    # TODO Consider cases for multiple exclusion zones
    union_poly::IGeometry{wkbPolygon} = reduce(AG.union, containing_polygons)

    exterior_ring::IGeometry{wkbLineString} = AG.getgeom(union_poly, 0)
    n_points::Int32 = AG.ngeom(exterior_ring)

    boundary_points = Vector{Point{2,Float64}}(undef, n_points)
    pt = AG.creategeom(wkbPoint)

    # iterate safely within all vertices on the exterior ring to populate pts
    @inbounds for i in 1:n_points
        pt = AG.getpoint(exterior_ring, i - 1) # 0-based indexing
        boundary_points[i] = Point{2,Float64}(pt[1], pt[2])
    end

    valid_boundary_points = Point{2,Float64}[]
    sizehint!(valid_boundary_points, length(boundary_points))
    @inbounds for p in boundary_points
        # If a point is not in the union (exclusion), it’s valid.
        if !point_in_exclusion(p, convex_hulls)
            push!(valid_boundary_points, p)
        end
    end

    # ignore `sqrt` for performance - gives same result
    closest_point::Point{2,Float64} = argmin(
        p -> ((p[1] - waypoint[1])^2 + (p[2] - waypoint[2])^2),
        valid_boundary_points
    )
    # Recursively call adjust_waypoint() until waypoint is outside exclusion zone
    if point_in_exclusion(closest_point, convex_hulls)
        return adjust_waypoint(closest_point, exclusions)
    end

    return closest_point
end

"""
    get_waypoints(sequence::DataFrame, exclusions::POLY_VEC)::DataFrame
    get_waypoints(
        current_point::Point{2,Float64},
        sequence::DataFrame,
        exclusions::POLY_VEC
    )::DataFrame

Calculate mothership waypoints between sequential clusters.
For each cluster, waypoint 1/3 dist before and after cluster centroid,
unless within exclusion zone, then adjust to closest boundary point.

# Arguments
- `sequence`: `id`, and centroid (`lat`, `long`) coordinates in sequence; including depot as
    the first and last rows with `id=0`.
- `exclusions`: Exclusion zone polygons.

# Returns
- A DataFrame containing waypoints and connecting clusters.
"""
function get_waypoints(sequence::DataFrame, exclusions::POLY_VEC)::DataFrame
    # TODO: Use graph to determine feasible paths and waypoints along path
    n_cluster_seqs = nrow(sequence)

    waypoints = Vector{Point{2,Float64}}(undef, 2 * (n_cluster_seqs - 2) + 2)
    connecting_clusters = Vector{NTuple{2,Int64}}(undef, 2 * (n_cluster_seqs - 2) + 2)

    waypoints[1] = Point{2,Float64}(sequence.lon[1], sequence.lat[1])
    connecting_clusters[1] = (sequence.id[1], sequence.id[1])

    for i in 2:(n_cluster_seqs-1)
        prev_lon, current_lon, next_lon = sequence.lon[i-1:i+1]
        prev_lat, current_lat, next_lat = sequence.lat[i-1:i+1]
        prev_clust, current_clust, next_clust = sequence.id[i-1:i+1]

        # Compute waypoints before and after the current cluster centroid
        prev_waypoint = Point{2,Float64}(
            2 / 3 * current_lon + 1 / 3 * prev_lon,
            2 / 3 * current_lat + 1 / 3 * prev_lat
        )
        next_waypoint = Point{2,Float64}(
            2 / 3 * current_lon + 1 / 3 * next_lon,
            2 / 3 * current_lat + 1 / 3 * next_lat
        )

        # Adjust waypoints if they are inside exclusion polygons
        prev_waypoint = point_in_exclusion(prev_waypoint, exclusions) ?
                        adjust_waypoint(prev_waypoint, exclusions) : prev_waypoint
        next_waypoint = point_in_exclusion(next_waypoint, exclusions) ?
                        adjust_waypoint(next_waypoint, exclusions) : next_waypoint

        waypoints[2*i-2] = prev_waypoint
        connecting_clusters[2*i-2] = (prev_clust, current_clust)

        waypoints[2*i-1] = next_waypoint
        connecting_clusters[2*i-1] = (current_clust, next_clust)
    end

    waypoints[2*(n_cluster_seqs-2)+2] = (sequence.lon[end], sequence.lat[end])
    connecting_clusters[2*(n_cluster_seqs-2)+2] = (sequence.id[end], sequence.id[end])

    return DataFrame(waypoint=waypoints, connecting_clusters=connecting_clusters)
end
function get_waypoints(
    current_point::Point{2,Float64},
    sequence::DataFrame,
    exclusions::POLY_VEC
)::DataFrame
    # TODO: Use graph to determine feasible paths and waypoints along path
    n_cluster_seqs = nrow(sequence)

    waypoints = Vector{Point{2,Float64}}(undef, 2 * (n_cluster_seqs - 2) + 2)
    connecting_clusters = Vector{NTuple{2,Int64}}(undef, 2 * (n_cluster_seqs - 2) + 2)

    waypoints[1] = current_point
    connecting_clusters[1] = (sequence.id[1], sequence.id[2])

    for i in 2:(n_cluster_seqs-1)
        prev_lon, current_lon, next_lon = sequence.lon[i-1:i+1]
        prev_lat, current_lat, next_lat = sequence.lat[i-1:i+1]
        prev_clust, current_clust, next_clust = sequence.id[i-1:i+1]

        # Compute waypoints before and after the current cluster centroid
        prev_waypoint = Point{2,Float64}(
            2 / 3 * current_lon + 1 / 3 * prev_lon,
            2 / 3 * current_lat + 1 / 3 * prev_lat
        )
        next_waypoint = Point{2,Float64}(
            2 / 3 * current_lon + 1 / 3 * next_lon,
            2 / 3 * current_lat + 1 / 3 * next_lat
        )
        #! prev_waypoint in exclusion calls adjust_waypoint() to find closest point outside exclusion
        # Adjust waypoints if they are inside exclusion polygons
        prev_waypoint = point_in_exclusion(prev_waypoint, exclusions) ?
                        adjust_waypoint(prev_waypoint, exclusions) :
                        prev_waypoint
        next_waypoint = point_in_exclusion(next_waypoint, exclusions) ?
                        adjust_waypoint(next_waypoint, exclusions) :
                        next_waypoint

        waypoints[2*i-2] = prev_waypoint
        connecting_clusters[2*i-2] = (prev_clust, current_clust)

        waypoints[2*i-1] = next_waypoint
        connecting_clusters[2*i-1] = (current_clust, next_clust)
    end

    waypoints[2*(n_cluster_seqs-2)+2] = (sequence.lon[end], sequence.lat[end])
    connecting_clusters[2*(n_cluster_seqs-2)+2] = (sequence.id[end], sequence.id[end])

    return DataFrame(waypoint=waypoints, connecting_clusters=connecting_clusters)
end

"""
    optimize_waypoints!(
        clusters::Vector{Cluster},
        ms_soln::MothershipSolution,
        tender_soln::Vector{TenderSolution},
        problem::Problem,
        waypoint_optim_method,
        time_limit::Float64,
        wpt_optim_plot_flag::Bool=false,
    )::Tuple{Vector{Cluster},MothershipSolution,Vector{TenderSolution}}

Optimize the waypoints of the mothership route using the provided optimization method.

# Arguments
- `clusters`: Current set of clusters
- `ms_soln`: Current mothership solution
- `tender_soln`: Current tender solutions
- `problem`: Problem instance
- `waypoint_optim_method`: Function to use in waypoint optimization.
- `time_limit`: Time limit for waypoint optimization, in seconds
- `wpt_optim_plot_flag`: Flag to enable plotting during waypoint optimization.

# Returns
- Updated clusters, mothership solution, and tender solutions after waypoint optimization
"""
function optimize_waypoints!(
    clusters::Vector{Cluster},
    ms_soln::MothershipSolution,
    tender_soln::Vector{TenderSolution},
    problem::Problem,
    waypoint_optim_method,
    time_limit::Float64,
    wpt_optim_plot_flag::Bool=false,
)::Tuple{Vector{Cluster},MothershipSolution,Vector{TenderSolution}}
    if isnothing(waypoint_optim_method)
        return clusters, ms_soln, tender_soln
    end

    @info "Optimizing full waypoint subset using $(waypoint_optim_method)"

    solution_temp::MSTSolution = optimize_waypoints(
        MSTSolution([clusters], [ms_soln], [tender_soln]),
        problem,
        waypoint_optim_method;
        time_limit=time_limit,
        plot_flag=wpt_optim_plot_flag,
    )

    clusters = solution_temp.cluster_sets[1]
    ms_soln = solution_temp.mothership_routes[1]
    tender_soln = solution_temp.tenders[1]

    return clusters, ms_soln, tender_soln
end

"""
    optimize_waypoints(
        soln::MSTSolution,
        problem::Problem,
        opt_method;
        cost_tol::Float64=0.0,
        gradient_tol::Float64=3e4,
        iterations::Int64=typemax(Int64),
        time_limit::Float64=200.0,
        plot_flag::Bool,
    )::MSTSolution
    optimize_waypoints(
        soln::MSTSolution,
        problem::Problem,
        opt_method,
        candidate_wpt_idxs::Union{AbstractVector{<:Integer},AbstractUnitRange{<:Integer}};
        cost_tol::Float64=0.0,
        gradient_tol::Float64=3e4,
        iterations::Int64=typemax(Int64),
        time_limit::Float64=200.0,
        plot_flag::Bool,
    )::MSTSolution

Optimize waypoint positions in the final mothership route using `opt_method` provided
to minimize the critical path time while avoiding exclusion zones.

This function optimizes the waypoints from the last mothership route (excluding depot start/
end points)by perturbing their positions within a constrained search space. The optimization
minimizes the critical path score while penalizing waypoints that fall within exclusion
zones.

Disturbance scenarios can be handled by optimizing only a subset of waypoints corresponding
to clusters affected by disturbances. `candidate_wpt_idxs` specifies the indices of
waypoints to optimize, allowing for partial optimization of the route.

# Arguments
- `soln::MSTSolution`: Current mothership and tender solution.
- `problem::Problem`: Definition of problem context (vessel specs and exclusion zones)
- `opt_method`: Optimization algorithm used to improve waypoint positions.
- `candidate_wpt_idxs`: Indices of waypoints to optimize. If not provided, all interior
(non-depot) waypoints will be optimized.
- `cost_tol`: Relative cost function tolerance. Stop optimization if function falls below
this value.
- `gradient_tol::Float64`: Gradient norm convergence tolerance. Stop optimization when
gradient norm falls below this value
- `iterations::Int64`: Maximum number of optimization iterations
- `time_limit::Float64`: Soft limit on the time spent optimizing
- `plot_flag::Bool`: Flag to enable plotting during waypoint optimization.

# Returns
Updated mothership and tender solution with optimized waypoint positions and regenerated
tender sorties
"""
function optimize_waypoints(
    soln::MSTSolution,
    problem::Problem,
    opt_method;
    cost_tol::Float64=0.0,
    gradient_tol::Float64=3e4,
    iterations::Int64=typemax(Int64),
    time_limit::Float64=200.0,
    plot_flag::Bool,
)::MSTSolution
    # Optimize all interior waypoints by default
    n_nodes = length(soln.mothership_routes[end].route.nodes)
    n_nodes <= 0 && return soln

    # Candidate waypoint indices to optimize (excluding depot start/end)
    candidate_wpt_idxs = 2:n_nodes-1

    # Pass to partial-optimization overload with all interior indices
    return optimize_waypoints(
        soln, problem, opt_method,
        candidate_wpt_idxs;
        cost_tol=cost_tol,
        gradient_tol=gradient_tol,
        iterations=iterations,
        time_limit=time_limit,
        plot_flag=plot_flag
    )
end
function optimize_waypoints(
    soln::MSTSolution,
    problem::Problem,
    opt_method,
    candidate_wpt_idxs::Union{AbstractVector{<:Integer},AbstractUnitRange{<:Integer}};
    cost_tol::Float64=0.0,
    gradient_tol::Float64=3e4,
    iterations::Int64=typemax(Int64),
    time_limit::Float64=200.0,
    plot_flag::Bool,
)::MSTSolution
    exclusions_mothership::POLY_VEC = problem.mothership.exclusion.geometry
    exclusions_tender::POLY_VEC = problem.tenders.exclusion.geometry
    vessel_weightings::NTuple{2,AbstractFloat} = (
        problem.mothership.weighting, problem.tenders.weighting
    )

    last_ms_route = soln.mothership_routes[end]

    # Flatten initial waypoints
    waypoints_initial::Vector{Point{2,Float64}} = last_ms_route.route.nodes
    n_wpts::Int64 = length(waypoints_initial)
    free_idxs::Vector{Int} = sort!(unique!(collect(candidate_wpt_idxs)))

    isempty(free_idxs) && return soln
    @assert all(1 .<= free_idxs .<= n_wpts) "Free waypoint indices must be within 1:$n_wpts"

    """
    Pack waypoints into optimization vector x, and unpack back into waypoints
    """
    @inline function pack_x(pts::Vector{Point{2,Float64}}, idxs::Vector{Int})::Vector{Float64}
        x = Vector{Float64}(undef, 2 * length(idxs))
        @inbounds for (j, i) in enumerate(idxs)
            x[2j-1] = pts[i][1]
            x[2j] = pts[i][2]
        end
        return x
    end

    best_soln::MSTSolution = soln
    best_score::Float64 = critical_path(soln, vessel_weightings)
    best_count::Int64 = 0
    soln_proposed::MSTSolution = soln

    fig_wpts = plot_flag ? Plot.debug_waypoints(problem, waypoints_initial) : nothing

    # Objective from x -> critical_path
    function obj(
        x::Vector{Float64};
    )::Float64
        # Rebuild waypoints
        @inbounds for (j, i) in enumerate(free_idxs)
            wpts[i] = Point{2,Float64}(x[2j-1], x[2j])
        end

        # Exit early here if any waypoints are inside exclusion zones
        has_bad_waypoint::Vector{Bool} = point_in_exclusion.(wpts, Ref(exclusions_mothership))
        exclusion_count::Int = count(has_bad_waypoint)
        if exclusion_count > 0
            naive_score = maximum(vessel_weightings) *
                          sum(haversine.(wpts[1:end-1], wpts[2:end]))
            return naive_score * (exclusion_count + 1)^2
        end

        soln_proposed = rebuild_solution_with_waypoints(
            soln,
            wpts,
            exclusions_mothership,
            exclusions_tender
        )

        score::Float64 = critical_path(soln_proposed, vessel_weightings)
        isinf(score) && throw(DomainError(score,
            "Critical path cost is infinite, indicating a waypoint in an exclusion zone."
        ))

        if score < best_score
            best_soln = soln_proposed
            best_score = score
            best_count = 0
        else
            # For debugging and tracking
            best_count += 1
        end
        plot_flag && Plot.scatter_by_id!(fig_wpts.current_axis[], wpts[2:end-1])

        return score
    end

    # Run Optim with the provided optimization method.
    opt_options = Optim.Options(
        iterations=iterations,
        f_reltol=cost_tol,
        g_abstol=gradient_tol,
        # show_trace=true, show_every=10,
        store_trace=true,
        allow_f_increases=false,  # allow or disallow objective function value to increase
        time_limit=time_limit
    )

    # Waypoints to optimize, formatted for Optim
    wpts::Vector{Point{2,Float64}} = waypoints_initial
    x0::Vector{Float64} = pack_x(wpts, free_idxs)

    # Ensure bounds are within lat/lon limits
    (lon_min, lon_max, lat_min, lat_max) = get_bbox_bounds([
        exclusions_mothership,
        exclusions_tender,
        problem.targets.points.geometry,
        [problem.depot]
    ])

    # Bound longitude and latitude bounding boxes to entire problem instance
    lb::Vector{Float64} = repeat([lon_min, lat_min], length(x0) ÷ 2)
    ub::Vector{Float64} = repeat([lon_max, lat_max], length(x0) ÷ 2)

    if opt_method isa Optim.ParticleSwarm
        opt_method = Optim.ParticleSwarm(
            lower=lb,
            upper=ub,
            n_particles=opt_method.n_particles,
        )
    end
    result = optimize(
        obj,
        lb,
        ub,
        x0,
        opt_method,
        opt_options
    )

    x_best_flat::Vector{Float64} = Optim.minimizer(result)
    x_best::Vector{Point{2,Float64}} = Point.(
        [
        waypoints_initial[1],  # first waypoint (depot)
        tuple.(x_best_flat[1:2:end], x_best_flat[2:2:end])...,
        waypoints_initial[end]  # last waypoint (depot)
    ]
    )
    best_soln = rebuild_solution_with_waypoints(
        best_soln,
        x_best,
        exclusions_mothership,
        exclusions_tender
    )

    @info "Type:" summary(result)
    @info "Minimum value:" minimum(result)
    if Optim.converged(result)
        @info "Converged at: "
        Optim.x_converged(result) && @info "x"
        Optim.f_converged(result) && @info "f(x)"
        Optim.g_converged(result) && @info "∇f(x)"
        @info "\nConverged after $(Optim.iterations(result)) iterations"
    else
        @info "Did not converge after $(Optim.iterations(result)) iterations"
    end
    @info "Best critical path score found: $best_score"

    if plot_flag
        Plot.route!(
            fig_wpts.current_axis[], best_soln.mothership_routes[end], labels=true, color=:black
        )
        target_points = problem.targets.points.geometry
        Plot.scatter!(
            fig_wpts.current_axis[], target_points, color=:black, markersize=10, marker=:x
        )
        display(fig_wpts)
        Plot.solution!(
            fig_wpts.current_axis[],
            problem,
            best_soln;
            highlight_critical_path=true,
        )
        display(fig_wpts)
        result_trace = Optim.trace(result)
        display(Plot.trace(result_trace, opt_method))
    end

    return best_soln
end

"""
    rebuild_solution_with_waypoints(
        solution_ex::MSTSolution,
        waypoints_proposed::Vector{Point{2,Float64}},
        exclusions_mothership::POLY_VEC,
        exclusions_tender::POLY_VEC
    )::MSTSolution

Rebuild full `MSTSolution` with updated waypoints and other attributes.

# Arguments
- `solution_ex`: The existing MSTSolution to be updated.
- `waypoints_proposed`: Vector of proposed waypoints to update the mothership route.
- `exclusions_mothership`: Exclusion zone polygons for the mothership.
- `exclusions_tender`: Exclusion zone polygons for the tenders.

# Returns
A new MSTSolution with the updated mothership route and tender sorties.
"""
function rebuild_solution_with_waypoints(
    solution_ex::MSTSolution,
    waypoints_proposed::Vector{Point{2,Float64}},
    exclusions_mothership::POLY_VEC,
    exclusions_tender::POLY_VEC
)::MSTSolution
    dist_vector_proposed, line_strings_proposed = get_feasible_vector(
        waypoints_proposed, exclusions_mothership
    )
    ms_route_new::Route = Route(
        copy(waypoints_proposed),
        dist_vector_proposed,
        vcat(line_strings_proposed...)
    )

    # Create a new MothershipSolution with the updated route
    ms_soln_new = MothershipSolution(
        solution_ex.mothership_routes[end].cluster_sequence,
        ms_route_new
    )

    # Update the tender solutions with the new waypoints
    tender_soln_new = generate_tender_sorties(
        solution_ex,
        waypoints_proposed,
        exclusions_tender
    )

    # Return a new MSTSolution with the updated mothership route
    return MSTSolution(
        solution_ex.cluster_sets,
        [solution_ex.mothership_routes[1:end-1]; ms_soln_new],
        [solution_ex.tenders[1:end-1]; [tender_soln_new]]
    )
end
