
struct Route
    nodes::Vector{Point{2,Float64}}
    dist_matrix::Matrix{Float64}
    line_strings::Vector{LineString{2,Float64}}
end

@kwdef struct MothershipSolution
    cluster_sequence::DataFrame
    route::Route
end

struct TenderSolution
    id::Int64
    start::Point{2,Float64}
    finish::Point{2,Float64}
    sorties::Vector{Route}
    dist_matrix::Matrix{Float64}
end

function TenderSolution(t::TenderSolution, sorties::Vector{Route})
    return TenderSolution(
        t.id,
        t.start,
        t.finish,
        sorties,
        t.dist_matrix
    )
end

struct MSTSolution
    cluster_sets::Vector{Vector{Cluster}}
    mothership_routes::Vector{MothershipSolution}
    tenders::Vector{Vector{TenderSolution}}
end

"""
    create_exclusion_zones(env_constraint::Raster, threshold::Float64)::Raster

Create exclusion zones based on environmental raster data and vessel threshold.

# Arguments
- `env_constraint`: Environmental constraint raster.
- `threshold`: Threshold for given vessel's environmental constraint.

# Returns
Exclusion zones for environmental constraint and vessel threshold provided.
"""
function create_exclusion_zones(env_constraint::Raster, threshold::Float64)::Raster
    # TODO: kwargs to specify max/min threshold values
    return (env_constraint .!== env_constraint.missingval) .&& (env_constraint .>= threshold)
end

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
    waypoint_geom = AG.createpoint(waypoint[1], waypoint[2])
    convex_hulls::POLY_VEC = AG.convexhull.(exclusions)
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
    pt = AG.creategeom(AG.wkbPoint)

    # iterate safely within all vertices on the exterior ring to populate pts
    @inbounds for i in 1:n_points
        pt = AG.getpoint(exterior_ring, i - 1) # 0-based indexing
        boundary_points[i] = Point(pt[1], pt[2])
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
    # TODO: Implement convex hull exclusion zones
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
    # TODO: Implement convex hull exclusion zones
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

using Optim

function optimize_waypoints(
    soln::MSTSolution,
    problem::Problem;
    tolerance::Float64=5e5, # stop when gradient norm below this
    iterations::Int=10  # treated as max iterations for Optim
)::MSTSolution
    exclusions_mothership::POLY_VEC = problem.mothership.exclusion.geometry
    exclusions_tender::POLY_VEC = problem.tenders.exclusion.geometry
    n_tenders::Int8 = problem.tenders.number
    t_cap::Int16 = problem.tenders.capacity
    vessel_weightings::NTuple{2,AbstractFloat} = (
        problem.mothership.weighting, problem.tenders.weighting
    )

    last_ms_route = soln.mothership_routes[end]

    # Flatten initial waypoints
    waypoints_initial::Vector{Point{2,Float64}} = last_ms_route.route.nodes[2:end-1]
    x0::Vector{Float64} = reinterpret(Float64, waypoints_initial)

    # Objective from x -> critical_path
    function obj(x::Vector{Float64}; penalty=100.0)::Float64
        # Rebuild waypoints
        wpts::Vector{Point{2,Float64}} = Point.(
            [
            last_ms_route.route.nodes[1],  # first waypoint (depot)
            tuple.(x[1:2:end], x[2:2:end])...,
            last_ms_route.route.nodes[end]  # last waypoint (depot)
        ])

        soln_proposed = rebuild_solution_with_waypoints(soln, wpts, exclusions_mothership)

        score = critical_path(soln_proposed, vessel_weightings)
        if any(point_in_exclusion.(wpts, Ref(exclusions_mothership)))
            # Penalise if any proposed waypoint is in an exclusion zone
            score *= penalty
        end

        return score
    end

    # Run Optim with simple gradient descent
    opt_options = Optim.Options(
        iterations=iterations,
        g_tol=tolerance,
        show_trace=true
    )

    result = optimize(
        obj,
        x0,
        GradientDescent(),
        opt_options
    )

    x_best_flat = Optim.minimizer(result)
    x_best::Vector{Point{2,Float64}} = Point.(
        [
        last_ms_route.route.nodes[1],  # first waypoint (depot)
        tuple.(x_best_flat[1:2:end], x_best_flat[2:2:end])...,
        last_ms_route.route.nodes[end]  # last waypoint (depot)
    ]
    )

    # Rebuild solution with the optimal waypoints
    soln_opt::MSTSolution = rebuild_solution_with_waypoints(soln, x_best, exclusions_mothership)

    # Regenerate tender sorties
    clust_seq::Vector{Int64} = last_ms_route.cluster_sequence.id[2:end-1]
    pairs::Vector{NTuple{2,Point{2,Float64}}} = tuple.(x_best[2:2:end-1], x_best[3:2:end-1])
    soln_opt.tenders[end] = tender_sequential_nearest_neighbour.(
        soln_opt.cluster_sets[end],
        pairs,
        Ref(n_tenders),
        Ref(t_cap),
        Ref(exclusions_tender)
    )

    return soln_opt
end

"""
    rebuild_solution_with_waypoints(
        solution_ex::MSTSolution,
        waypoints_proposed::Vector{Point{2,Float64}},
        exclusions_mothership::POLY_VEC
    )::MSTSolution

Rebuild full `MSTSolution` with updated waypoints and other attributes.

# Arguments
- `solution_ex`: The existing MSTSolution to be updated.
- `waypoints_proposed`: Vector of proposed waypoints to update the mothership route.
- `exclusions_mothership`: Exclusion zone polygons for the mothership.

# Returns
A new MSTSolution with the updated mothership route and tender sorties.
"""
function rebuild_solution_with_waypoints(
    solution_ex::MSTSolution,
    waypoints_proposed::Vector{Point{2,Float64}},
    exclusions_mothership::POLY_VEC
)::MSTSolution
    # Update the waypoints in the mothership route
    dist_vector_proposed, line_strings_proposed = get_feasible_vector(
        waypoints_proposed, exclusions_mothership
    )
    n = length(dist_vector_proposed)
    # Ensure the distance matrix is square and matches the number of waypoints
    dist_matrix_proposed = zeros(Float64, n + 1, n + 1)
    diag_idxs = CartesianIndex.(1:n, 2:n+1)
    dist_matrix_proposed[diag_idxs] .= dist_vector_proposed
    ms_route_new::Route = Route(
        waypoints_proposed,
        dist_matrix_proposed,
        vcat(line_strings_proposed...)
    )

    # Create a new MothershipSolution with the updated route
    ms_soln_new = MothershipSolution(
        solution_ex.mothership_routes[end].cluster_sequence,
        ms_route_new
    )

    # Update the tender solutions with the new waypoints
    tender_soln_new = generate_tender_sorties(solution_ex, waypoints_proposed)

    # Return a new MSTSolution with the updated mothership route
    return MSTSolution(
        solution_ex.cluster_sets,
        [ms_soln_new],
        [tender_soln_new]
    )
end

"""
    generate_tender_sorties(
        soln::MSTSolution,
        tmp_wpts::Vector{Point{2,Float64}}
    )::Vector{TenderSolution}

Generate tender sorties based on the updated waypoints, including all updated attributes.

# Arguments
- `soln`: The existing MSTSolution containing tenders.
- `tmp_wpts`: Vector of temporary waypoints to update the tender sorties.

# Returns
Updated TenderSolution objects with new waypoints.
"""
function generate_tender_sorties(
    soln::MSTSolution,
    tmp_wpts::Vector{Point{2,Float64}}
)::Vector{TenderSolution}
    # Update the tender solutions with the new waypoints
    tender_soln_ex = soln.tenders[end]

    ids = getfield.(tender_soln_ex, :id)
    starts_ex = getfield.(tender_soln_ex, :start)
    finishes_ex = getfield.(tender_soln_ex, :finish)
    sorties = getfield.(tender_soln_ex, :sorties)
    dist_mats = getfield.(tender_soln_ex, :dist_matrix)

    js = 1:length(tender_soln_ex)

    starts_new = tmp_wpts[2 .* js]
    finishes_new = tmp_wpts[2 .* js.+1]

    # mask which tenders need updating
    mask = (starts_ex .!= starts_new) .|| (finishes_ex .!= finishes_new)

    sorties_new::Vector{Vector{Route}} = Vector{Vector{Route}}(undef, length(sorties[mask]))
    for (i, s) in enumerate(sorties[mask])
        new_routes = Vector{Route}(undef, length(s))
        for (j, r) in enumerate(s)
            line_strings_updated = deepcopy(r.line_strings)
            # Recalculate line strings from waypoint to first node and last node to finish
            if r.line_strings[1].points[1] != starts_new[i] # start changes
                continue_from = findfirst(
                    ==(r.nodes[1]),
                    first.(getproperty.(line_strings_updated, :points))
                )

                # re-calc linestrings for the first segment: start to first node
                new_segment_flattened = [LineString{2,Float64}(
                    [starts_new[i], r.nodes[1]]
                )]
                line_strings_updated = [
                    new_segment_flattened;
                    line_strings_updated[continue_from:end]
                ]
            end
            if r.line_strings[end].points[end] != finishes_new[i] # finish changes
                keep_until = findlast(
                    ==(r.nodes[end]),
                    last.(getproperty.(line_strings_updated, :points))
                )
                # re-calc linestrings for the last segment: last node to finish
                new_segment_flattened = [LineString{2,Float64}(
                    [r.nodes[end], finishes_new[i]]
                )]
                line_strings_updated = [
                    line_strings_updated[1:keep_until];
                    new_segment_flattened
                ]
            end
            new_routes[j] = Route(
                r.nodes,
                r.dist_matrix,
                line_strings_updated
            )
        end
        sorties_new[i] = new_routes
    end

    tender_soln_new = copy(tender_soln_ex)
    tender_soln_new[mask] = TenderSolution.(
        ids[mask],
        starts_new[mask],
        finishes_new[mask],
        sorties_new,
        dist_mats[mask]
    )
    return tender_soln_new
end

"""
    nearest_neighbour(
        cluster_centroids::DataFrame,
        exclusions_mothership::POLY_VEC,
        exclusions_tender::POLY_VEC,
    )::MothershipSolution
    nearest_neighbour(
        cluster_centroids::DataFrame,
        exclusions_mothership::POLY_VEC,
        exclusions_tender::POLY_VEC,
        current_point::Point{2,Float64},
        ex_ms_route::MothershipSolution,
        cluster_seq_idx::Int64
    )::MothershipSolution

Apply the nearest neighbor algorithm:
- starting from the depot (1st row/col) and returning to the depot, or
- starting from the current point and returning to the depot.

# Arguments
- `cluster_centroids`: DataFrame containing id, lat, lon. Depot has `id=0` in row 1.
- `exclusions_mothership`: Exclusion zone polygons for mothership.
- `exclusions_tender`: Exclusion zone polygons for tenders.
- `current_point`: Point{2, Float64} representing the current location of the mothership.
- `ex_ms_route`: MothershipSolution object containing the existing route.
- `cluster_seq_idx`: Index denoting the mothership position by cluster sequence index.

# Returns
MothershipSolution object containing:
- `cluster_sequence`: Centroid sequence DataFrame (containing id, lat, lon).
- `route`: Route object containing waypoints, distance matrix, and line strings.
    - `nodes`: Vector of Point{2, Float64} waypoints.
    - `dist_matrix`: Distance matrix between centroids.
    - `line_strings`: Vector of LineString objects for each path.
"""
function nearest_neighbour(
    cluster_centroids::DataFrame,
    exclusions_mothership::POLY_VEC,
    exclusions_tender::POLY_VEC,
)::MothershipSolution
    # TODO: Use vector rather than DataFrame for cluster_centroids
    # adjust_waypoints to ensure not within exclusion zones - allows for feasible path calc
    feasible_centroids::Vector{Point{2,Float64}} = adjust_waypoint.(
        Point{2,Float64}.(cluster_centroids.lon, cluster_centroids.lat),
        Ref(exclusions_mothership)
    )

    cluster_centroids[!, :lon] = [pt[1] for pt in feasible_centroids]
    cluster_centroids[!, :lat] = [pt[2] for pt in feasible_centroids]

    # Create distance matrix between feasible nodes - cluster centroids
    dist_matrix = get_feasible_matrix(
        feasible_centroids,
        exclusions_mothership
    )[1]

    tour_length = size(dist_matrix, 1)
    visited = falses(tour_length)
    tour = Vector{Int64}(undef, tour_length + 1)
    total_distance = 0.0

    idx = 1
    current_location = 1
    tour[1] = current_location
    visited[current_location] = true

    while idx < tour_length
        idx += 1
        distances = dist_matrix[current_location, :]
        distances[visited] .= Inf
        nearest_idx = argmin(distances)

        tour[idx] = nearest_idx
        total_distance += dist_matrix[current_location, nearest_idx]
        visited[nearest_idx] = true
        current_location = nearest_idx
    end

    # Return to the depot and adjust cluster_sequence to zero-based indexing
    tour[end] = 1
    total_distance += dist_matrix[current_location, 1]
    tour .-= 1

    ordered_centroids = cluster_centroids[[findfirst(==(id), cluster_centroids.id) for id in tour], :]

    # combine exclusions for mothership and tenders
    exclusions_all = vcat(exclusions_mothership, exclusions_tender)
    waypoints = get_waypoints(ordered_centroids, exclusions_all)

    # Calc feasible path between waypoints.
    waypoint_dist_vector, waypoint_path_vector = get_feasible_vector(
        waypoints.waypoint, exclusions_mothership
    )
    path = vcat(waypoint_path_vector...)

    return MothershipSolution(
        cluster_sequence=ordered_centroids,
        route=Route(waypoints.waypoint, dist_matrix, path)
    )
end
function nearest_neighbour(
    cluster_centroids::DataFrame,
    exclusions_mothership::POLY_VEC,
    exclusions_tender::POLY_VEC,
    current_point::Point{2,Float64},
    ex_ms_route::MothershipSolution,
    cluster_seq_idx::Int64
)::MothershipSolution
    # TODO: Use vector rather than DataFrame for cluster_centroids
    # adjust_waypoints to ensure not within exclusion zones - allows for feasible path calc
    feasible_centroids::Vector{Point{2,Float64}} = adjust_waypoint.(
        Point{2,Float64}.(cluster_centroids.lon, cluster_centroids.lat),
        Ref(exclusions_mothership)
    )
    cluster_centroids[!, :lon] = [pt[1] for pt in feasible_centroids]
    cluster_centroids[!, :lat] = [pt[2] for pt in feasible_centroids]

    # Create distance matrix between start, end, and feasible cluster centroids
    dist_matrix = get_feasible_matrix(
        vcat([current_point], feasible_centroids),
        exclusions_mothership
    )[1]
    centroid_matrix = dist_matrix[3:end, 3:end]
    current_dist_vector = dist_matrix[1, 3:end]
    end_dist_vector = dist_matrix[2, 3:end]

    tour_length = size(centroid_matrix, 1)
    visited = falses(tour_length)
    tour = Vector{Int64}(undef, tour_length)

    current_location = argmin(current_dist_vector)
    total_distance = current_dist_vector[current_location]
    tour[1] = current_location
    visited[current_location] = true
    idx = 1

    while !all(visited)
        idx += 1
        distances = centroid_matrix[current_location, :]
        distances[visited] .= Inf
        nearest_idx = argmin(distances)

        tour[idx] = nearest_idx
        total_distance += distances[nearest_idx]
        visited[nearest_idx] = true
        current_location = nearest_idx
    end

    # Return to the depot and adjust cluster_sequence to zero-based indexing
    push!(tour, 0)
    total_distance += end_dist_vector[current_location]

    sequenced_remaining_centroids = cluster_centroids[tour.+1, :]

    # combine exclusions for mothership and tenders
    exclusions_all = vcat(exclusions_mothership, exclusions_tender)
    waypoints = get_waypoints(
        current_point,
        vcat(DataFrame(ex_ms_route.cluster_sequence[cluster_seq_idx, :]), sequenced_remaining_centroids),
        exclusions_all
    )

    cluster_sequence = vcat(
        ex_ms_route.cluster_sequence[1:cluster_seq_idx, :],
        sequenced_remaining_centroids
    )
    ordered_clusters = sort(cluster_sequence[1:end-1, :], :id)
    ordered_nodes = Point{2,Float64}.(ordered_clusters.lon, ordered_clusters.lat)
    updated_full_dist_matrix = get_feasible_matrix(
        ordered_nodes,
        exclusions_mothership
    )[1]

    # Calc feasible path between waypoints.
    _, waypoint_path_vector = get_feasible_vector(
        waypoints.waypoint, exclusions_mothership
    )
    ex_path = ex_ms_route.route.line_strings[
        1:findfirst(
        x -> x == ex_ms_route.route.nodes[2*cluster_seq_idx-1],
        getindex.(getproperty.(ex_ms_route.route.line_strings, :points), 2)
    )
    ]
    new_path = vcat(waypoint_path_vector...)
    full_path = vcat(ex_path, new_path...)

    return MothershipSolution(
        cluster_sequence=cluster_sequence,
        route=Route(
            vcat(
                ex_ms_route.route.nodes[1:2*cluster_seq_idx-2],
                waypoints.waypoint
            ),
            updated_full_dist_matrix,
            full_path
        )
    )
end

"""
    two_opt(
        ms_soln_current::MothershipSolution,
        exclusions_mothership::POLY_VEC,
        exclusions_tender::POLY_VEC,
    )::MothershipSolution
    two_opt(
        ms_soln_current::MothershipSolution,
        exclusions_mothership::POLY_VEC,
        exclusions_tender::POLY_VEC,
        cluster_seq_idx::Int64
    )::MothershipSolution
    two_opt(
        tender_soln_current::TenderSolution,
        exclusions_tender::POLY_VEC,
    )::TenderSolution

Apply the 2-opt heuristic to improve current routes by uncrossing crossed links between
waypoints for the whole route or between current location and depot.

# Arguments
- `ms_soln_current`: Current MothershipSolution to improve.
- `tender_soln_current`: Current TenderSolution to improve.
- `exclusions_mothership`: Exclusion zone polygons for mothership.
- `exclusions_tender`: Exclusion zone polygons for tenders.
- `cluster_seq_idx`: Index denoting the mothership position by cluster sequence index.

# Returns
MothershipSolution object containing:
    - `cluster_sequence`: DataFrame of cluster centroids in sequence (containing id, lon, lat).
    - `route`: Route object containing waypoints, distance matrix, and line strings.
        - `nodes`: Vector of Point{2, Float64} waypoints.
        - `dist_matrix`: Distance matrix between centroids.
        - `line_strings`: Vector of LineString objects for each path.
"""
function two_opt(
    ms_soln_current::MothershipSolution,
    exclusions_mothership::POLY_VEC,
    exclusions_tender::POLY_VEC,
)::MothershipSolution

    cluster_centroids = ms_soln_current.cluster_sequence
    dist_matrix = ms_soln_current.route.dist_matrix

    # If depot is last row, remove
    if cluster_centroids.id[1] == cluster_centroids.id[end]
        cluster_centroids = cluster_centroids[1:end-1, :]
    end

    # Initialize route as ordered waypoints
    best_route = cluster_centroids.id .+ 1
    best_distance = return_route_distance(best_route, dist_matrix)
    improved = true

    while improved
        improved = false
        for j in 1:(length(best_route)-1)
            for i in (j+1):length(best_route)
                new_route = two_opt_swap(best_route, j, i)
                new_distance = return_route_distance(new_route, dist_matrix)

                if new_distance < best_distance
                    best_route = new_route
                    best_distance = new_distance
                    improved = true
                end
            end
        end
    end

    # Re-orient route to start from and end at the depot, and adjust to zero-based indexing
    best_route = orient_route(best_route)
    push!(best_route, best_route[1])
    best_route .-= 1

    ordered_nodes = cluster_centroids[[findfirst(==(id), cluster_centroids.id) for id in best_route], :]
    exclusions_all = vcat(exclusions_mothership, exclusions_tender)
    waypoints = get_waypoints(ordered_nodes, exclusions_all)

    waypoint_dist_vector, waypoint_path_vector = get_feasible_vector(
        waypoints.waypoint,
        exclusions_mothership
    )

    path = vcat(waypoint_path_vector...)

    return MothershipSolution(
        cluster_sequence=ordered_nodes,
        route=Route(waypoints.waypoint, dist_matrix, path)
    )
end
function two_opt(
    ms_soln_current::MothershipSolution,
    exclusions_mothership::POLY_VEC,
    exclusions_tender::POLY_VEC,
    cluster_seq_idx::Int64
)::MothershipSolution
    # TODO: only apply to the route between the current cluster and the depot
    # TODO: Then concatenate the new route with the existing route
    cluster_centroids = ms_soln_current.cluster_sequence
    dist_matrix = ms_soln_current.route.dist_matrix

    # If depot is last row, remove
    if cluster_centroids.id[1] == cluster_centroids.id[end]
        cluster_centroids = cluster_centroids[1:end-1, :]
    end

    # Initialize route as ordered waypoints
    best_route = cluster_centroids.id .+ 1
    best_distance = return_route_distance(best_route, dist_matrix)
    improved = true

    while improved
        improved = false
        for j in cluster_seq_idx:(length(best_route)-1)
            for i in (j+1):length(best_route)
                new_route = two_opt_swap(best_route, j, i)
                new_distance = return_route_distance(new_route, dist_matrix)

                if new_distance < best_distance
                    best_route = new_route
                    best_distance = new_distance
                    improved = true
                    @debug "Improved by swapping $(j) and $(i)"
                end
            end
        end
    end

    # Re-orient route to start from and end at the depot, and adjust to zero-based indexing
    best_route = orient_route(best_route)
    push!(best_route, best_route[1])
    best_route .-= 1
    best_route = orient_route(
        best_route,
        ms_soln_current.cluster_sequence.id[1:cluster_seq_idx]
    )
    ordered_nodes = cluster_centroids[
        [findfirst(==(id), cluster_centroids.id) for id in best_route], :
    ]
    exclusions_all = vcat(exclusions_mothership, exclusions_tender)
    waypoints = get_waypoints(ordered_nodes, exclusions_all)
    updated_waypoints = vcat(
        ms_soln_current.route.nodes[1:2*cluster_seq_idx-1],
        waypoints.waypoint[2*cluster_seq_idx:end]
    )

    _, waypoint_path_vector = get_feasible_vector(
        updated_waypoints,
        exclusions_mothership
    )

    path = vcat(waypoint_path_vector...)

    return MothershipSolution(
        cluster_sequence=ordered_nodes,
        route=Route(updated_waypoints, dist_matrix, path)
    )
end
function two_opt(
    tender_soln_current::TenderSolution,
    exclusions_tender::POLY_VEC,
)::TenderSolution
    sorties_current = tender_soln_current.sorties
    sorties_new = Vector{Route}(undef, length(sorties_current))

    for (idx, sortie) in enumerate(sorties_current)
        nodes = vcat([tender_soln_current.start], sortie.nodes, [tender_soln_current.finish])
        n = length(nodes)

        if n ≤ 3
            # If there is only 1 node between start and finish, no change is possible
            sorties_new[idx] = sortie
            continue
        end

        seq_best = collect(1:n)
        dist_best = tender_sortie_dist(seq_best, sortie.dist_matrix)
        improved = true

        while improved
            improved = false
            for j in 2:(n-2)
                for i in (j+1):(n-1)
                    seq_proposed = two_opt_swap(seq_best, j, i)
                    dist_proposed = tender_sortie_dist(seq_proposed, sortie.dist_matrix)
                    if dist_proposed < dist_best
                        seq_best, dist_best = seq_proposed, dist_proposed
                        improved = true
                    end
                end
            end
        end

        # build new sortie ordering, ignoring start and finish
        nodes_new = nodes[seq_best[2:end-1]]
        dist_matrix_new = sortie.dist_matrix[seq_best, seq_best]

        # recompute feasible path through exclusions
        new_path = vcat(get_feasible_vector(nodes[seq_best], exclusions_tender)[2]...)

        sorties_new[idx] = Route(nodes_new, dist_matrix_new, new_path)
    end

    #? recompute or delete distance matrix for all TenderSolution with new routes

    return TenderSolution(
        tender_soln_current.id,
        tender_soln_current.start,
        tender_soln_current.finish,
        sorties_new,
        tender_soln_current.dist_matrix,
    )
end

"""
    two_opt_swap(route::Vector{Int64}, i::Int, j::Int)::Vector{Int64}

Swap two points in a route.

# Arguments
- `route`: Vector of cluster indices.
- `i`: Index of the first point to swap.
- `j`: Index of the second point to swap.

# Returns
Route with swapped points.
"""
function two_opt_swap(route::Vector{Int64}, i::Int, j::Int)::Vector{Int64}
    return vcat(route[1:(i-1)], reverse(route[i:j]), route[(j+1):end])
end

"""
    orient_route(route::Vector{Int64})::Vector{Int64}

Orient the route such that it starts with the depot (1).

# Arguments
- `route`: Vector of cluster indices.

# Returns
Route starting with the depot.
"""
function orient_route(route::Vector{Int64})::Vector{Int64}
    idx = findfirst(==(1), route)
    return vcat(route[idx:end], route[1:idx-1])
end
function orient_route(route::Vector{Int64}, start_segment::Vector{Int64})::Vector{Int64}
    if length(start_segment) < 2
        throw(ArgumentError("start_segment must have at least 2 elements"))
    end
    if start_segment[2] == route[2]
        return route
    elseif start_segment[2] == route[end-1]
        return reverse(route)
    else
        throw(ArgumentError("start_segment does not match route"))
    end
    return
end

"""
    tender_sequential_nearest_neighbour(
        cluster::Cluster,
        waypoints::NTuple{2,Point{2,Float64}},
        n_tenders::Int8,
        t_cap::Int16,
        exclusions::Vector{IGeometry{wkbPolygon}}
    )::TenderSolution

Assign nodes to tenders sequentially (stop-by-stop) based on nearest neighbor.

# Arguments
- `cluster`: Cluster object containing nodes.
- `waypoints`: Tuple of start and end waypoints.
- `n_tenders`: Number of tenders.
- `t_cap`: Tender capacity.
- `exclusions`: Exclusion zone polygons.

# Returns
- `solution`: TenderRoutingSolution object containing:
    - `cluster_id`: Cluster ID.
    - `start`: Start waypoint.
    - `finish`: End waypoint.
    - `sorties`: Vector of Sortie objects containing nodes and sortie distance.
        - `nodes`: Vector of Point{2, Float64} waypoints.
        - `dist_matrix`: Distance matrix between nodes.
        - `line_strings`: Vector of LineString objects for each path.
"""
function tender_sequential_nearest_neighbour(
    cluster::Cluster,
    waypoints::NTuple{2,Point{2,Float64}},
    n_tenders::Int8,
    t_cap::Int16,
    exclusions::POLY_VEC
)::TenderSolution
    nodes = [[waypoints[1]]; cluster.nodes]
    full_nodes = vcat(nodes, [waypoints[2]])

    # Compute the full feasible matrix and associated paths once.
    dist_matrix, path_matrix = get_feasible_matrix(full_nodes, exclusions)

    tender_tours = [Vector{Int}(undef, t_cap) for _ in 1:n_tenders]
    tour_lengths = zeros(Int, n_tenders)
    n_nodes = length(nodes)

    visited = falses(n_nodes)
    visited[1] = true
    visited[2:end] .= isinf.(dist_matrix[1, 2:n_nodes])

    # sequentially assign closest un-visited nodes stop-by-stop, tender-by-tender
    for _ in 1:t_cap
        for t in 1:n_tenders
            if all(visited)
                break
            end

            current_node = tour_lengths[t] == 0 ? 1 : tender_tours[t][tour_lengths[t]]

            #? Will this change the original matrix?
            distances = dist_matrix[current_node, 1:n_nodes]
            distances[visited] .= Inf
            nearest_idx = argmin(distances)

            tour_lengths[t] += 1
            tender_tours[t][tour_lengths[t]] = nearest_idx

            visited[nearest_idx] = true
        end
    end

    # delete excess elements and remove empty tours by slicing tender_tours to tour_lengths
    tender_tours = filter(
        !isempty,
        [tour[1:tour_lengths[t]] for (t, tour) in enumerate(tender_tours)]
    )

    routes = Vector{Route}(undef, length(tender_tours))
    for (t, tour) in enumerate(tender_tours)
        # Route indices incl waypoints: (1), tour nodes, end (last index of full_nodes)
        route_indices = Vector{Int}(undef, length(tour) + 2)
        route_indices[1] = 1
        route_indices[2:end-1] = tour
        route_indices[end] = length(full_nodes)

        # Submatrix for this route
        route_matrix = dist_matrix[route_indices, route_indices]

        route_paths = [
            (i < j ? path_matrix[i, j] : path_matrix[j, i])
            for (i, j) in zip(route_indices[1:end-1], route_indices[2:end])
        ]

        routes[t] = Route(
            cluster.nodes[tour.-1],
            route_matrix,
            vcat(route_paths...)
        )
    end

    initial_sortie = TenderSolution(
        cluster.id,
        waypoints[1],
        waypoints[2],
        routes,
        dist_matrix
    )
    improved_sortie = two_opt(
        initial_sortie,
        exclusions
    )

    return improved_sortie
end
