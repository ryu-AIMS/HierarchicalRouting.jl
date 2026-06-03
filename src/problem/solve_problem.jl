
"""
    initial_solution(
        problem::Problem;
        k::Int=1,
        disturbance_clusters::Set{Int64}=Set{Int64}(),
        waypoint_optim_method=nothing,
        seed::Union{Nothing,Int64}=nothing,
        rng::AbstractRNG=Random.GLOBAL_RNG,
    )::MSTSolution

Generate an initial solution to the problem for:
- the mothership by using the nearest neighbour heuristic
and then improving using the 2-opt heuristic, and
-  for tenders using the sequential nearest neighbour heuristic.
Optionally, optimize waypoints using a provided optimization method.

# Arguments
- `problem`: Problem instance to solve
- `k`: Number of clusters to generate
- `disturbance_clusters`: Set of sequenced clusters to simulate disturbances before.
- `waypoint_optim_method`: Function to use in waypoint optimization.
- `seed`: Optional seed for random number generation
- `rng`: Random number generator

# Returns
Best total MSTSolution found
"""
function initial_solution(
    problem::Problem;
    k::Int=1,
    disturbance_clusters::Set{Int64}=Set{Int64}(),
    waypoint_optim_method=nothing,
    seed::Union{Nothing,Int64}=nothing,
    rng::AbstractRNG=Random.GLOBAL_RNG,
)::MSTSolution
    return solve(
        problem;
        k,
        disturbance_clusters,
        seed,
        rng,
        waypoint_optim_method,
        do_improve=false
    )
end

"""
    solve(
        problem::Problem;
        seed::Union{Nothing,Int64}=nothing,
        rng::AbstractRNG=Random.GLOBAL_RNG,
        k::Int=1,
        cluster_iterations::Int=1000,
        cluster_restarts::Int=20,
        disturbance_clusters::Set{Int64}=Set{Int64}(),
        do_improve::Bool=true,
        info_log::Bool=true,
        soln_progress_plot_flag::Bool=false,
        waypoint_optim_method=nothing,
        time_limit::Real=200.0,
        wpt_optim_plot_flag::Bool=false,
        temp_init::Float64=0.25,
        cooling_rate::Float64=0.9,
        min_iters::Int=1000,
        static_limit::Int=1,
        max_iterations::Int=typemax(Int),
        sa_improve_plot_flag::Bool=true,
        output_dir::String="",
    )::MSTSolution

Generate a solution to the problem for:
- the mothership by using the nearest neighbour heuristic to generate an initial solution,
    and then improving using the 2-opt heuristic, and
-  for tenders using the sequential nearest neighbour heuristic.
Optimize solution using simulated annealing as default optimization method.
Optionally, optimize waypoints using a set or provided optimization method.

# Arguments
- `problem`: Problem instance to solve
- `seed`: Optional seed for random number generation
- `rng`: AbstractRNG for random number generation
- `k`: Number of clusters to generate
- `cluster_iterations`: Number of iterations to perform in clustering step
- `cluster_restarts`: Number of restarts to perform in clustering step
- `disturbance_clusters`: Set of sequenced clusters to simulate disturbances before.
- `do_improve`: Whether to improve the initial solution by optimization tender sorties
- `info_log::Bool`: Flag to switch info statement logging
- `soln_progress_plot_flag`: Flag to plot solution progress for debugging/visualization
- `waypoint_optim_method`: Function to use in waypoint optimization.
- `time_limit`: Time limit for waypoint optimization, in seconds
- `wpt_optim_plot_flag`: Flag to plot waypoint optimization for debugging/visualization
- `temp_init`: Initial temperature for simulated annealing
- `cooling_rate`: Cooling rate for simulated annealing
- `min_iters`: Minimum number of iterations to perform before allowing early exit
- `static_limit`: Number of iterations to allow stagnation before early exit
- `max_iterations`: Maximum number of iterations. Default = typemax(Int).
- `sa_improve_plot_flag`: Flag to plot simulated annealing iteratively improved solutions
- `output_dir::String`: Path to output directory. If empty, do not save outputs.

# Returns
Best total MSTSolution found
"""
function solve(
    problem::Problem;
    seed::Union{Nothing,Int64}=nothing,
    rng::Union{Nothing,AbstractRNG}=nothing,
    k::Int=1,
    cluster_iterations::Int=1000,
    cluster_restarts::Int=20,
    disturbance_clusters::Set{Int64}=Set{Int64}(),
    do_improve::Bool=true,
    info_log::Bool=true,
    soln_progress_plot_flag::Bool=false,
    waypoint_optim_method=nothing,
    time_limit::Real=200.0,
    wpt_optim_plot_flag::Bool=false,
    temp_init::Float64=0.25,
    cooling_rate::Float64=0.9,
    min_iters::Int=1000,
    static_limit::Int=1,
    max_iterations::Int=typemax(Int),
    sa_improve_plot_flag::Bool=true,
    output_dir::String="",
)::MSTSolution
    rng = !isnothing(seed) ? Random.MersenneTwister(seed) :
          (!isnothing(rng) ? rng : Random.GLOBAL_RNG)
    output_to_file::Bool = output_dir == "" ? false : true

    if output_to_file
        output_dir = output_dir == "" ? "$seed" : output_dir
        mkpath("$output_dir")
    end

    # Cluster the problem data
    clusters::Vector{Cluster} = cluster_problem(
        problem;
        k,
        max_iter=cluster_iterations,
        n_restarts=cluster_restarts,
    )
    cluster_centroids_df::DataFrame = generate_cluster_df(clusters, problem.depot)

    n_clusters::Int = length(clusters)
    ordered_disturbances::Vector{Int64} = sort(
        [d for d in disturbance_clusters if d <= n_clusters]
    )

    # Route the mothership using nearest neighbour and 2-opt
    ms_route::MothershipSolution = optimize_mothership_route(problem, cluster_centroids_df)
    clust_seq::Vector{Int64} = filter(
        c -> c != 0 && c <= n_clusters,
        ms_route.cluster_sequence.id
    )

    # Generate initial tender solutions using sequential nearest neighbour
    initial_tenders::Vector{TenderSolution} = [
        tender_sequential_nearest_neighbour(
            clusters[clust_seq][j],
            (ms_route.route.nodes[2j], ms_route.route.nodes[2j+1]),
            problem.tenders.number,
            problem.tenders.capacity,
            problem.tenders.exclusion.geometry
        )
        for j in 1:length(clust_seq)
    ]

    solution::MSTSolution = MSTSolution([clusters], [ms_route], [initial_tenders])
    info_log && (@info "$(output_dir) Initial solution generated with objective value: " *
                       "$(critical_path(solution, problem))")
    if soln_progress_plot_flag || output_to_file
        fig_initial = Plot.solution(
            problem,
            solution;
            highlight_critical_path_flag=true,
            title="Initial",
            size=(700, 875)
        )
        soln_progress_plot_flag && display(fig_initial)
        output_to_file && CairoMakie.save("$output_dir/1_initial_solution.png", fig_initial)
    end

    if do_improve
        info_log && @info "$(output_dir) Improving solution by perturbing initial sorties"
        # Optimize the initial tenders solution up to the first disturbance
        next_cluster_idx::Int64 = !isempty(ordered_disturbances) ?
                                  ordered_disturbances[1] :
                                  n_clusters

        solution, _ = improve_solution(
            solution,
            problem,
            1,
            next_cluster_idx;
            temp_init,
            cooling_rate,
            min_iters,
            static_limit,
            max_iterations,
            output_dir,
            info_log,
            sa_improve_plot_flag,
            rng,
        )

        # Update cluster sequence after improvement
        clust_seq = filter(
            c -> c != 0 && c <= length(solution.cluster_sets[end]),
            solution.mothership_routes[end].cluster_sequence.id
        )
        if soln_progress_plot_flag || output_to_file
            fig_sa_opt = Plot.solution(
                problem,
                solution;
                highlight_critical_path_flag=true,
                title="SA Optimized",
                size=(700, 875)
            )
            soln_progress_plot_flag && display(fig_sa_opt)
            output_to_file &&
                CairoMakie.save("$output_dir/2_sa_optimized_solution.png", fig_sa_opt)
        end
    end

    # Apply solution to the first set of clusters pre-disturbance
    info_log && @info "$(output_dir) Optimizing waypoints using PSO"
    solution = optimize_waypoints(
        solution,
        problem,
        waypoint_optim_method;
        time_limit=Float64(time_limit),
        plot_flag=wpt_optim_plot_flag,
        output_dir,
        info_log
    )
    if soln_progress_plot_flag || output_to_file
        fig_pso_opt = Plot.solution(
            problem,
            solution;
            highlight_critical_path_flag=true,
            title="PSO Optimized",
            size=(700, 875)
        )
        soln_progress_plot_flag && display(fig_pso_opt)
        output_to_file &&
            CairoMakie.save("$output_dir/3_pso_optimized_solution.png", fig_pso_opt)
    end

    isempty(disturbance_clusters) && return solution

    # Simulate disturbance events
    total_tender_capacity::Int = Int(problem.tenders.number * problem.tenders.capacity)
    time_limit_pso_disturbed = Float64(time_limit / length(ordered_disturbances))
    solution = _apply_disturbance_events!(
        solution,
        clust_seq,
        ordered_disturbances,
        problem,
        total_tender_capacity,
        time_limit_pso_disturbed;
        rng,
        temp_init,
        cooling_rate,
        min_iters,
        static_limit,
        max_iterations,
        output_dir,
        info_log,
        sa_improve_plot_flag,
        do_improve,
        waypoint_optim_method,
        wpt_optim_plot_flag,
    )

    return solution
end

"""
    optimize_mothership_route(
        problem::Problem,
        cluster_centroids_df::DataFrame
    )::MothershipSolution

Generate an optimized mothership route using the nearest neighbour heuristic and 2-opt for:
- the whole mothership route to/from the depot, or
- the remaining/partial mothership route to clusters based on current position.

# Arguments
- `problem`: Problem instance to solve
- `cluster_centroids_df`: DataFrame containing cluster centroids

# Returns
- The optimized mothership route as a `MothershipSolution` object.
"""
function optimize_mothership_route(
    problem::Problem,
    cluster_centroids_df::DataFrame
)::MothershipSolution
    # Nearest Neighbour to generate initial mothership route & matrix
    ms_soln_NN::MothershipSolution = nearest_neighbour(
        cluster_centroids_df,
        problem.mothership.exclusion.geometry,
        problem.tenders.exclusion.geometry
    )

    # 2-opt to improve the NN soln
    ms_soln_2opt::MothershipSolution = two_opt(
        ms_soln_NN,
        problem.mothership.exclusion.geometry,
        problem.tenders.exclusion.geometry
    )
    return ms_soln_2opt
end

"""
    update_mothership_waypoints(
        problem::Problem,
        cluster_centroids_df::DataFrame,
        disturb_clust_idx::Int64,
        ms_route::MothershipSolution,
    )::MothershipSolution

Update unvisited mothership waypoints after a disturbance event, wrt fixed visited segment.

# Arguments
- `problem`: Problem instance to solve
- `cluster_centroids_df`: DataFrame containing cluster centroids
- `disturb_clust_idx`: Index of the disturbance cluster
- `ms_route`: Current mothership route

# Returns
- The updated mothership route as a `MothershipSolution` object.
"""
function update_mothership_waypoints(
    problem::Problem,
    cluster_centroids_df::DataFrame,
    disturb_clust_idx::Int64,
    ms_route::MothershipSolution,
)::MothershipSolution
    # Rebuild waypoints from existing cluster sequence with updated centroids
    cluster_seq_ids::Vector{Int64} = filter(!=(0), ms_route.cluster_sequence.id)

    # Construct vector of Clusters
    clusters::Vector{Cluster} = [
        Cluster(r.id, Point{2,Float64}(r.lon, r.lat), [Point{2,Float64}(r.lon, r.lat)])
        for r in eachrow(cluster_centroids_df) if r.id != 0
    ]
    cluster_seq_df::DataFrame = get_cluster_sequence_df(
        problem.depot,
        cluster_seq_ids,
        clusters
    )

    exclusions_mothership::POLY_VEC = problem.mothership.exclusion.geometry
    exclusions_tender::POLY_VEC = problem.tenders.exclusion.geometry
    exclusions_all::POLY_VEC = vcat(exclusions_mothership, exclusions_tender)

    waypoints::DataFrame = get_waypoints(cluster_seq_df, exclusions_all)

    # Preserve fixed visited route, recompute rest
    n_fixed_wpts::Int64 = 2 * (disturb_clust_idx - 1)
    final_wpts::Vector{Point{2,Float64}} = vcat(
        ms_route.route.nodes[1:n_fixed_wpts],
        waypoints.waypoint[n_fixed_wpts+1:end]
    )
    wpt_dists, wpt_paths = get_feasible_vector(
        final_wpts, exclusions_mothership
    )

    return MothershipSolution(
        cluster_seq_df,
        Route(final_wpts, wpt_dists, vcat(wpt_paths...))
    )
end

"""
    improve_solution(
        initial_solution::MSTSolution,
        problem::Problem,
        current_cluster_idx::Int,
        next_cluster_idx::Int;
        temp_init::Float64,
        cooling_rate::Float64,
        min_iters::Int,
        static_limit::Int,
        rng::AbstractRNG,
        max_iterations::Int=typemax(Int),
        opt_function::Function=simulated_annealing,
        objective_function::Function=critical_path,
        output_dir::String="",
        info_log::Bool=true,
        sa_improve_plot_flag::Bool=true,
    )::Tuple{MSTSolution,Float64}
    improve_solution(
        init_solution::MSTSolution,
        problem::Problem;
        temp_init::Float64,
        cooling_rate::Float64,
        min_iters::Int,
        static_limit::Int,
        rng::AbstractRNG,
        max_iterations::Int=typemax(Int),
        opt_function::Function=simulated_annealing,
        objective_function::Function=critical_path,
        output_dir::String="",
        info_log::Bool=true,
        sa_improve_plot_flag::Bool=true,
    )::Tuple{MSTSolution,Float64}

Improve the solution using the optimization function `opt_function` with the objective
function `objective_function` to improve full and partial solutions.

# Arguments
- `initial_solution`: Initial solution to improve
- `problem`: Problem instance to solve
- `current_cluster_idx`: Index of the current cluster in the sequence
- `next_cluster_idx`: Index of the next cluster in the sequence
- `temp_init`: Initial temperature for simulated annealing
- `cooling_rate`: Cooling rate for simulated annealing
- `min_iters`: Minimum number of iterations to perform before allowing early exit
- `static_limit`: Number of iterations to allow stagnation before early exit
- `max_iterations`: Maximum number of iterations. Default = typemax(Int).
- `opt_function`: Optimization function to improve the solution.
    Default = `simulated_annealing`
- `objective_function`: Objective function to quantify and evaluate the solution.
    Default = `critical_path`
- `output_dir::String`: Path to output directory. If empty, do not save outputs.
- `info_log::Bool`: Flag to switch info statement logging
- `sa_improve_plot_flag`: Flag to plot simulated annealing solution perturbations

# Returns
- `soln_best`: Solution with lowest objective value found
- `z_best`: Objective value associated with the best solution
"""
function improve_solution(
    initial_solution::MSTSolution,
    problem::Problem,
    current_cluster_idx::Int,
    next_cluster_idx::Int;
    temp_init::Float64,
    cooling_rate::Float64,
    min_iters::Int,
    static_limit::Int,
    rng::AbstractRNG,
    max_iterations::Int=typemax(Int),
    opt_function::Function=simulated_annealing,
    objective_function::Function=critical_path,
    output_dir::String="",
    info_log::Bool=true,
    sa_improve_plot_flag::Bool=true,
)::Tuple{MSTSolution,Float64}
    current_mothership_route::MothershipSolution = initial_solution.mothership_routes[end]

    cluster_seq_ids::Vector{Int64} = current_mothership_route.cluster_sequence.id

    exclusions_all = vcat(
        problem.mothership.exclusion.geometry,
        problem.tenders.exclusion.geometry,
    )

    soln_best_partial, z_best = opt_function(
        problem,
        initial_solution,
        objective_function,
        max_iterations,
        temp_init,
        cooling_rate,
        min_iters,
        static_limit,
        rng;
        perturb_idxs=current_cluster_idx:(next_cluster_idx-1),
        output_dir,
        info_log,
        plot_flag=sa_improve_plot_flag,
    )

    merged_clusters = soln_best_partial.cluster_sets[end]
    sort!(merged_clusters, by=c -> c.id)

    # Update full mothership route with the optimized partial route
    depot = current_mothership_route.route.nodes[1]
    interior_ids = @view cluster_seq_ids[2:end-1]
    cluster_seq_df = get_cluster_sequence_df(
        depot, collect(interior_ids), merged_clusters
    )
    all_wpts_df = get_waypoints(cluster_seq_df, exclusions_all)

    # Preserve pre-visited waypoints exactly; only recompute from current_cluster_idx onwards
    n_fixed = 2 * (current_cluster_idx - 1) + 1
    fixed_waypoints = current_mothership_route.route.nodes[1:n_fixed]
    new_waypoints = all_wpts_df.waypoint[2*current_cluster_idx:end]

    final_waypoints = current_cluster_idx > 1 ?
                      vcat(fixed_waypoints, new_waypoints) :
                      all_wpts_df.waypoint

    wpt_dists, wpt_paths = get_feasible_vector(
        final_waypoints, problem.mothership.exclusion.geometry
    )
    updated_ms_soln = MothershipSolution(
        cluster_seq_df,
        Route(final_waypoints, wpt_dists, vcat(wpt_paths...))
    )

    starts = updated_ms_soln.route.nodes[2 .* eachindex(interior_ids)]
    finishes = updated_ms_soln.route.nodes[2 .* eachindex(interior_ids).+1]

    ordered_tenders = [
        _reconcile_tender(
            soln_best_partial.tenders[end][k], starts[k], finishes[k],
            problem.tenders.exclusion.geometry
        )
        for k in eachindex(interior_ids)
    ]

    soln_best = MSTSolution(
        [merged_clusters],
        [updated_ms_soln],
        [ordered_tenders]
    )

    return soln_best, z_best
end
function improve_solution(
    init_solution::MSTSolution,
    problem::Problem;
    temp_init::Float64,
    cooling_rate::Float64,
    min_iters::Int,
    static_limit::Int,
    rng::AbstractRNG,
    max_iterations::Int=typemax(Int),
    opt_function::Function=simulated_annealing,
    objective_function::Function=critical_path,
    output_dir::String="",
    info_log::Bool=true,
    sa_improve_plot_flag::Bool=true,
)::Tuple{MSTSolution,Float64}
    current_cluster_idx::Int64 = 1
    next_cluster_idx::Int64 = length(init_solution.cluster_sets[end]) + 1

    return improve_solution(
        init_solution,
        problem,
        current_cluster_idx,
        next_cluster_idx;
        opt_function,
        objective_function,
        max_iterations,
        temp_init,
        cooling_rate,
        min_iters,
        static_limit,
        output_dir,
        info_log,
        sa_improve_plot_flag,
        rng,
    )
end

function _reconcile_tender(
    t::TenderSolution,
    start::Point{2,Float64},
    finish::Point{2,Float64},
    exclusions::POLY_VEC
)::TenderSolution
    t.start == start && t.finish == finish && return t
    return TenderSolution(t.id, start, finish,
        [_build_sortie_route(s.nodes, start, finish, exclusions) for s in t.sorties]
    )
end
