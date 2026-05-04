
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
        k::Int=1,
        disturbance_clusters::Set{Int64}=Set{Int64}(),
        seed::Union{Nothing,Int64}=nothing,
        rng::AbstractRNG=Random.GLOBAL_RNG,
        waypoint_optim_method=nothing,
        do_improve::Bool=true,
        time_limit::Real=200.0,
        wpt_optim_plot_flag::Bool=false,
        cross_cluster_flag::Bool=false,
        soln_progress_plot_flag::Bool=false,
    )::MSTSolution

Generate a solution to the problem for:
- the mothership by using the nearest neighbour heuristic to generate an initial solution,
    and then improving using the 2-opt heuristic, and
-  for tenders using the sequential nearest neighbour heuristic.
Optimize solution using simulated annealing as default optimization method.
Optionally, optimize waypoints using a set or provided optimization method.

# Arguments
- `problem`: Problem instance to solve
- `k`: Number of clusters to generate
- `disturbance_clusters`: Set of sequenced clusters to simulate disturbances before.
- `seed`: Optional seed for random number generation
- `rng`: AbstractRNG for random number generation
- `waypoint_optim_method`: Function to use in waypoint optimization.
- `do_improve`: Whether to improve the initial solution by optimization tender sorties
- `time_limit`: Time limit for waypoint optimization, in seconds
- `wpt_optim_plot_flag`: Flag to plot waypoint optimization for debugging/visualization
- `cross_cluster_flag`: Flag to allow perturbations across clusters in solution improvement
- `soln_progress_plot_flag`: Flag to plot solution progress for debugging/visualization

# Returns
Best total MSTSolution found
"""
function solve(
    problem::Problem;
    k::Int=1,
    disturbance_clusters::Set{Int64}=Set{Int64}(),
    seed::Union{Nothing,Int64}=nothing,
    rng::AbstractRNG=Random.GLOBAL_RNG,
    waypoint_optim_method=nothing,
    do_improve::Bool=true,
    time_limit::Real=200.0,
    wpt_optim_plot_flag::Bool=false,
    cross_cluster_flag::Bool=false,
    soln_progress_plot_flag::Bool=false,
)::MSTSolution
    if !isnothing(seed)
        Random.seed!(rng, seed)
    end

    # Cluster the problem data
    clusters::Vector{Cluster} = cluster_problem(problem; k)
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
    @info "Initial solution generated with objective value: $(critical_path(solution, problem))"
    soln_progress_plot_flag && display(Plot.solution(
        problem,
        solution;
        highlight_critical_path_flag=true,
        title="Initial",
        size=(700, 875)
    ))

    if do_improve
        @info "Improving initial solution using simulated annealing"
        # Optimize the initial tenders solution up to the first disturbance
        next_cluster_idx::Int64 = !isempty(ordered_disturbances) ?
                                  ordered_disturbances[1] :
                                  n_clusters

        solution, _ = improve_solution(
            solution,
            problem,
            1,
            next_cluster_idx;
            cross_cluster_flag
        )

        # Update cluster sequence after improvement
        clust_seq = filter(
            c -> c != 0 && c <= length(solution.cluster_sets[end]),
            solution.mothership_routes[end].cluster_sequence.id
        )
        soln_progress_plot_flag && display(Plot.solution(
            problem,
            solution;
            highlight_critical_path_flag=true,
            title="SA Optimized",
            size=(700, 875)
        ))
    end

    # Apply solution to the first set of clusters pre-disturbance
    @info "Optimizing waypoints using PSO"
    solution = optimize_waypoints(
        solution,
        problem,
        waypoint_optim_method;
        time_limit=Float64(time_limit),
        plot_flag=wpt_optim_plot_flag
    )
    soln_progress_plot_flag && display(Plot.solution(
        problem,
        solution;
        highlight_critical_path_flag=true,
        title="PSO Optimized",
        size=(700, 875)
    ))

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
    optimize_mothership_route(
        problem::Problem,
        cluster_centroids_df::DataFrame,
        cluster_seq_idx::Int64,
        ms_route::MothershipSolution,
        cluster_ids_visited::Vector{Int64}
    )::MothershipSolution

Generate an optimized mothership route using the nearest neighbour heuristic and 2-opt for:
- the whole mothership route to/from the depot, or
- the remaining/partial mothership route to clusters based on current position.

# Arguments
- `problem`: Problem instance to solve
- `cluster_centroids_df`: DataFrame containing cluster centroids
- `cluster_seq_idx`: Index of the cluster sequence
- `ms_route`: Current mothership route
- `cluster_ids_visited`: Vector of cluster IDs that have been visited

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
function optimize_mothership_route(
    problem::Problem,
    cluster_centroids_df::DataFrame,
    cluster_seq_idx::Int64,
    ms_route::MothershipSolution,
    cluster_ids_visited::Vector{Int64}
)::MothershipSolution
    start_point::Point{2,Float64} = ms_route.route.nodes[2*cluster_seq_idx-1]

    remaining_clusters_df::DataFrame = filter(
        row -> row.id ∉ cluster_ids_visited,
        cluster_centroids_df
    )

    # Nearest Neighbour to generate initial mothership route & matrix
    ms_soln_NN::MothershipSolution = nearest_neighbour(
        remaining_clusters_df,
        problem.mothership.exclusion.geometry,
        problem.tenders.exclusion.geometry,
        start_point,
        ms_route,
        cluster_seq_idx
    )

    # 2-opt to improve the NN soln
    ms_soln_2opt::MothershipSolution = two_opt(
        ms_soln_NN,
        problem.mothership.exclusion.geometry,
        problem.tenders.exclusion.geometry,
        cluster_seq_idx
    )

    return ms_soln_2opt
end

"""
    improve_solution(
        initial_solution::MSTSolution,
        problem::Problem,
        current_cluster_idx::Int,
        next_cluster_idx::Int;
        opt_function::Function=simulated_annealing,
        objective_function::Function=critical_path,
        cross_cluster_flag::Bool=true,
        max_iterations::Int=1_000,
        temp_init::Float64=2.0,
        cooling_rate::Float64=0.9,
        min_iters::Int=50,
        static_limit::Int=20,
    )::Tuple{MSTSolution,Float64}
    improve_solution(
        init_solution::MSTSolution,
        problem::Problem;
        opt_function::Function=simulated_annealing,
        objective_function::Function=critical_path,
        cross_cluster_flag::Bool=true,
        max_iterations::Int=1_000,
        temp_init::Float64=2.0,
        cooling_rate::Float64=0.9,
        min_iters::Int=50,
        static_limit::Int=20,
    )

Improve the solution using the optimization function `opt_function` with the objective
function `objective_function` to improve full and partial solutions.

# Arguments
- `initial_solution`: Initial solution to improve
- `problem`: Problem instance to solve
- `current_cluster_idx`: Index of the current cluster in the sequence
- `next_cluster_idx`: Index of the next cluster in the sequence
- `opt_function`: Optimization function to improve the solution
- `objective_function`: Objective function to quantify and evaluate the solution
- `cross_cluster_flag`: Boolean flag to indicate if perturbation across clusters should be
    considered. Default = true.
- `max_iterations`: Maximum number of iterations
- `temp_init`: Initial temperature for simulated annealing
- `cooling_rate`: Cooling rate for simulated annealing
- `min_iters`: Minimum number of iterations to perform before allowing early exit
- `static_limit`: Number of iterations to allow stagnation before early exit

# Returns
- `soln_best`: Solution with lowest objective value found
- `z_best`: Objective value associated with the best solution
"""
function improve_solution(
    initial_solution::MSTSolution,
    problem::Problem,
    current_cluster_idx::Int,
    next_cluster_idx::Int;
    opt_function::Function=simulated_annealing,
    objective_function::Function=critical_path,
    cross_cluster_flag::Bool=true,
    max_iterations::Int=1_000,
    temp_init::Float64=2.0,
    cooling_rate::Float64=0.9,
    min_iters::Int=50,
    static_limit::Int=20,
)::Tuple{MSTSolution,Float64}
    current_mothership_route::MothershipSolution = initial_solution.mothership_routes[end]

    cluster_seq_ids::Vector{Int64} = current_mothership_route.cluster_sequence.id
    final_cluster_idx = next_cluster_idx == length(cluster_seq_ids) - 1 ?
                        next_cluster_idx - 1 :
                        next_cluster_idx

    clust_seq_current::Vector{Int64} = cluster_seq_ids[
        current_cluster_idx+1:final_cluster_idx+1
    ]
    clust_seq_noncurrent::Vector{Int64} = setdiff(cluster_seq_ids, clust_seq_current)
    filter!(!=(0), clust_seq_noncurrent)

    cluster_set::Vector{Cluster} = initial_solution.cluster_sets[end]
    current_clusters::Vector{Cluster} = cluster_set[clust_seq_current]
    sort!(current_clusters, by=c -> c.id)

    tender_set::Vector{TenderSolution} = initial_solution.tenders[end]
    current_tender_routes::Vector{TenderSolution} = tender_set[current_cluster_idx:final_cluster_idx]
    noncurrent_tender_routes::Vector{TenderSolution} = setdiff(
        tender_set,
        current_tender_routes
    )

    current_solution = MSTSolution(
        [current_clusters],
        [current_mothership_route],
        [current_tender_routes]
    )
    soln_best_partial, z_best = opt_function(
        problem,
        current_solution,
        objective_function,
        max_iterations,
        temp_init,
        cooling_rate,
        min_iters,
        static_limit;
        cross_cluster_flag,
    )

    merged_clusters = vcat(
        cluster_set[clust_seq_noncurrent],
        soln_best_partial.cluster_sets[end]
    )
    sort!(merged_clusters, by=c -> c.id)

    merged_tenders = vcat(soln_best_partial.tenders[end], noncurrent_tender_routes)
    sort!(merged_tenders, by=t -> t.id)
    interior_ids = @view cluster_seq_ids[2:end-1]
    ordered_tenders = merged_tenders[interior_ids]

    soln_best = MSTSolution(
        [merged_clusters],
        [soln_best_partial.mothership_routes[end]],
        [ordered_tenders]
    )

    return soln_best, z_best
end
function improve_solution(
    init_solution::MSTSolution,
    problem::Problem;
    opt_function::Function=simulated_annealing,
    objective_function::Function=critical_path,
    cross_cluster_flag::Bool=true,
    max_iterations::Int=1_000,
    temp_init::Float64=2.0,
    cooling_rate::Float64=0.9,
    min_iters::Int=50,
    static_limit::Int=20,
)
    current_cluster_idx::Int64 = 1
    next_cluster_idx::Int64 = length(init_solution.cluster_sets[end])

    return improve_solution(
        init_solution,
        problem,
        current_cluster_idx,
        next_cluster_idx;
        cross_cluster_flag=cross_cluster_flag,
        opt_function,
        objective_function,
        max_iterations,
        temp_init,
        cooling_rate,
        min_iters,
        static_limit
    )
end
