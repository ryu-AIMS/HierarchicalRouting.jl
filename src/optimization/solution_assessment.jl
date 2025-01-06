
"""
    perturb_swap_solution(soln::MSTSolution, clust_idx::Int=-1, exclusions::DataFrame = DataFrame())

Perturb the solution by swapping two nodes between two sorties in a cluster.

# Arguments
- `soln` : Solution to perturb.
- `clust_idx` : Index of the cluster to perturb. Default = -1.
- `exclusions` : DataFrame of exclusions. Default = DataFrame().

# Returns
- `new_soln` : Perturbed solution.
"""
function perturb_swap_solution(soln::MSTSolution, clust_idx::Int=-1, exclusions::DataFrame = DataFrame())
    new_soln = deepcopy(soln)

    if clust_idx == -1
        # Choose ONE random cluster - assume fixed clustering
        clust_idx = rand(1:length(new_soln.tenders))
    end
    cluster = new_soln.tenders[clust_idx]

    if length(cluster.sorties) < 2
        # No perturbation possible if a cluster has less than 2 sorties
        return new_soln
    end

    # Choose TWO random sorties from the cluster
    sortie_a_idx, sortie_b_idx = rand(1:length(cluster.sorties), 2)
    # TODO: Allow same sortie if not applying two-opt
    while sortie_a_idx == sortie_b_idx
        sortie_b_idx = rand(1:length(cluster.sorties))
    end
    sortie_a, sortie_b = cluster.sorties[sortie_a_idx], cluster.sorties[sortie_b_idx]

    if isempty(sortie_a.nodes) || isempty(sortie_b.nodes)
        # No perturbation possible if a sortie has no nodes
        return new_soln
    end

    node_a_idx, node_b_idx = rand(1:length(sortie_a.nodes)), rand(1:length(sortie_b.nodes))
    node_a, node_b = sortie_a.nodes[node_a_idx], sortie_b.nodes[node_b_idx]

    # Swap the nodes between the two sorties
    new_soln.tenders[clust_idx].sorties[sortie_a_idx].nodes[node_a_idx] = node_b
    new_soln.tenders[clust_idx].sorties[sortie_b_idx].nodes[node_b_idx] = node_a

    # TODO:Re-run two-opt on the modified sorties
    # Recompute the feasible paths for the modified sorties
    updated_tender_tours = [
        [
            [soln.tenders[clust_idx].start];
            sortie.nodes;
            [soln.tenders[clust_idx].finish]
        ]
        for sortie in [
            new_soln.tenders[clust_idx].sorties[sortie_a_idx]; new_soln.tenders[clust_idx].sorties[sortie_b_idx]
        ]
    ]

    feasible_paths = Matrix{
        Tuple{Dict{Int64, Point{2, Float64}}, Vector{SimpleWeightedEdge{Int64, Float64}}}
    }[
        get_feasible_matrix(s, exclusions)[2] for s in updated_tender_tours
    ]

    paths = [get_linestrings(feasible_paths[s], updated_tender_tours[s]) for s in 1:length(feasible_paths)]

    # Update linestrings for the modified sorties
    new_soln.tenders[clust_idx].line_strings[sortie_a_idx] = get_linestrings(feasible_paths[1], updated_tender_tours[1]) # paths[1]
    new_soln.tenders[clust_idx].line_strings[sortie_b_idx] = get_linestrings(feasible_paths[2], updated_tender_tours[2])
      # paths[2]

    # TODO: Recompute the cost for each modified sortie
    # This assumes a function `compute_sortie_cost` exists
    # sortie1.cost = compute_sortie_cost(sortie1)
    # sortie2.cost = compute_sortie_cost(sortie2)

    # Update the tender's total cost if necessary
    # tender.cost = sum(sortie.cost for sortie in tender.sorties)

    return new_soln
end

"""
    simulated_annealing(
        soln_init::MSTSolution,
        objective_function::Function,
        perturb_function::Function,
        exclusions::DataFrame = DataFrame(),
        max_iterations::Int = 5_000,
        temp_init::Float64 = 500.0,
        cooling_rate::Float64 = 0.99,
        static_limit::Int = 250
    )

Simulated Annealing optimization algorithm to optimize the solution.

# Arguments
- `soln_init` : Initial solution.
- `objective_function` : Function to evaluate the solution.
- `perturb_function` : Function to perturb the solution.
- `exclusions` : DataFrame of exclusions. Default = DataFrame().
- `max_iterations` : Maximum number of iterations. Default = 5_000.
- `temp_init` : Initial temperature. Default = 500.0.
- `cooling_rate` : Rate of cooling to guide acceptance probability for SA algorithm. Default = 0.99 = 99%.
- `static_limit` : Number of iterations to allow stagnation before early exit. Default = 250.

# Returns
- `soln_best` : Best solution::MSTSolution found.
- `obj_best` : Objective value of the best solution.
"""
function simulated_annealing(
    soln_init::MSTSolution,
    objective_function::Function,
    perturb_function::Function,
    exclusions::DataFrame = DataFrame(),
    max_iterations::Int = 5_000,
    temp_init::Float64 = 500.0,
    cooling_rate::Float64 = 0.99,
    static_limit::Int = 250
)

    # Initialize best solution as initial
    soln_best = soln_init
    obj_best = objective_function(soln_init)

    for clust_idx in 1:length(soln_init.clusters)
        # Initialize current solution as best, reset temp
        temp = temp_init
        soln_current = soln_best
        obj_current = obj_best
        static_ctr = 0

        @info "\nCluster $(soln_init.clusters[clust_idx].id)"
        @info "Iteration \tBest Value \t\tTemp"
        @info "0\t\t$obj_best\t$temp"

        for iteration in 1:max_iterations
            soln_proposed = perturb_function(soln_current, clust_idx, exclusions)
            obj_proposed = objective_function(soln_proposed)

            static_ctr += 1
            improvement = obj_current - obj_proposed

            # If the new solution is improved OR meets acceptance prob criteria
            if improvement > 0 || exp(improvement / temp) > rand()
                soln_current = soln_proposed
                obj_current = obj_proposed

                if obj_current < obj_best
                    static_ctr = 0
                    soln_best = soln_current
                    obj_best = obj_current
                    @info "$iteration\t\t$obj_best\t$temp"
                end
            end

            temp *= cooling_rate

            if iteration % 100 == 0
                @info "$iteration\t\t$obj_best\t$temp"
            end

            if static_ctr >= static_limit
                @info "Early exit at iteration $iteration due to stagnation. temp=$temp"
                break
            end
        end
    end

    @info "Final Value: $obj_best"
    @info "Δ: $(objective_function(soln_init) - obj_best)"
    return soln_best, obj_best
end
