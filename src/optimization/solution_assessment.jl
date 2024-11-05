
function perturb_solution(soln::MSTSolution, clust_idx::Int=-1)
    new_soln = deepcopy(soln)

    if clust_idx == -1
        # Choose ONE random cluster - assume fixed clustering
        clust_idx = rand(1:length(new_soln.tenders))
    end
    cluster = new_soln.tenders[clust_idx]

    if length(cluster.sorties) < 2
        return new_soln # No perturbation possible if a cluster has less than 2 sorties
    end

    # Choose TWO random sorties from the cluster
    sortie_a_idx, sortie_b_idx = rand(1:length(cluster.sorties), 2)
    # TODO: Allow same sortie if not applying two-opt
    while sortie_a_idx == sortie_b_idx
        sortie_b_idx = rand(1:length(cluster.sorties))
    end
    sortie_a, sortie_b = cluster.sorties[sortie_a_idx], cluster.sorties[sortie_b_idx]

    if isempty(sortie_a.nodes) || isempty(sortie_b.nodes)
        return new_soln # No perturbation possible if a sortie has no nodes
    end

    node_a_idx, node_b_idx = rand(1:length(sortie_a.nodes)), rand(1:length(sortie_b.nodes))
    node_a, node_b = sortie_a.nodes[node_a_idx], sortie_b.nodes[node_b_idx]

    # Swap the nodes between the two sorties
    new_soln.tenders[clust_idx].sorties[sortie_a_idx].nodes[node_a_idx] = node_b
    new_soln.tenders[clust_idx].sorties[sortie_b_idx].nodes[node_b_idx] = node_a

    # TODO:Re-run two-opt on the modified sorties

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
        max_iterations::Int = 100_000,
        temp_init::Float64 = 500.0,
        cooling_rate::Float64 = 0.99_99
    )

Simulated Annealing optimization algorithm to optimize the solution.

# Arguments
- `soln_init` : Initial solution.
- `objective_function` : Function to evaluate the solution.
- `perturb_function` : Function to perturb the solution.
- `max_iterations` : Maximum number of iterations. Default = 10_000.
- `temp_init` : Initial temperature. Default = 500.0.
- `cooling_rate` : Rate of cooling to guide acceptance probability for SA algorithm. Default = 0.995 = 99.95%.
- `static_limit` : Number of iterations to allow stagnation before early exit. Default = 5000.

# Returns
- `soln_best` : Best solution::MSTSolution found.
- `obj_best` : Objective value of the best solution.
"""
function simulated_annealing(
    soln_init::MSTSolution,
    objective_function::Function,
    perturb_function::Function,
    max_iterations::Int = 100_000,
    temp_init::Float64 = 500.0,
    cooling_rate::Float64 = 0.99_99,
    static_limit::Int = 5000
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
            soln_proposed = perturb_function(soln_current, clust_idx)
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
