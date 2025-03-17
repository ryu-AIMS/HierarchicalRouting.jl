
"""
    perturb_swap_solution(
        soln::MSTSolution;
        clust_idx::Int=-1,
        exclusions::DataFrame = DataFrame()
    )

Perturb the solution by swapping two nodes in a cluster.

# Arguments
- `soln` : Solution to perturb.
- `clust_idx` : Index of the cluster to perturb. Default = -1.
- `exclusions` : DataFrame of exclusions. Default = DataFrame().

# Returns
Perturbed full solution.
"""
function perturb_swap_solution(
    soln::MSTSolution;
    clust_idx::Int=-1,
    exclusions::DataFrame = DataFrame()
)
    # Choose a random cluster (assume fixed clustering) if none provided
    clust_idx = clust_idx == -1 ? rand(1:length(soln.tenders)) : clust_idx

    # Copy the tender solution to perturb
    tender = deepcopy(soln.tenders[clust_idx])

    # If < 2 sorties in cluster, no perturbation possible
    if length(tender.sorties) < 2
        return soln
    end

    # Choose two random sorties from the cluster
    sortie_a_idx, sortie_b_idx = rand(1:length(tender.sorties), 2)
    sortie_a, sortie_b = tender.sorties[sortie_a_idx], tender.sorties[sortie_b_idx]

    # No perturbation possible if a sortie has no nodes
    if isempty(sortie_a.nodes) || isempty(sortie_b.nodes)
        return soln
    end

    # Determine node indices to swap
    if sortie_a_idx == sortie_b_idx
        # Ensure two distinct node indices are selected from same sortie
        sortie_length = length(sortie_a.nodes)
        if sortie_length < 2
            return soln
        elseif sortie_length == 2
            node_a_idx, node_b_idx = 1, 2
        else
            node_a_idx, node_b_idx = rand(1:sortie_length, 2)
            while node_a_idx == node_b_idx
                node_b_idx = rand(1:sortie_length)
            end
        end
    else
        # Chose two random node indices from different sorties
        node_a_idx = rand(1:length(sortie_a.nodes))
        node_b_idx = rand(1:length(sortie_b.nodes))
    end

    # Swap the nodes between the two sorties
    node_a, node_b = sortie_a.nodes[node_a_idx], sortie_b.nodes[node_b_idx]
    tender.sorties[sortie_a_idx].nodes[node_a_idx] = node_b
    tender.sorties[sortie_b_idx].nodes[node_b_idx] = node_a

    # TODO:Re-run two-opt on the modified sorties
    # Recompute the feasible paths for the modified sorties
    updated_tender_tours = [
        [[tender.start]; sortie.nodes; [tender.finish]]
        for sortie in [tender.sorties[sortie_a_idx], tender.sorties[sortie_b_idx]]
    ]

    # Update linestrings for the modified sorties
    tender.sorties[sortie_a_idx] = Route(
        tender.sorties[sortie_a_idx].nodes,
        tender.sorties[sortie_a_idx].dist_matrix,
        vcat(get_feasible_vector(updated_tender_tours[1], exclusions)[2]...)
    )
    tender.sorties[sortie_b_idx] = Route(
        tenders[clust_idx].sorties[sortie_b_idx].nodes,
        tenders[clust_idx].sorties[sortie_b_idx].dist_matrix,
        vcat(get_feasible_vector(updated_tender_tours[2], exclusions)[2]...)
    )

    tenders_all = [i == clust_idx ? tender : soln.tenders[i] for i in 1:length(soln.tenders)]
    new_soln = MSTSolution(soln.clusters, soln.mothership, tenders_all)

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
        cooling_rate::Float64 = 0.95,
        static_limit::Int = 150
    )

Simulated Annealing optimization algorithm to optimize the solution.

# Arguments
- `soln_init` : Initial solution.
- `objective_function` : Function to evaluate the solution.
- `perturb_function` : Function to perturb the solution.
- `exclusions` : DataFrame of exclusions. Default = DataFrame().
- `max_iterations` : Maximum number of iterations. Default = 5_000.
- `temp_init` : Initial temperature. Default = 500.0.
- `cooling_rate` : Rate of cooling to guide acceptance probability for SA algorithm. Default = 0.95 = 95%.
- `static_limit` : Number of iterations to allow stagnation before early exit. Default = 150.

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
    cooling_rate::Float64 = 0.95,
    static_limit::Int = 150
)

    # Initialize best solution as initial
    soln_best = soln_init
    obj_best = objective_function(soln_init)

    # TODO: Display cost values for each cluster, rather than total solution cost
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
                @info "$iteration\t\t$obj_best\t$temp"
                @info "Early exit at iteration $iteration due to stagnation."
                break
            end
        end
    end

    @info "\nFinal Value: $obj_best"
    @info "Δ: $(objective_function(soln_init) - obj_best)"
    return soln_best, obj_best
end
