
"""
    perturb_swap_solution(
        soln::MSTSolution,
        clust_idx::Int64=-1,
        exclusions::DataFrame = DataFrame();
        enforce_diff_sortie::Bool = false
    )::MSTSolution

Perturb the solution by swapping two nodes in a cluster.

# Arguments
- `soln`: Solution to perturb.
- `clust_idx`: Index of the cluster to perturb. Default = -1 randomly selects a cluster.
- `exclusions`: DataFrame of exclusions. Default = DataFrame().
- `enforce_diff_sortie`: Boolean flag to enforce different sorties for node swaps.

# Returns
Perturbed full solution.
"""
function perturb_swap_solution(
    soln::MSTSolution,
    clust_idx::Int64=-1,
    exclusions::DataFrame = DataFrame();
    enforce_diff_sortie::Bool = false
)::MSTSolution
    # Choose a random cluster (assume fixed clustering) if none provided
    clust_idx = clust_idx == -1 ? rand(1:length(soln.tenders)) : clust_idx

    # Copy the tender solution to perturb
    tender = deepcopy(soln.tenders[clust_idx])
    sorties = tender.sorties

    # If < 2 sorties in cluster, no perturbation possible
    if length(sorties) < 2
        return soln
    end

    # Choose two random sorties to swap nodes between
    # if enforce_diff_sortie is true: ensure different sorties are selected
    sortie_a_idx, sortie_b_idx = enforce_diff_sortie ?
        shuffle(1:length(sorties))[1:2] :
        rand(1:length(sorties), 2)

    sortie_a, sortie_b = sorties[sortie_a_idx], sorties[sortie_b_idx]

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
            node_a_idx, node_b_idx = shuffle(1:sortie_length)[1:2]
        end
    else
        # Chose two random node indices from different sorties
        node_a_idx = rand(1:length(sortie_a.nodes))
        node_b_idx = rand(1:length(sortie_b.nodes))
    end

    # Swap the nodes between the two sorties
    node_a, node_b = sortie_a.nodes[node_a_idx], sortie_b.nodes[node_b_idx]
    sortie_a.nodes[node_a_idx] = node_b
    sortie_b.nodes[node_b_idx] = node_a

    # TODO:Re-run two-opt on the modified sorties
    # Recompute the feasible paths for the modified sorties
    updated_tender_tours::Vector{Vector{Point{2, Float64}}} = [
        [[tender.start]; sortie.nodes; [tender.finish]]
        for sortie in [sortie_a, sortie_b]
    ]

    # Get new linestrings
    linestrings::Vector{Vector{Vector{LineString{2, Float64}}}} = getindex.(
        get_feasible_vector.(updated_tender_tours, Ref(exclusions)),
        2
    )

    # Update the modified sorties
    sorties[sortie_a_idx] = Route(
        sortie_a.nodes,
        sortie_a.dist_matrix,
        vcat(linestrings[1]...)
    )
    sorties[sortie_b_idx] = Route(
        sortie_b.nodes,
        sortie_b.dist_matrix,
        vcat(linestrings[2]...)
    )

    tenders_all::Vector{TenderSolution} = soln.tenders
    tenders_all[clust_idx] = tender

    return MSTSolution(soln.cluster_sets, soln.mothership_routes, tenders_all)
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
- `soln_init`: Initial solution.
- `objective_function`: Function to evaluate the solution.
- `perturb_function`: Function to perturb the solution.
- `exclusions`: DataFrame of exclusions. Default = DataFrame().
- `max_iterations`: Maximum number of iterations. Default = 5_000.
- `temp_init`: Initial temperature. Default = 500.0.
- `cooling_rate`: Rate of cooling to guide acceptance probability for SA algorithm. Default = 0.95 = 95%.
- `static_limit`: Number of iterations to allow stagnation before early exit. Default = 150.

# Returns
- `soln_best`: Best solution::MSTSolution found.
- `obj_best`: Objective value of the best solution.
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
    soln_best = deepcopy(soln_init)
    obj_best = objective_function(soln_init)
    # TODO: Display cost values for each cluster, rather than total solution cost
    for clust_idx in 1:length(soln_init.cluster_sets[end])
        # Initialize current solution as best, reset temp
        temp = temp_init
        soln_current = deepcopy(soln_best)
        obj_current = obj_best
        static_ctr = 0

        @info "\nCluster $(soln_init.cluster_sets[end][clust_idx].id)"
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
                    soln_best = deepcopy(soln_current)
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
