
"""
    perturb_swap_solution(
        soln::MSTSolution,
        clust_seq_idx_target::Int64=-1,
        exclusions::DataFrame = DataFrame();
        enforce_diff_sortie::Bool = false
    )::MSTSolution
    perturb_swap_solution(
        soln::MSTSolution,
        cluster_pair::Tuple{Int, Int},
        exclusions::DataFrame = DataFrame()
    )::MSTSolution

Perturb the solution by swapping two nodes:
- within a cluster, or
- between two clusters.

# Arguments
- `soln`: Solution to perturb.
- `clust_seq_idx_target`: Sequence index of cluster to perturb. Default = -1: randomly
    selects cluster.
- `exclusions`: DataFrame of exclusions. Default = DataFrame().
- `enforce_diff_sortie`: Boolean flag to enforce different sorties for node swaps.

# Returns
Perturbed full solution.
"""
function perturb_swap_solution(
    soln::MSTSolution,
    clust_seq_idx_target::Int64=-1,
    exclusions::DataFrame = DataFrame();
    enforce_diff_sortie::Bool = false
)::MSTSolution
    clust_seq_idx = clust_seq_idx_target == -1 ?
        rand(1:length(soln.tenders)) :
        clust_seq_idx_target

    tender = deepcopy(soln.tenders[clust_seq_idx])
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
        end
        node_a_idx, node_b_idx = sortie_length == 2 ? (1, 2) : shuffle(1:sortie_length)[1:2]
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
        [[tender.start]; sortie_a.nodes; [tender.finish]],
        [[tender.start]; sortie_b.nodes; [tender.finish]]
    ]

    # Get new linestrings
    linestrings::Vector{Vector{Vector{LineString{2, Float64}}}} = getindex.(
        get_feasible_vector.(updated_tender_tours, Ref(exclusions)),
        2
    )

    # Get the new feasible distance matrices for the modified sorties
    updated_tender_matrices::Vector{Matrix{Float64}} = getindex.(get_feasible_matrix.(
        updated_tender_tours,
        Ref(exclusions)
    ), 1)

    # Update the modified sorties
    sorties[sortie_a_idx] = Route(
        sortie_a.nodes,
        updated_tender_matrices[1],
        vcat(linestrings[1]...)
    )
    sorties[sortie_b_idx] = Route(
        sortie_b.nodes,
        updated_tender_matrices[2],
        vcat(linestrings[2]...)
    )

    tenders_all::Vector{TenderSolution} = copy(soln.tenders)
    tenders_all[clust_seq_idx] = TenderSolution(
        tender.id,
        tender.start,
        tender.finish,
        sorties,
        tender.dist_matrix  #? recompute
    )

    return MSTSolution(soln.cluster_sets, soln.mothership_routes, tenders_all)
end
function perturb_swap_solution(
    soln::MSTSolution,
    cluster_pair::Tuple{Int, Int},
    exclusions::DataFrame = DataFrame()
)::MSTSolution
    clust_a_seq_idx, clust_b_seq_idx = cluster_pair

    tender_a = soln.tenders[clust_a_seq_idx]
    tender_b = soln.tenders[clust_b_seq_idx]

    # Pick random sorties and ensure both have nodes
    sortie_a_idx = rand(1:length(tender_a.sorties))
    sortie_b_idx = rand(1:length(tender_b.sorties))
    sortie_a = deepcopy(tender_a.sorties[sortie_a_idx])
    sortie_b = deepcopy(tender_b.sorties[sortie_b_idx])

    if isempty(sortie_a.nodes) || isempty(sortie_b.nodes)
        return soln
    end

    # Swap two random nodes across the two sorties
    node_a_idx, node_b_idx = rand(1:length(sortie_a.nodes)), rand(1:length(sortie_b.nodes))
    node_a, node_b = sortie_a.nodes[node_a_idx], sortie_b.nodes[node_b_idx]
    sortie_a.nodes[node_a_idx] = node_b
    sortie_b.nodes[node_b_idx] = node_a

    # Update routes for modified sorties
    tours = [
        [[tender_a.start]; sortie_a.nodes; [tender_a.finish]],
        [[tender_b.start]; sortie_b.nodes; [tender_b.finish]]
    ]
    updated_linestrings = getindex.(get_feasible_vector.(tours, Ref(exclusions)), 2)
    updated_sortie_matrices = getindex.(get_feasible_matrix.(tours, Ref(exclusions)), 1)

    sortie_a = Route(
        sortie_a.nodes,
        updated_sortie_matrices[1],
        vcat(updated_linestrings[1]...)
    )
    sortie_b = Route(
        sortie_b.nodes,
        updated_sortie_matrices[2],
        vcat(updated_linestrings[2]...)
    )

    tenders_all::Vector{TenderSolution} = copy(soln.tenders)
    tenders_all[clust_a_seq_idx] = TenderSolution(
        tender_a.id,
        tender_a.start,
        tender_a.finish,
        sorties_a,
        tender_a.dist_matrix
    )
    tenders_all[clust_b_seq_idx] = TenderSolution(
        tender_b.id,
        tender_b.start,
        tender_b.finish,
        sorties_b,
        tender_b.dist_matrix
    )

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
)::Tuple{MSTSolution, Float64}
    # Initialize best solution as initial
    soln_best = deepcopy(soln_init)
    obj_best = objective_function(soln_init)
    cluster_set::Vector{Cluster} = soln_init.cluster_sets[end]
    # TODO: Display cost values for each cluster, rather than total solution cost
    for clust_idx in 1:length(cluster_set)
        # Initialize current solution as best, reset temp
        temp = temp_init
        soln_current = deepcopy(soln_best)
        obj_current = obj_best
        static_ctr = 0

        @info "\n\tCluster\t\t$(cluster_set[clust_idx].id)\n\t" *
            "Iteration \tBest Value \t\tTemp\n\t0\t\t$obj_best\t$temp"

        for iteration in 1:max_iterations
            if rand() < 0.5
                # swap two nodes within the same cluster
                soln_proposed = perturb_function(soln_current, clust_idx, exclusions)
            else
                # swap two nodes between two different random clusters
                clust_swap_idx = shuffle(setdiff(1:length(cluster_set), clust_idx))[1]
                soln_proposed = perturb_function(
                    soln_current,
                    (clust_idx, clust_swap_idx),
                    exclusions
                )
            end
            obj_proposed = objective_function(soln_proposed)
            improvement = obj_current - obj_proposed
            static_ctr += 1

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
                @info "$iteration\t\t$obj_best\t$temp\n\t" *
                    "Early exit at iteration $iteration due to stagnation."
                break
            end
        end
    end

    @info "\nFinal Value: $obj_best\nΔ: $(objective_function(soln_init) - obj_best)"
    return soln_best, obj_best
end
