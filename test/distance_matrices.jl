if !@isdefined(TEST_SOLUTION)
    const (_, TEST_SOLUTION) = initialize_instances()
end

@testset "Ensure consistent distance matrices in sorties and tenders (clusters)" begin
    tenders = TEST_SOLUTION.tenders[end]

    # Get tender distances from sorties
    sortie_dists::Vector{Vector{Float64}} = [
        vcat(getfield.(t.sorties, :dist_matrix)...)
        for t in tenders
    ]

    # Get tender distance matrices for whole clusters
    dist_matrices::Vector{Matrix{Float64}} = getfield.(tenders, :dist_matrix)

    # Check sortie distances in tender distance matrices
    sortie_dist_in_matrices::BitVector = map(
        (u, v) -> all(u .∈ Ref(v)),
        sortie_dists,
        dist_matrices
    )

    @test all(sortie_dist_in_matrices)
end
