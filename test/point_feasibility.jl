if !@isdefined(TEST_PROBLEM) || !@isdefined(TEST_SOLUTION)
    const (TEST_PROBLEM, TEST_SOLUTION) = initialize_instances()
end

@testset "Ensure point feasibility" begin
    exclusions_mothership = TEST_PROBLEM.mothership.exclusion.geometry
    exclusions_tenders = TEST_PROBLEM.tenders.exclusion.geometry

    waypoints = TEST_SOLUTION.mothership_routes[end].route.nodes

    # No waypoints in mothership exclusion zones
    @test all(.!HierarchicalRouting.point_in_exclusion.(
        waypoints,
        Ref(exclusions_mothership)
    ))

    # No waypoints in tender exclusion zones
    @test all(.!HierarchicalRouting.point_in_exclusion.(
        waypoints,
        Ref(exclusions_tenders)
    ))
end

@testset "Ensure mothership route nodes match linestrings" begin
    ms_route = TEST_SOLUTION.mothership_routes[end].route
    @test all(
        ms_route.nodes[1:end-1] .∈ Ref(
            getindex.(getfield.(ms_route.line_strings, :points), 1)
        )
    )
end
