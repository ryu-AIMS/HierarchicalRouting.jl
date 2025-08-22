if !@isdefined(TEST_PROBLEM) || !@isdefined(TEST_SOLUTION)
    const (TEST_PROBLEM, TEST_SOLUTION) = initialize_instances()
end

@testset "Ensure point feasibility" begin
    exclusions_mothership = TEST_PROBLEM.mothership.exclusions.geometry
    exclusions_tenders = TEST_PROBLEM.tenders.exclusions.geometry

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
