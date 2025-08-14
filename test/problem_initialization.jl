if !@isdefined(TEST_PROBLEM)
    const TEST_PROBLEM, _ = initialize_instances()
end

@testset "Initialize problem" begin
    @testset "Initialize with geojson file" begin
        @test TEST_PROBLEM isa HierarchicalRouting.Problem
        @test TEST_PROBLEM.targets.path == TARGET_SCENARIO_PATH
    end
end
