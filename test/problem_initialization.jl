@testset "Initialize problem" begin
    @testset "Initialize with geojson file" begin
        _target_scenario_folder = "../data/targets/scenarios"
        _target_scenario_name = "output_slopes_3-10m"
        _target_scenario_path=joinpath(
            _target_scenario_folder,
            _target_scenario_name*".geojson"
        )

        problem = HierarchicalRouting.load_problem(_target_scenario_path)
        @test problem isa HierarchicalRouting.Problem
        @test problem.scenario_name == _target_scenario_name
    end
end
