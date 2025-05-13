@testset "Initialize problem" begin
    @testset "Initialize with geojson file" begin
        _target_scenario_path="data/test_target_locations.geojson"
        _subset_path="data/test_target_area.gpkg"
        _bathy_path="data/test_bathy.tif"
        _wave_disturbance_path="data/test_env_disturbances.geojson"

        using GeometryBasics

        depot = GeometryBasics.Point{2, Float64}(146.175, -16.84)
        draft_ms = Float64(-10.0)
        draft_t = Float64(-5.0)
        weight_ms = Float16(5.0)
        weight_t = Float16(2.0)
        n_tenders = Int8(3)
        t_cap = Int16(2)

        problem = HierarchicalRouting.load_problem(
        _target_scenario_path,
        _subset_path,
        _bathy_path,
        _wave_disturbance_path,
        depot,
        draft_ms,
        draft_t,
        weight_ms,
        weight_t,
        n_tenders,
        t_cap;
        debug_mode=true
    );

        @test problem isa HierarchicalRouting.Problem
        @test problem.targets.path == _target_scenario_path
    end
end
