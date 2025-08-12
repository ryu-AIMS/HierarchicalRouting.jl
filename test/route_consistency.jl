@testset "Ensure route consistency" begin
    _target_scenario_path = "data/test_target_locations.geojson"
    _subset_path = "data/test_target_area.gpkg"
    _bathy_path = "data/test_bathy.tif"
    _wave_disturbance_path = "data/test_env_disturbances.geojson"

    depot = (146.175, -16.84)
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
    )

    test_solution = initial_solution(problem; k=8)

    # Helpers
    ## Grab the start/end point of the first/last LineString in a Route
    first_pt(r::HierarchicalRouting.Route) = r.line_strings[1].points[1]
    last_pt(r::HierarchicalRouting.Route) = r.line_strings[end].points[end]

    ## Check if all tender sorties start/end at the same point
    uniform_start(ts::HierarchicalRouting.TenderSolution, ms_nodes) = all(
        first_pt.(ts.sorties) .∈ Ref(ms_nodes)
    )
    uniform_end(ts::HierarchicalRouting.TenderSolution, ms_nodes) = all(
        last_pt.(ts.sorties) .∈ Ref(ms_nodes)
    )

    # Variables
    tenders = test_solution.tenders[end]
    ms_route = test_solution.mothership_routes[end].route
    ms_route_nodes = ms_route.nodes
    ms_line_string_pts = unique(vcat(getfield.(ms_route.line_strings, :points)...))

    # Tests
    ## 1 All mothership route nodes are covered in linestrings
    @test all(ms_route_nodes .∈ Ref(ms_line_string_pts))

    ## 2 All tender sortie linestrings start in the same place for each cluster
    @test all(length(unique(first_pt.(ts.sorties))) == 1 for ts in tenders)

    ## 3 All tender sortie linestrings end in the same place for each cluster
    @test all(length(unique(last_pt.(ts.sorties))) == 1 for ts in tenders)

    ## 4 All tender start points coincide with mothership nodes
    @test all(uniform_start.(tenders, Ref(ms_route_nodes)))

    ## 5 All tender end points coincide with mothership nodes
    @test all(uniform_end.(tenders, Ref(ms_route_nodes)))
end
