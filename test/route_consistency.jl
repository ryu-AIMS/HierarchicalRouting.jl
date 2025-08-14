if !@isdefined(TEST_PROBLEM) || !@isdefined(TEST_SOLUTION)
    const (TEST_PROBLEM, TEST_SOLUTION) = initialize_instances()
end

@testset "Ensure route consistency" begin
    # Grab the start/end point of the first/last LineString in a Route
    first_pt(r::HierarchicalRouting.Route) = r.line_strings[1].points[1]
    last_pt(r::HierarchicalRouting.Route) = r.line_strings[end].points[end]

    # Check that the first/last point of each sortie is a point on the mothership route
    uniform_start(ts::HierarchicalRouting.TenderSolution, ms_nodes) = all(
        first_pt.(ts.sorties) .∈ Ref(ms_nodes)
    )
    uniform_end(ts::HierarchicalRouting.TenderSolution, ms_nodes) = all(
        last_pt.(ts.sorties) .∈ Ref(ms_nodes)
    )

    tenders = TEST_SOLUTION.tenders[end]
    ms_route = TEST_SOLUTION.mothership_routes[end].route
    ms_route_nodes = ms_route.nodes
    ms_line_string_pts = unique(vcat(getfield.(ms_route.line_strings, :points)...))

    # All mothership route nodes are visited by linestring points
    @test all(ms_route_nodes .∈ Ref(ms_line_string_pts))

    # Tender sorties start at the same point for each sortie within each cluster
    sortie_start_points = [first_pt.(ts.sorties) for ts in tenders]
    @test all(length.(unique.(sortie_start_points)) .== 1)

    # Tender sorties end at the same point for each sortie within each cluster
    sortie_end_points = [last_pt.(ts.sorties) for ts in tenders]
    @test all(length.(unique.(sortie_end_points)) .== 1)

    # All tender sorties start at a point on the mothership route
    @test all(uniform_start.(tenders, Ref(ms_route_nodes)))

    # All tender sorties end at a point on the mothership route
    @test all(uniform_end.(tenders, Ref(ms_route_nodes)))
end
