@testset "Ensure route consistency" begin
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
    tenders = TEST_SOLUTION.tenders[end]
    ms_route = TEST_SOLUTION.mothership_routes[end].route
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
