
@testset "initial_solution()" begin
    problem = HierarchicalRouting.load_problem("output_slopes_3-10m.geojson")

    @testset "Mothership - initial solution: nearest_neighbour()" begin
        clusters, cluster_centroids_df = HierarchicalRouting.process_problem(problem)
        ms_soln_NN = HierarchicalRouting.nearest_neighbour(
            cluster_centroids_df,
            problem.mothership.exclusion,
            problem.tenders.exclusion
        )
        @test typeof(ms_soln_NN) == MSTSolution
    end

    @testset "two_opt()" begin
        clusters, cluster_centroids_df = HierarchicalRouting.process_problem(problem)
        ms_soln_NN = HierarchicalRouting.nearest_neighbour(
            cluster_centroids_df,
            problem.mothership.exclusion,
            problem.tenders.exclusion
        )
        ms_soln_2opt = HierarchicalRouting.two_opt(
            ms_soln_NN,
            problem.mothership.exclusion,
            problem.tenders.exclusion
        )
        @test typeof(ms_soln_2opt) == MSTSolution
    end

    @testset "tender_sequential_nearest_neighbour()" begin
        solution_init = initial_solution(problem);
        typeof(solution_init) == MSTSolution

        @test all(
            [length(s.sorties) for s in solution_init.tenders[end]] .<= problem.tenders.number
        )
        @test typeof(tender_soln) == Vector{HierarchicalRouting.TenderSolution}
    end
end

@testset "Clustering: process_problem()" begin
    clusters, cluster_centroids_df = HierarchicalRouting.process_problem(problem)
end

# @testset "Feasible path logic" begin

#     # TODO Check points not contained within exclusion zone polygons
#     # TODO Check points contained within exclusion zone polygons
#     # TODO Check points in convex hull
#     # TODO Check paths between points are feasible - no intersecting polygons
#     # TODO Check paths between points are shortest
#     # ? Check paths between points are symmetrical

#     # TODO Check closest_crossed_polygon returns correct polygon
#     # TODO Check find_widest_points returns correct points

#     @testset "find_widest_points()" begin
#         widest_points, polygon_idxs = HierarchicalRouting.find_widest_points(
#             problem.mothership.depot,
#             problem.tenders.depot,
#             problem.mothership.exclusion
#         )
#     end

#     @testset "closest_crossed_polygon()" begin
#         closest_polygon, polygon_idx = HierarchicalRouting.closest_crossed_polygon(
#             problem.mothership.depot,
#             problem.tenders.depot,
#             problem.mothership.exclusion
#         )
#     end

#     @testset "is_visible()" begin
#         # TODO: Check points not contained within exclusion zone polygons
#         # TODO: Check points contained within exclusion zone polygons
#         # TODO: Check points in convex hull
#         # TODO: Check points along 'touching' exclusion zone boundaries
#         # TODO: Check points on vertices
#         # TODO: Check points on edges

#         is_visible = HierarchicalRouting.is_visible(
#             problem.mothership.depot,
#             problem.tenders.depot,
#             problem.mothership.exclusion
#         )
#     end
# end

# @testset "Graph/Network building" begin
#     # TODO Check build_graph() contains start and end points
#     # TODO Check for scenario with final point in exclusion zone convex hull
#     # TODO Check for scenario with initial point in exclusion zone convex hull
#     # TODO Check for scenario with both (initial, final) points in exclusion zone convex hull

#     @test "build_graph()" begin
#     end

#     @test "build_network()" begin
#     end

#     @test "a_star()" begin
#     end

# end
