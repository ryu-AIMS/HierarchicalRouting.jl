
@testset "initial_solution()" begin

    problem = HierarchicalRouting.load_problem("output_slopes_3-10m.geojson")

    @testset "Clustering: process_problem()" begin

        clusters, cluster_centroids_df = HierarchicalRouting.process_problem(problem)

        @test "Check ..." begin

        end

        @test "Check ..." begin

        end

    end

    @testset "Feasible path logic" begin

        # Check points not contained within exclusion zone polygons
        # Check points contained within exclusion zone polygons
        # Check points in convex hull

        # Check paths between points are feasible - no intersecting polygons
        # Check paths between points are shortest
        # ? Check paths between points are symmetrical

        # Check closest_crossed_polygon returns correct polygon
        # Check find_widest_points returns correct points

        @test "find_widest_points()" begin
            widest_points, polygon_idxs = HierarchicalRouting.find_widest_points(
                problem.mothership.depot,
                problem.tenders.depot,
                problem.mothership.exclusion
            )
        end

        @test "closest_crossed_polygon()" begin
            closest_polygon, polygon_idx = HierarchicalRouting.closest_crossed_polygon(
                problem.mothership.depot,
                problem.tenders.depot,
                problem.mothership.exclusion
            )
        end

        @test "is_visible()" begin
            # Check points not contained within exclusion zone polygons
            # Check points contained within exclusion zone polygons
            # Check points in convex hull
            # Check points along 'touching' exclusion zone boundaries
            # Check points on vertices
            # Check points on edges

            is_visible = HierarchicalRouting.is_visible(
                problem.mothership.depot,
                problem.tenders.depot,
                problem.mothership.exclusion
            )
        end
    end

    @testset "Graph/Network building" begin
        # Check build_graph() contains start and end points
        # Check for scenario with final point in exclusion zone convex hull
        # Check for scenario with initial point in exclusion zone convex hull
        # Check for scenario with both (initial, final) points in exclusion zone convex hull

        @test "build_graph()" begin
        end

        @test "build_network()" begin
        end

        @test "a_star()" begin
        end

    end

    @testset "Mothership - initial solution: nearest_neighbour()" begin

        clusters, cluster_centroids_df = HierarchicalRouting.process_problem(problem)

        ms_soln_NN = HierarchicalRouting.nearest_neighbour(cluster_centroids_df, problem.mothership.exclusion, problem.tenders.exclusion)

        @test typeof(ms_soln_NN) == MSTSolution

    end

    @testset "two_opt()" begin

        clusters, cluster_centroids_df = HierarchicalRouting.process_problem(problem)

        ms_soln_NN = HierarchicalRouting.nearest_neighbour(cluster_centroids_df, problem.mothership.exclusion, problem.tenders.exclusion)

        ms_soln_2opt = HierarchicalRouting.two_opt(ms_soln_NN, problem.mothership.exclusion, problem.tenders.exclusion)

        @test typeof(ms_soln_2opt) == MSTSolution

    end

    @testset "tender_sequential_nearest_neighbour()" begin

        clusters, cluster_centroids_df = HierarchicalRouting.process_problem(problem)

        ms_soln_NN = HierarchicalRouting.nearest_neighbour(cluster_centroids_df, problem.mothership.exclusion, problem.tenders.exclusion)

        ms_soln_2opt = HierarchicalRouting.two_opt(ms_soln_NN, problem.mothership.exclusion, problem.tenders.exclusion)

        clust_seq = [i for i in ms_soln_2opt.cluster_sequence.id if i!==0 && i <= length(clusters)]
        tender_soln = HierarchicalRouting.TenderSolution[]

        for (i, cluster_id) in enumerate(clust_seq)
            start_waypoint =  ms_soln_2opt.route.nodes[2 * i]
            end_waypoint =  ms_soln_2opt.route.nodes[2 * i + 1]
            @info "$(i): Clust $(cluster_id) from $(start_waypoint) to $(end_waypoint)"

            t_solution = HierarchicalRouting.tender_sequential_nearest_neighbour(
                clusters[cluster_id],
                (start_waypoint, end_waypoint),
                problem.tenders.number, problem.tenders.capacity, problem.tenders.exclusion
            )

            push!(tender_soln, t_solution)
        end

        @test "Check the type of tender_soln" begin
            @test typeof(tender_soln) == Vector{HierarchicalRouting.TenderSolution}
        end

    end

end
