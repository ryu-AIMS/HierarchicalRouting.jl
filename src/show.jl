
function Base.show(io::IO, mime::MIME"text/plain", p::Problem)
    ntargets = nrow(p.targets.points)
    ntendr = p.tenders.number
    tcap = p.tenders.capacity
    vessel_weighting = (p.mothership.weighting, p.tenders.weighting)
    println(
        io,
        """
        Problem
        -------

        Depot: $(round.(p.depot; digits=4))
        Deployment Locations: $ntargets
        Exclusion zones:
            $(nrow(p.mothership.exclusion)) mothership
            $(nrow(p.tenders.exclusion)) tender

        Tenders: $ntendr boats
        Capacity: $tcap
        Vessel weightings:
            Mothership: $(vessel_weighting[1])
            Tenders:\t$(vessel_weighting[2])
        """)
end

function Base.show(io::IO, mime::MIME"text/plain", c::Cluster)
    pts = length(c.nodes)
    println(io, "Cluster(id=$(c.id), centroid=$(round.(c.centroid; digits=4)), points=$pts)")
end
function Base.show(io::IO, mime::MIME"text/plain", clusters::Vector{Cluster})
    n = length(clusters)
    println(io, "$n clusters:")
    for c in clusters
        show(io, mime, c)
    end
end

function Base.show(io::IO, mime::MIME"text/plain", mothership::MothershipSolution)
    sequence::DataFrame = mothership.cluster_sequence

    println(
        io,
        """
        MothershipSolution:
            $(length(mothership.route.nodes)) waypoints
            (one entry and exit per cluster, including initial/end)

            Cluster visitation sequence: $(sequence.id)
        """
    )
end

function Base.show(io::IO, mime::MIME"text/plain", tender::TenderSolution)
    println(
        io,
        """
        TenderSolution:
            id: $(tender.id)
            Start location: $(tender.start)
            Finish location: $(tender.finish)
            Number of sorties: $(length(tender.sorties))
        """
    )
end
function Base.show(io::IO, mime::MIME"text/plain", tenders::Vector{TenderSolution})
    n = length(tenders)
    println(io, "$n tenders:")

    for tender in tenders
        show(io, mime, tender)
    end
end

function Base.show(io::IO, mime::MIME"text/plain", soln::MSTSolution)
    nclust = length(soln.cluster_sets[end])
    ntend = length(soln.tenders[end])

    println(
        io,
        """
        MSTSolution:
            Clusters: $nclust
            Tenders: $ntend sorties
        """
    )

    Base.show(io, mime, soln.mothership_routes[end])
end

function Base.show(io::IO, mime::MIME"text/plain", routes::Vector{Route})
    header = "$(length(routes)) Routes"
    println(
        io,
        """
        $header
        $("="^length(header))
        """
    )

    for (i, r) in enumerate(routes)
        println(io, "Route $i")
        show(io, mime, r)

        println("------")
    end
end

function Base.show(io::IO, mime::MIME"text/plain", route::Route)

    waypoints = route.nodes

    for (i, w) in enumerate(waypoints)
        println("Deployment Location $i: $w")
    end

    println(io, "Route to/from deployment locations:")

    # Get output of linestrings exactly as how it would normally be shown in the REPL
    io_capture = IOBuffer()
    show(
        IOContext(io_capture, :color => true, :limit => true),
        MIME("text/plain"),
        GI.coordinates.(route.line_strings)
    )
    println(String(take!(io_capture)))
end
