
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
    print(io, "Cluster(id=$(c.id), centroid=$(round.(c.centroid; digits=4)), points=$pts)")
end
function Base.show(io::IO, mime::MIME"text/plain", clusters::Vector{Cluster})
    n = length(clusters)
    println(io, "$n clusters:")
    for c in clusters
        show(io, mime, c)
        print(io, "\n")
    end
end

function Base.show(io::IO, mime::MIME"text/plain", mothership::MothershipSolution)
    # nclust = nrows(mothership.cluster_sequence)
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
