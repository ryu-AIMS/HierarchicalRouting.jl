
function show(io::IO, ::MIME"text/plain", p::Problem)
    ntargets = nrow(p.targets.points)
    ntendr = p.tenders.number
    tcap = p.tenders.capacity
    vessel_weighting = (p.mothership.weighting, p.tenders.weighting)
    println(io, "Problem:")
    println(io, "\tDepot:\t\t$(round.(p.depot; digits=4))")
    println(io, "\tTargets:\t$ntargets locations")
    println(io, "\tExclusion zones:\n\t\t" *
                "\t$(nrow(p.mothership.exclusion)) mothership\n\t\t" *
                "\t$(nrow(p.tenders.exclusion)) tender"
    )
    println(io, "\tTenders:\t$ntendr boats × capacity $tcap")
    println(io, "\tVessel weightings:\n\t\t\t" *
                "Mothership:\t$(vessel_weighting[1])\n\t\t\t" *
                "Tenders:\t$(vessel_weighting[2])"
    )
end

function show(io::IO, ::MIME"text/plain", c::Cluster)
    pts = length(c.nodes)
    print(io, "Cluster(id=$(c.id), centroid=$(round.(c.centroid; digits=4)), points=$pts)")
end
function show(io::IO, ::MIME"text/plain", clusters::Vector{Cluster})
    n = length(clusters)
    println(io, "$n clusters:")
    for c in clusters
        show(io, "text/plain", c)
        println(io)
    end
end

function show(io::IO, ::MIME"text/plain", mothership::MothershipSolution)
    # nclust = nrows(mothership.cluster_sequence)
    sequence::DataFrame = mothership.cluster_sequence

    println(io, "MothershipSolution:")
    println(io, "\tWaypoints: $(length(mothership.route.nodes)) nodes")
    println(io, "\tMothership sequence: $(sequence.id)")
end

function show(io::IO, ::MIME"text/plain", soln::MSTSolution)
    nclust = length(soln.cluster_sets[end])
    ntend = length(soln.tenders[end])

    println(io, "MSTSolution:" *
                "\n\tClusters: $nclust" *
                "\n\tTenders: $ntend sorties"
    )

    show(stdout, "text/plain", soln.mothership_routes[end])
end
