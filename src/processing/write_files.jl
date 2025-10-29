
"""
    export_points(
        clusters::Vector{Cluster},
        output_path::String="outputs"
    )::DataFrame

Export points to a GeoPackage file, saved in the output directory.

# Arguments
- `clusters`: Vector of clusters.
- `output_path`: Path to the output directory.

# Returns
- `df`: DataFrame with id, geometry, and cluster_id columns.
"""
function export_points(
    clusters::Vector{Cluster},
    output_path::String="outputs"
)::DataFrame
    total_points = sum(length.(getfield.(clusters, :nodes)))
    ids = collect(1:total_points)
    geoms = vcat(getfield.(clusters, :nodes)...)
    cluster_ids = vcat(
        [fill.(getfield(cluster, :id), length(cluster.nodes)) for cluster in clusters]...
    )

    df = DataFrame(id=ids, geometry=geoms, cluster_id=cluster_ids)
    GDF.write(joinpath(output_path, "points.gpkg"), df)
    return df
end

"""
    export_clusters(
        cluster_sequence::DataFrame,
        include_depot::Bool=false,
        output_path::String="outputs",
    )::DataFrame

Export clusters to a GeoPackage file, saved in the output directory.

# Arguments
- `cluster_sequence`: DataFrame with clusters in sequence visited: columns = id, lon, lat.
- `output_path`: Path to the output directory.
- `include_depot`: Boolean indicating whether to include depot as a cluster in the output.

# Returns
- `df`: Cluster sequence with id, geometry, and order_id columns.
"""
function export_clusters(
    cluster_sequence::DataFrame,
    include_depot::Bool=false,
    output_path::String="outputs",
)::DataFrame
    if !include_depot
        cluster_sequence = cluster_sequence[cluster_sequence.id.!=0, :]
    end
    df = DataFrame(
        order_id=include_depot ?
                 (0:size(cluster_sequence, 1)-1) :
                 (1:size(cluster_sequence, 1)),
        cluster_id=cluster_sequence.id,
        geometry=AG.createpoint.(cluster_sequence.lon, cluster_sequence.lat)
    )

    GDF.write(joinpath(output_path, "clusters.gpkg"), df)
    return df
end

"""
    export_exclusions(
        problem::Problem,
        output_path::String="outputs"
    )::Tuple{DataFrame,DataFrame}

Export mothership and tender exclusions (nested in Problem struct) to GeoPackage files,
    saved in the output directory.

# Arguments
- `problem`: Problem instance.
- `output_path`: Path to the output directory.

# Returns
- `exclusions_ms`: Mothership exclusions.
- `exclusions_tenders`: Tenders exclusions.
"""
function export_exclusions(
    problem::Problem,
    output_path::String="outputs"
)::Tuple{DataFrame,DataFrame}
    exclusions_ms = problem.mothership.exclusion
    exclusions_tenders = problem.tenders.exclusion

    # Ensure DataFrames have primary key column (id)
    exclusions_ms.id = collect(1:nrow(exclusions_ms))
    exclusions_tenders.id = collect(1:nrow(exclusions_tenders))

    GDF.write(joinpath(output_path, "exclusions_ms.gpkg"), exclusions_ms)
    GDF.write(joinpath(output_path, "exclusions_tenders.gpkg"), exclusions_tenders)

    return exclusions_ms, exclusions_tenders
end

function write_exclusions(
    ms_exclusions::DataFrame,
    t_exclusions::DataFrame,
    draft_ms::Float64,
    draft_t::Float64,
    subset_path::String,
    output_dir::String="outputs",
)
    !isdir(output_dir) && mkpath(output_dir)
    case_study_name = splitext(basename(subset_path))[1]
    if !isfile(joinpath(output_dir, "ms_exclusions_$(draft_ms)m_$(case_study_name).gpkg"))
        GDF.write(
            joinpath(output_dir, "ms_exclusions_$(draft_ms)m_$(case_study_name).gpkg"),
            ms_exclusions
        )
    end
    if !isfile(joinpath(output_dir, "t_exclusions_$(draft_t)m_$(case_study_name).gpkg"))
        GDF.write(
            joinpath(output_dir, "t_exclusions_$(draft_t)m_$(case_study_name).gpkg"),
            t_exclusions
        )
    end
end

"""
    export_mothership_routes(
    line_strings::Vector{LineString{2,Float64}},
    output_path::String="outputs"
)::DataFrame

Export mothership routes to a GeoPackage file, saved in the output directory.

# Arguments
- `line_strings`: LineStrings.
- `output_path`: Path to the output directory.

# Returns
- `df`: DataFrame with id and geometry columns.
"""
function export_mothership_routes(
    line_strings::Vector{LineString{2,Float64}},
    output_path::String="outputs"
)::DataFrame
    ids::Vector{Int64} = collect(1:length(line_strings))
    geometries::Vector{IGeometry{wkbLineString}} = _process_line.(line_strings)

    df = DataFrame(id=ids, geometry=geometries)

    GDF.write(joinpath(output_path, "routes_ms.gpkg"), df)
    return df
end

"""
    _process_line(line::LineString{2,Float64})::IGeometry{wkbLineString}

Process a LineString to create a geometry.

# Arguments
- `line`: LineString.

# Returns
- LineString geometry.
"""
function _process_line(line::LineString{2,Float64})::IGeometry{wkbLineString}
    return AG.createlinestring([(p[1], p[2]) for p in line.points])
end

"""
    export_tender_routes(
    tender_soln::Vector{TenderSolution},
    output_path::String="outputs"
)::DataFrame

Export tender routes to a GeoPackage file, saved in the output directory.

# Arguments
- `tender_soln`: Tender solutions.
- `output_path`: Path to the output directory.

# Returns
- `df`: DataFrame with id, cluster_id, sortie_id, and geometry columns.
"""
function export_tender_routes(
    tender_soln::Vector{TenderSolution},
    output_path::String="outputs"
)::DataFrame
    total_routes::Int64 = sum(length.(getfield.(tender_soln, :sorties)))
    ids::Vector{Int64} = collect(1:total_routes)

    cluster_ids::Vector{Int} = vcat(map(t -> fill(t.id, length(t.sorties)), tender_soln)...)
    sortie_ids::Vector{Int} = vcat(map(t -> 1:length(t.sorties), tender_soln)...)

    # Build geometry column as a linestring for each sortie in each cluster
    coords::Vector{Vector{NTuple{2,Float64}}} = [
        vcat(
            [(node[1], node[2]) for line in sortie.line_strings for node in line.points],
            [(tender.finish[1], tender.finish[2])]
        )
        for tender in tender_soln for sortie in tender.sorties
    ]
    geometries::Vector{IGeometry{wkbLineString}} = AG.createlinestring.(coords)

    df = DataFrame(
        id=ids,
        cluster_id=cluster_ids,
        sortie_id=sortie_ids,
        geometry=geometries
    )

    GDF.write(joinpath(output_path, "routes_tenders.gpkg"), df)

    return df
end
