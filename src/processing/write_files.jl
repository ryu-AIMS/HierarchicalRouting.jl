
"""
    _get_output_details()::Tuple{String, Int}

Get the output path and EPSG code from the config file.

# Returns
- `output_path`: Path to the output directory.
- `EPSG_code`: EPSG code for the output coordinate reference system.
"""
function _get_output_details()::Tuple{String, Int}
    config = TOML.parsefile(joinpath("src",".config.toml"))

    output_path = config["output_dir"]["path"]
    EPSG_code = config["parameters"]["EPSG_code"]
    return output_path, EPSG_code
end

"""
    export_points(clusters::Vector{Cluster})::DataFrame

Export points to a GeoPackage file.
Saved in the output directory, using the EPSG code from the config file.

# Arguments
- `clusters`: Clusters.

# Returns
- `df`: DataFrame with id, geometry, and cluster_id columns.
"""
function export_points(clusters::Vector{Cluster})::DataFrame
    total_points = sum(length.(clusters.nodes))

    ids = Vector{Int}(undef, total_points)
    geoms = Vector{NTuple{2, Float64}}(undef, total_points)
    cluster_ids = Vector{Int}(undef, total_points)

    i = 1
    for (cluster_id, cluster) in enumerate(clusters)
        for point in cluster.nodes
            ids[i] = i
            geoms[i] = (point[1], point[2])
            cluster_ids[i] = cluster_id
            i += 1
        end
    end

    df = DataFrame(id = ids, geometry = geoms, cluster_id = cluster_ids)

    output_path, EPSG_code = _get_output_details()
    # TODO: crs based on lat/lons - but current coords as row/col refs
    GDF.write(joinpath(output_path, "points.gpkg"), df, crs=EPSG(EPSG_code))

    return df
end

"""
    export_clusters(cluster_sequence::DataFrame)::DataFrame

Export clusters to a GeoPackage file.
Saved in the output directory, using the EPSG code from the config file.

# Arguments
- `cluster_sequence`: Cluster sequence.

# Returns
- `df`: Cluster sequence with id, geometry, and order_id columns.
"""
function export_clusters(cluster_sequence::DataFrame)::DataFrame
    df = DataFrame(
        id = cluster_sequence.id,
        geometry = [AG.createpoint(lon, lat) for (lon, lat) in zip(cluster_sequence.lon, cluster_sequence.lat)],
        order_id = 0:size(cluster_sequence, 1)-1
    )

    output_path, EPSG_code = _get_output_details()

    # TODO: crs based on lat/lons - but current coords as row/col refs
    GDF.write(joinpath(output_path, "clusters.gpkg"), df, crs=EPSG(EPSG_code))

    return df
end

"""
    export_exclusions(exclusions::DataFrame, filepath::String)::Tuple{DataFrame, DataFrame}

Export exclusions to a GeoPackage file.
Saved in the output directory, using the EPSG code from the config file.

# Arguments
- `problem`: Problem instance.

# Returns
- `exclusions_ms`: Mothership exclusions.
- `exclusions_tenders`: Tenders exclusions.
"""
function export_exclusions(problem::Problem)::Tuple{DataFrame, DataFrame}
    exclusions_ms = problem.mothership.exclusion
    exclusions_tenders = problem.tenders.exclusion

    # Ensure DataFrames have primary key column (id)
    exclusions_ms.id = collect(1:nrow(exclusions_ms))
    exclusions_tenders.id = collect(1:nrow(exclusions_tenders))

    output_path, EPSG_code = _get_output_details()

    # TODO: crs based on lat/lons - but current coords as row/col refs
    GDF.write(joinpath(output_path, "exclusions_ms.gpkg"), exclusions_ms, crs=EPSG(EPSG_code))
    GDF.write(joinpath(output_path, "exclusions_tenders.gpkg"), exclusions_tenders, crs=EPSG(EPSG_code))

    return exclusions_ms, exclusions_tenders
end

"""
    export_mothership_routes(line_strings::Vector{LineString{2, Float64}})::DataFrame

Export mothership routes to a GeoPackage file.
Saved in the output directory, using the EPSG code from the config file.

# Arguments
- `line_strings`: LineStrings.

# Returns
- `df`: DataFrame with id and geometry columns.
"""
function export_mothership_routes(line_strings::Vector{LineString{2, Float64}})::DataFrame
    ids::Vector{Int64} = collect(1:length(line_strings))
    geometries::Vector{AG.IGeometry{AG.wkbLineString}} = process_line.(line_strings)

    df = DataFrame(id = ids, geometry = geometries)

    output_path, EPSG_code = _get_output_details()
    GDF.write(joinpath(output_path, "routes_ms.gpkg"), df, crs=EPSG(EPSG_code))

    return df
end

"""
    process_line(line::LineString{2,Float64})::AG.IGeometry{AG.wkbLineString}

Process a LineString to create a geometry.

# Arguments
- `line`: LineString.

# Returns
- LineString geometry.
"""
function process_line(line::LineString{2,Float64})::AG.IGeometry{AG.wkbLineString}
    coords = [
        ((p[1][1], p[1][2]), (p[2][1], p[2][2]))
        for p in line.points
    ]

    return AG.createlinestring(collect(Iterators.flatten(coords)))
end

"""
    export_tender_routes(tender_soln::Vector{TenderSolution})::DataFrame

Export tender routes to a GeoPackage file.
Saved in the output directory, using the EPSG code from the config file.

# Arguments
- `tender_soln`: Tender solutions.

# Returns
- `df`: DataFrame with id, cluster_id, sortie_id, and geometry columns.
"""
function export_tender_routes(tender_soln::Vector{TenderSolution})::DataFrame

    total_routes::Int64 = sum(length.(tender_soln.line_strings))
    ids::vector{Int64} = collect(1:total_routes)

    # Build cluster_id and sortie_id columns using nested comprehensions
    cluster_ids::Vector{Int} = [
        tender.id
        for tender in tender_soln
            for _ in tender.line_strings
    ]
    sortie_ids::Vector{Int64} = [
        sortie_id
        for tender in tender_soln
            for sortie_id in range(1:length(tender.line_strings))
    ]

    # Build geometry column: for each tender and sortie, create the linestring
    geometries::Vector{AG.IGeometry{AG.wkbLineString}} = [
        AG.createlinestring(
            vcat(
                [(node[1][1], node[1][2]) for segment in sortie_lines for node in segment],
                [(tender.finish[1], tender.finish[2])]
            )
        )
        for tender in tender_soln for sortie_lines in tender.line_strings
    ]

    df = DataFrame(
        id = ids,
        cluster_id = cluster_ids,
        sortie_id = sortie_ids,
        geometry = geometries
    )

    output_path, EPSG_code = _get_output_details()
    GDF.write(joinpath(output_path, "routes_tenders.gpkg"), df, crs=EPSG(EPSG_code))

    return df
end
