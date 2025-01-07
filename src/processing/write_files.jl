
"""
    _get_output_details()

Get the output path and EPSG code from the config file.

# Returns
- `output_path::String`: Path to the output directory.
- `EPSG_code::Int`: EPSG code for the output coordinate reference system.
"""
function _get_output_details()
    config = TOML.parsefile(joinpath("src",".config.toml"))

    output_path = config["output_dir"]["path"]
    EPSG_code = config["parameters"]["EPSG_code"]
    return output_path, EPSG_code
end

"""
    export_points(clusters::Vector{HierarchicalRouting.Cluster})

Export points to a GeoPackage file.
Saved in the output directory, using the EPSG code from the config file.

# Arguments
- `clusters::Vector{HierarchicalRouting.Cluster}`: Clusters.

# Returns
- `df` : DataFrame with id, geometry, and cluster_id columns.
"""
function export_points(clusters::Vector{HierarchicalRouting.Cluster})
    df = DataFrame(
        id = Int[],
        geometry = NTuple{2, Float64}[],
        cluster_id = Int[]
    )

    for (cluster_id, cluster) in enumerate(clusters)
        for point in cluster.nodes
            push!(df, (size(df, 1) + 1, (point[1], point[2]), cluster_id))
        end
    end
    output_path, EPSG_code = _get_output_details()

    # TODO: crs based on lat/lons - but current coords as row/col refs
    GDF.write(joinpath(output_path, "points.gpkg"), df, crs=EPSG(EPSG_code))

    return df
end

"""
    export_clusters(cluster_sequence::DataFrame)

Export clusters to a GeoPackage file.
Saved in the output directory, using the EPSG code from the config file.

# Arguments
- `cluster_sequence::DataFrame`: Cluster sequence.

# Returns
- `df::DataFrame`: Cluster sequence with id, geometry, and order_id columns.
"""
function export_clusters(cluster_sequence::DataFrame) #, filepath::String)
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
    export_exclusions(exclusions::DataFrame, filepath::String)

Export exclusions to a GeoPackage file.
Saved in the output directory, using the EPSG code from the config file.

# Arguments
- `problem::HierarchicalRouting.Problem`: Problem instance.

# Returns
- `exclusions_ms::DataFrame`: Mothership exclusions.
- `exclusions_tenders::DataFrame`: Tenders exclusions.
"""
function export_exclusions(problem::HierarchicalRouting.Problem)
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
    export_mothership_routes(line_strings::Vector{LineString{2, Float64}})

Export mothership routes to a GeoPackage file.
Saved in the output directory, using the EPSG code from the config file.

# Arguments
- `line_strings::Vector{LineString{2, Float64}}`: LineStrings.

# Returns
- `df::DataFrame`: DataFrame with id and geometry columns.
"""
function export_mothership_routes(line_strings::Vector{LineString{2, Float64}})
    df = DataFrame(
        id = Int[],
        geometry = ArchGDAL.IGeometry{ArchGDAL.wkbLineString}[]
    )

    for (route_id, line) in enumerate(line_strings)
        coords = [((p[1][1], p[1][2]), (p[2][1], p[2][2])) for p in line.points]

        # Convert to ArchGDAL geometry
        flat_coords = Iterators.flatten(coords) |> collect
        # Pass flat coordinates
        line = AG.createlinestring(flat_coords)

        push!(df, (route_id, line))
    end

    output_path, EPSG_code = _get_output_details()

    GDF.write(joinpath(output_path, "routes_ms.gpkg"), df, crs=EPSG(EPSG_code))

    return df
end

"""
    export_tender_routes(tender_soln::Vector{HierarchicalRouting.TenderSolution})

Export tender routes to a GeoPackage file.
Saved in the output directory, using the EPSG code from the config file.

# Arguments
- `tender_soln::Vector{HierarchicalRouting.TenderSolution}`: Tender solutions.

# Returns
- `df::DataFrame`: DataFrame with id, cluster_id, sortie_id, and geometry columns.
"""
function export_tender_routes(tender_soln::Vector{HierarchicalRouting.TenderSolution})
    df = DataFrame(
        id = Int[],
        cluster_id = Int[],
        sortie_id = Int[],
        geometry = ArchGDAL.IGeometry{ArchGDAL.wkbLineString}[]
    )

    id = 1
    for tender in tender_soln
        for (sortie_id, sortie_lines) in enumerate(tender.line_strings)
            nodes = [
                [
                    (node[1][1],node[1][2])
                    for segment in sortie_lines
                    for node in segment
                ];
                (tender.finish[1], tender.finish[2])
            ]

            line_geometry = AG.createlinestring(nodes)
            push!(df, (id, tender.id, sortie_id, line_geometry))
            id += 1
        end
    end

    output_path, EPSG_code = _get_output_details()

    GDF.write(joinpath(output_path, "routes_tenders.gpkg"), df, crs=EPSG(EPSG_code))

    return df
end
