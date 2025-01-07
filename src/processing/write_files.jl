
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
    GDF.write(joinpath(output_path, "test_points.gpkg"), df, crs=EPSG(EPSG_code))

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
    GDF.write(joinpath(output_path, "test_clusters.gpkg"), df, crs=EPSG(EPSG_code))

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
    GDF.write(joinpath(output_path, "test_exclusions_ms.gpkg"), exclusions_ms, crs=EPSG(EPSG_code))
    GDF.write(joinpath(output_path, "test_exclusions_tenders.gpkg"), exclusions_tenders, crs=EPSG(EPSG_code))

    return exclusions_ms, exclusions_tenders
end
