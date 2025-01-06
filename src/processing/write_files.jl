
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
- `clusters::Vector{HierarchicalRouting.Cluster}`: Clusters.
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

    return clusters
end
