
"""
    crop_to_subset(spatial_dataset::Raster, subset::DataFrame)

Extract a subset from the given `spatial_dataset` raster based on the specified `subset`.

# Arguments
- `spatial_dataset::Raster`: The spatial dataset from which to extract the subset.
- `subset::DataFrame`: The subset DataFrame containing geometries in area.

# Returns
- A cropped raster based on the subset area.
"""
function crop_to_subset(spatial_dataset::Raster, subset::DataFrame)
    return Rasters.crop(spatial_dataset, to=subset.geom)
end

"""
    target_threshold(targets::Raster, threshold::Float64)

Mask the target raster based on the specified threshold value.
i.e. set all values below the threshold to 0.

# Arguments
- `targets::Raster`: The target raster.
- `threshold::Float64`: The threshold value.

# Returns
- A masked raster.
"""
function target_threshold(targets::Raster, threshold::Float64)
    # TODO: add min/max argument flag, or use a range
    # TODO: include masks for other environmental constraints
    masked_raster = targets .* (targets .>= threshold)
    return masked_raster
end

"""
    get_linestrings(graph_matrix::Matrix{
    Tuple{
        Vector{Point{2, Float64}},
        Vector{SimpleWeightedGraphs.SimpleWeightedEdge{Int64, Float64}}
    }},
    waypoints::Vector{Point{2, Float64}})

Between each sequential pair of waypoints: create a single LineString across all points defining its path.

# Arguments
- `graph_matrix::Matrix{Tuple{Dict{Int64, Point{2, Float64}}, Vector{SimpleWeightedGraphs.SimpleWeightedEdge{Int64, Float64}}}}`: The graph matrix.
- `waypoints::Vector{Point{2, Float64}}`: The waypoints.

# Returns
- A vector of LineStrings.
"""
function get_linestrings(graph_matrix::Matrix{
    Tuple{
        Vector{Point{2, Float64}},
        Vector{SimpleWeightedGraphs.SimpleWeightedEdge{Int64, Float64}}
    }},
    waypoints::Vector{Point{2, Float64}}
)

    line_strings = Vector{LineString{2, Float64}}()

    for i in 1:(length(waypoints) - 1)
        start_point = waypoints[i]
        end_point = waypoints[i + 1]

        # Find the graph entry corresponding to this waypoint pair
        for graph_dict_path in graph_matrix
            point_vector, shortest_path = graph_dict_path

            # Check if the path starts and ends at the desired waypoints
            if !isempty(shortest_path) && (
                (
                    point_vector[shortest_path[1].src] == start_point &&
                    point_vector[shortest_path[end].dst] == end_point
                )
                ||
                (
                    point_vector[shortest_path[1].src] == end_point &&
                    point_vector[shortest_path[end].dst] == start_point
                )
            )

                # Extract points along path
                path_points = [point_vector[e.src] for e in shortest_path]
                push!(path_points, point_vector[shortest_path[end].dst])

                push!(line_strings, LineString(path_points))
                break
            end
        end
    end

    return line_strings
end
