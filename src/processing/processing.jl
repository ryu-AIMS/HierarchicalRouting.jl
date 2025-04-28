"""
    target_threshold(targets::Raster, threshold::Float64)::Raster

Mask the target raster based on the specified threshold value.
i.e. set all values below the threshold to 0.

# Arguments
- `targets`: The target raster.
- `threshold`: The threshold value.

# Returns
- A masked raster.
"""
function target_threshold(targets::Raster, threshold::Float64)::Raster
    # TODO: add min/max argument flag, or use a range
    # TODO: include masks for other environmental constraints
    masked_raster = targets .* (targets .>= threshold)
    return masked_raster
end

"""
    get_disturbance_value(pt::Point{2, Float64}, disturbance_data::DataFrame)::Float64

Return the value for a given point based on environmental disturbance data.

# Arguments
- `pt`: The location point for which to get the environmental data value.
- `disturbance_data`: The disturbance data containing environmental values. Given in the
    form of a DataFrame with columns `geometry` and `Hs_MEAN`.

# Returns
- The environmental data value for the point.
"""
function get_disturbance_value(pt::Point{2, Float64}, disturbance_data::DataFrame)::Float64
    for row in eachrow(disturbance_data)
        if HierarchicalRouting.point_in_exclusion(pt, row.geometry)
            return row.Hs_MEAN
        end
    end
    return 0.0 # Default value if no data
end

"""
    assign_disturbance_data(
        nodes::Vector{Point{2, Float64}},
        disturbance_data::DataFrame
    )::Vector{Float64}

Assign disturbance data to a vector of nodes.

# Arguments
- `nodes`: A vector of nodes (points).
- `disturbance_data`: The disturbance data containing wave values.

# Returns
- A vector of wave values corresponding to the nodes.
"""
function assign_disturbance_data(
    nodes::Vector{Point{2, Float64}},
    disturbance_data::DataFrame
)::Vector{Float64}
    return [get_disturbance_value(pt, disturbance_data) for pt in nodes]
end

"""
    create_disturbance_data_dataframe(
        nodes::Vector{Point{2, Float64}},
        disturbance_data::DataFrame
    )::DataFrame

Return a DataFrame containing nodes and their corresponding disturbance values.

# Arguments
- `nodes`: A vector of nodes (points).
- `disturbance_data`: The disturbance data containing wave values.

# Returns
- A DataFrame with two columns: `node` and `wave_value`.
"""
function create_disturbance_data_dataframe(
    nodes::Vector{Point{2, Float64}},
    disturbance_data::DataFrame
)::DataFrame
    wave_values = [get_disturbance_value(pt, disturbance_data) for pt in nodes]
    return DataFrame(node = nodes, wave_value = wave_values)
end
