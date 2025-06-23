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
function get_disturbance_value(pt::Point{2,Float64}, disturbance_data::DataFrame)::Float64
    for row in eachrow(disturbance_data)
        if HierarchicalRouting.point_in_exclusion(pt, row.geometry)
            return row.Hs_MEAN
        end
    end
    return 0.0 # Default value if no data
end

"""
    create_disturbance_data_dataframe(
        nodes::Vector{Point{2, Float64}},
        disturbance_data::DataFrame
    )::DataFrame

Return a DataFrame containing nodes and their corresponding disturbance values.

# Arguments
- `nodes`: A vector of nodes (points).
- `disturbance_data`: The disturbance data containing disturbance values.

# Returns
- A DataFrame with two columns: `node` and `disturbance_value`.
"""
function create_disturbance_data_dataframe(
    nodes::Vector{Point{2,Float64}},
    disturbance_data::DataFrame
)::DataFrame
    disturbance_values = get_disturbance_value.(nodes, Ref(disturbance_data))
    return DataFrame(node=nodes, disturbance_value=disturbance_values)
end
