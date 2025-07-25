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
    get_disturbance_value(pt::Point{2,Float64}, disturbance_data::DataFrame)::Float64

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
        nodes::Vector{Point{2,Float64}},
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

"""
    linestring_segment_to_keep(
        section::Symbol,
        point::Point{2,Float64},
        line_strings::Vector{LineString{2,Float64}},
    )::Vector{LineString{2,Float64}}

Return the segment of line strings that contains the specified point, either from the start
or to the end of the line strings.

# Arguments
- `section`: A symbol indicating whether to keep the segment from the start (`:from`) or to the end (`:to`).
- `point`: The point to check against the line strings.
- `line_strings`: A vector of line strings to search through.

# Returns
- A vector of line strings that contain the specified point in the specified section.
"""
function linestring_segment_to_keep(
    section::Symbol,
    point::Point{2,Float64},
    line_strings::Vector{LineString{2,Float64}},
)::Vector{LineString{2,Float64}}
    linestring_points = getfield.(line_strings, :points)
    if section == :from
        start_seg_points = getindex.(linestring_points, 1)
        leg_start_idx = findfirst(==(point), start_seg_points)
        segment = line_strings[leg_start_idx:end]
        return segment
    elseif section == :to
        end_seg_points = getindex.(linestring_points, 2)
        leg_end_idx = findfirst(==(point), end_seg_points)
        segment = line_strings[1:leg_end_idx]
        return segment
    else
        error("Invalid section specified. Use `:from` or `:to`")
    end
end

"""
    make_superdiag_matrix(v::Vector{Float64})::Matrix{Float64}

Make a square distance matrix from a vector of distances. The matrix will have zeros on the
diagonal and distance values on the first superdiagonal.
- Superdiagonal is the first diagonal above the main diagonal.
- **NOTE**: The distance matrix is filled with **Inf** other than the diagonal and superdiagonal.

"""
function make_superdiag_matrix(v::Vector{Float64})::Matrix{Float64}
    n = length(v) + 1
    M = fill(Inf, n, n)
    M[CartesianIndex.(1:n, 1:n)] .= 0.0
    M[CartesianIndex.(1:n-1, 2:n)] = v
    return M
end

"""
    get_superdiag_vals(M::Matrix{Float64})::Vector{Float64}

Return the values from the first superdiagonal of a square matrix.
- Superdiagonal is the first diagonal above the main diagonal.
"""
function get_superdiag_vals(M::Matrix{Float64})::Vector{Float64}
    return [M[i, i+1] for i in 1:(size(M, 1)-1)]
end
