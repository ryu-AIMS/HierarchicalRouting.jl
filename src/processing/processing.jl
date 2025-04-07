"""
    target_threshold(targets::Raster, threshold::Float64)::Raster

Mask the target raster based on the specified threshold value.
i.e. set all values below the threshold to 0.

# Arguments
- `targets::Raster`: The target raster.
- `threshold::Float64`: The threshold value.

# Returns
- A masked raster.
"""
function target_threshold(targets::Raster, threshold::Float64)::Raster
    # TODO: add min/max argument flag, or use a range
    # TODO: include masks for other environmental constraints
    masked_raster = targets .* (targets .>= threshold)
    return masked_raster
end

function assign_wave_data(
    nodes::Vector{Point{2, Float64}},
    disturbance_data::DataFrame
)::Vector{Float64}
    disturbance_values = Vector{Union{Missing, Float64}}(undef, length(nodes))
    wave_val::Float64 = 0.0

    for (i, pt) in enumerate(nodes)
        wave_val = 0.0

        for row in eachrow(disturbance_data)
            if HierarchicalRouting.point_in_exclusion(pt, row.geometry)
                wave_val = row.Hs_MEAN
                break
            end
        end
        disturbance_values[i] = wave_val
    end

    return disturbance_values
end

