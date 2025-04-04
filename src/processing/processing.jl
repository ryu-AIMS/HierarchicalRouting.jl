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
