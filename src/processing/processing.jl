
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
    target_threshold(target_bathy::Raster, ms_depth)

Mask the target environmental raster based on the minimum threshold value.

# Arguments
- `target_bathy::Raster`: The target bathymetry raster.
- `ms_depth`: The minimum depth threshold.
"""
function target_threshold(targets::Raster, threshold::Float64)
    # TODO: add min/max argument flag, or use a range
    # TODO: include masks for other environmental constraints
    masked_raster = targets .* (targets .>= threshold)
    return masked_raster
end
