import ArchGDAL as AG
using DataFrames
using Statistics

"""
    extract_subset(spatial_dataset::Raster, subset::GeoDataFrame)

Extract a subset of a raster dataset based on a GeoDataFrame.
"""
function extract_subset(spatial_dataset::Raster, subset)
    result_raster = Rasters.trim(mask(spatial_dataset; with=subset.geom))
    return result_raster
end
