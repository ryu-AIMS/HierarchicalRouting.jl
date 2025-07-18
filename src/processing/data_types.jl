"""
    polygonize_binary(rast::Raster)::DataFrame

Polygonize a bit/binary Raster dataset, returning only the `true` regions.

Implementation based on:

- https://github.com/jl-pkgs/GeoArrays.jl/blob/master/src/gdal_polygonize.jl#L3
- Reading the examples and code in [GDAL.jl](https://github.com/JuliaGeo/GDAL.jl)

# Arguments
- `rast` : Raster

# Returns
DataFrame of geometries
"""
function polygonize_binary(rast::Raster)::DataFrame
    # Convert to ArchGDAL dataset
    w, h = size(rast)
    d = AG.create(AG.getdriver("MEM");
        width=w,
        height=h,
        nbands=1,
        dtype=Int8
    )

    # Apply affine transformation
    x, y = dims(rast)
    y0 = maximum([y[1], y[end]])
    transform = [
        x[1],       # x0: left edge in CRS units
        step(x),    # dx: pixel width in CRS units
        0.0,        # row rotation (typically 0)
        y0,         # y0: top edge in CRS units
        0.0,        # column rotation (typically 0)
        step(y)     # dy: pixel height in CRS units (negative = north-up orientation)
    ]
    AG.setgeotransform!(d, transform)

    # Create raster dataset and write data
    dataset = AG.RasterDataset(d)
    band = AG.getband(dataset, 1)
    AG.write!(band, convert.(Int8, rast.data))

    # Set projection
    crs = Rasters.crs(rast)
    AG.setproj!(dataset, crs.val)

    # Create in-memory store for geometries
    drv = AG.GDAL.gdalgetdriverbyname("Memory")
    geom_ds = AG.GDAL.gdalcreate(drv, "", 0, 0, 0, AG.GDAL.GDT_Unknown, C_NULL)

    # Define CRS
    gdal_ref = AG.GDAL.osrnewspatialreference(C_NULL)
    AG.GDAL.osrimportfromepsg(gdal_ref, convert(EPSG, crs).val[1])

    # Create layer
    layer = AG.GDAL.gdaldatasetcreatelayer(
        geom_ds,
        "out",
        gdal_ref,
        AG.GDAL.wkbPolygon,
        C_NULL
    )

    # Create field and store polygonized data
    field_defn = AG.GDAL.ogr_fld_create("exclusion", AG.GDAL.OFTInteger)
    field = AG.GDAL.ogr_l_createfield(layer, field_defn, AG.GDAL.TRUE)
    AG.GDAL.gdalpolygonize(
        band,
        C_NULL,
        layer,
        field,
        C_NULL,
        C_NULL,
        C_NULL
    )

    # Collect geometries
    # Can't easily pre-allocate the vector as we don't know how many true-valued geoms
    # there are.
    n_features = AG.GDAL.ogr_l_getfeaturecount(layer, 1)
    layer_geoms = Vector{IGeometry}()
    sizehint!(layer_geoms, n_features)
    for i in 0:(n_features-1)
        feat = AG.Feature(AG.GDAL.ogr_l_getfeature(layer, i))
        field_val = AG.getfield(feat, 0)  # Get first field (C is a zero-indexed language)

        if field_val == 1
            push!(layer_geoms, AG.getgeom(feat))
        end

        AG.destroy(feat)
    end

    # Close handles
    AG.destroy(band)
    AG.destroy(d)
    AG.GDAL.gdalclose(geom_ds)
    AG.GDAL.ogr_fld_destroy(field_defn)
    AG.GDAL.osrrelease(gdal_ref)

    return DataFrame(geometry=layer_geoms)
end
