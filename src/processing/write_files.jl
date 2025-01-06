
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
