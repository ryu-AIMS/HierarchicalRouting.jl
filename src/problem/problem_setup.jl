
using Glob
using TOML

# struct Threshold
#     # value::Float64
#     # type::Symbol

#     min::Float64
#     max::Float64

#     Threshold(; min::Float64=-Inf, max::Float64=Inf) = new(min, max)
# end

# struct Constraint
#     name::String
#     threshold::Threshold
#     file_path::String
# end

@kwdef struct Vessel
    exclusion::DataFrame
    capacity::Int16 = typemax(Int16)
    number::Int8 = Int8(1)
    weighting::Float16 = Float16(1.0)

    function Vessel(exclusion::DataFrame, capacity::Int16, number::Int8, weighting::Float16)
        errors = validate_vessel(capacity, number, weighting)
        if isempty(errors)
            return new(exclusion, capacity, number, weighting)
        else
            throw(ArgumentError(join(errors, "\n")))
        end
    end
end

function validate_vessel(capacity::Int16, number::Int8, weighting::Float16)::Vector{String}
    errors::Vector{String} = []
    (weighting <= 0.0) ? push!(errors, "Weighting must be greater than 0") : 0
    (capacity <= 0) ? push!(errors, "Capacity must be greater than 0") : 0
    (number <= 0) ? push!(errors, "Number of vessels must be greater than 0") : 0
    return errors
end

struct Targets
    points::DataFrame
    path::String
    disturbance_gdf::DataFrame
end

struct Problem
    depot::Point{2, Float64}
    targets::Targets
    mothership::Vessel
    tenders::Vessel
end

"""
    load_problem(target_path::String)::Problem

Load the problem data from the configuration file and return a Problem object.

# Arguments
- `target_path`: The name of the target scenario to load.

# Returns
The problem object.
"""
function load_problem(target_path::String)::Problem
    config = TOML.parsefile(".config.toml")

    draft_ms::Float64 = config["parameters"]["depth_ms"]
    draft_t::Float64 = config["parameters"]["depth_t"]
    n_tenders::Int8 = config["parameters"]["n_tenders"]
    t_cap::Int16 = config["parameters"]["t_cap"]
    EPSG_code::Int16 = config["parameters"]["EPSG_code"]

    scenario_name = try
        split(split(split(target_path, "/")[end], "\\")[end], ".")[1]
    catch
        ""
    end

    depot::Point{2, Float64} = Point{2, Float64}(
        config["parameters"]["depot_x"], config["parameters"]["depot_y"]
    )

    wave_disturbance_dir = config["data_dir"]["wave_disturbances"]
    wave_disturbance = GDF.read(first(glob("*.geojson", wave_disturbance_dir)))

    env_constraints_dir = config["data_dir"]["env_constraints"]
    env_subfolders = readdir(env_constraints_dir)
    env_paths = Dict(
        subfolder => joinpath(env_constraints_dir, subfolder)
        for subfolder in env_subfolders
    )
    env_data = Dict(
        subfolder => (
            env_dir = path,
            rast_file = isempty(glob("*.tif", path)) ? nothing : first(glob("*.tif", path))
        ) for (subfolder, path) in env_paths
    )
    # TODO: Process all environmental constraints to create single cumulative exclusion zone

    site_dir = config["data_dir"]["site"]
    subset = GDF.read(first(glob("*.gpkg", site_dir)))
    subset_bbox = get_bbox_bounds_from_df(subset)
    target_gdf_subset = filter_within_bbox(
        GDF.read(target_path), subset_bbox
    )

    output_dir = config["output_dir"]["path"]
    suitable_targets_filename = splitext(basename(target_path))[1]
    target_subset_path::String = joinpath(
        output_dir,
        "target_subset_$(suitable_targets_filename).tif"
    )

    target_raster = process_targets(
        target_gdf_subset.geometry,
        target_subset_path,
        subset,
        EPSG_code
    )
    targets_gdf = raster_to_gdf(target_raster)

    disturbance_data_subset = filter_within_bbox(
        wave_disturbance, subset_bbox
    )
    suitable_targets_subset = process_geometry_targets(
        target_gdf_subset.geometry, EPSG_code
    )
    indices::Vector{CartesianIndex{2}} = findall(
        x -> x != suitable_targets_subset.missingval,
        suitable_targets_subset
    )
    coords::Vector{Point{2, Float64}} = [
        Point(
            suitable_targets_subset.dims[1][idx[1]],
            suitable_targets_subset.dims[2][idx[2]]
        )
        for idx in indices
    ]
    disturbance_df = create_disturbance_data_dataframe(
        coords,
        disturbance_data_subset
    )
    targets = Targets(targets_gdf, target_path, disturbance_df)

    bathy_subset_path = joinpath(output_dir, "bathy_subset.tif")

    # Process exclusions
    if !(config["DEBUG"]["debug_mode"])
        ms_exclusion_zones_df = read_and_polygonize_exclusions(
            env_data["bathy"].rast_file,
            draft_ms,
            subset,
            EPSG_code,
            bathy_subset_path,
            joinpath(output_dir, "ms_exclusion.gpkg"),
            joinpath(output_dir, "ms_exclusion.tif")
        )
    else
        ms_exclusion_zones_df = read_and_polygonize_exclusions(
            env_data["bathy"].rast_file,
            draft_ms,
            subset,
            EPSG_code
        )
    end
    ms_exclusions::DataFrame = ms_exclusion_zones_df |> filter_and_simplify_exclusions! |>
        buffer_exclusions! |> unionize_overlaps! |> filter_and_simplify_exclusions!

    if !(config["DEBUG"]["debug_mode"])
        t_exclusions = read_and_polygonize_exclusions(
            env_data["bathy"].rast_file,
            draft_t,
            subset,
            EPSG_code,
            bathy_subset_path,
            joinpath(output_dir, "t_exclusion.gpkg"),
            joinpath(output_dir, "t_exclusion.tif")
        )
    else
        t_exclusions = read_and_polygonize_exclusions(
            env_data["bathy"].rast_file,
            draft_t,
            subset,
            EPSG_code
        )
    end

    filter_and_simplify_exclusions!(t_exclusions, min_area=1E-7, simplify_tol=5E-4)
    buffer_exclusions!(t_exclusions, buffer_dist=0.0)
    unionize_overlaps!(t_exclusions)
    #? Redundant 2nd call to simplify_exclusions!()
    filter_and_simplify_exclusions!(t_exclusions, min_area=1E-7, simplify_tol=5E-4)
    t_exclusions = adjust_exclusions(
        targets.points.geometry,
        t_exclusions
    )

    mothership = Vessel(
        exclusion = ms_exclusions,
        weighting = Float16(config["parameters"]["weight_ms"]) #! user-defined
    )
    tenders = Vessel(
        exclusion = t_exclusions,
        capacity = t_cap,
        number = n_tenders,
        weighting = Float16(config["parameters"]["weight_t"]) #! user-defined
    )

    return Problem(depot, targets, mothership, tenders)
end

"""
    is_within_bbox(geom, min_x, max_x, min_y, max_y)::Bool

Checks if the geometry is within the bounding box defined by the given coordinates.

# Arguments
- `geom`: The geometry to check.
- `min_x`: The minimum x-coordinate of the bounding box.
- `max_x`: The maximum x-coordinate of the bounding box.
- `min_y`: The minimum y-coordinate of the bounding box.
- `max_y`: The maximum y-coordinate of the bounding box.

# Returns
True if the geometry is within the bounding box, false otherwise.
"""
function is_within_bbox(geom, min_x, max_x, min_y, max_y)::Bool
    env = AG.envelope(geom)
    return env.MinX ≥ min_x && env.MaxX ≤ max_x && env.MinY ≥ min_y && env.MaxY ≤ max_y
end

"""
    get_bbox_bounds_from_df(df::DataFrame)::Tuple{Float64, Float64, Float64, Float64}

Retrieves the bounding box bounds from the given GeoDataFrame.

# Arguments
- `df`: The GeoDataFrame to retrieve the bounding box from.

# Returns
A tuple containing the bounding box coordinates:
- min_x,
- max_x,
- min_y,
- max_y.
"""
function get_bbox_bounds_from_df(df::DataFrame)::Tuple{Float64, Float64, Float64, Float64}
    min_x = minimum(getfield.(AG.envelope.(df.geom), 1))
    max_x = maximum(getfield.(AG.envelope.(df.geom), 2))
    min_y = minimum(getfield.(AG.envelope.(df.geom), 3))
    max_y = maximum(getfield.(AG.envelope.(df.geom), 4))

    return (min_x, max_x, min_y, max_y)
end

"""
    filter_within_bbox(
        gdf::DataFrame,
        bbox_bounds::Tuple{Float64, Float64, Float64, Float64}
    )::DataFrame

Filters the geometries in the given GeoDataFrame `gdf` that are within the bounding box
    defined by the geometries in `subset`.

# Arguments
- `gdf`: The GeoDataFrame to filter.
- `bbox_bounds`: A tuple containing the bounding box coordinates:
    - min_x,
    - max_x,
    - min_y,
    - max_y.

# Returns
A GeoDataFrame containing only the geometries within the bounding box defined by `subset`.
"""
function filter_within_bbox(
    gdf::DataFrame,
    bbox_bounds::Tuple{Float64, Float64, Float64, Float64}
)::DataFrame
    filtered_gdf = filter(
        row -> is_within_bbox(
            row.geometry,
            bbox_bounds[1], bbox_bounds[2], bbox_bounds[3], bbox_bounds[4]
        ),
        gdf
    )
    return filtered_gdf
end

function adjust_exclusions(
    points::Vector{Point{2, Float64}},
    exclusions::DataFrame
)::DataFrame
    exclusion_geometries = exclusions.geometry
    poly_adjustment_mask = point_in_exclusion.(points, Ref(exclusions))
    contained_points = points[poly_adjustment_mask]
    containing_poly_ids = containing_exclusion.(contained_points, Ref(exclusions))

    # Dictionary to map exclusion polygons to contained points
    polygon_points_map = Dict{Int, Vector{Point{2, Float64}}}()
    for (point, poly_id) in zip(contained_points, containing_poly_ids)
        if !haskey(polygon_points_map, poly_id)
            polygon_points_map[poly_id] = Vector{Point{2, Float64}}()
        end
        push!(polygon_points_map[poly_id], point)
    end

    updated_vertices = Dict{Int, Vector{Tuple{Float64, Float64}}}()
    for (poly_id, points_to_insert) in polygon_points_map
        @info "Processing polygon ID: $poly_id"
        polygon = exclusion_geometries[poly_id]
        polygon_points = AG.getgeom(polygon, 0)
        n_vertices = AG.ngeom(polygon_points)
        polygon_vertices = [AG.getpoint(polygon_points, i)[1:2] for i in 0:n_vertices-1]

        for point in points_to_insert
            min_dist = Inf
            insert_index = n_vertices

            for i in 1:n_vertices-1
                p1 = Point{2,Float64}(polygon_vertices[i])
                p2 = Point{2,Float64}(polygon_vertices[i+1])

                # TODO: Use a more robust method to find the closest point on the edge
                dist = (
                    abs(p1[1] - point[1]) + abs(p1[2] - point[2])
                    +
                    abs(p2[1] - point[1]) + abs(p2[2] - point[2])
                )

                if dist < min_dist
                    min_dist = dist
                    insert_index = i
                end
            end
            # Insert the point at index
            pre_points = polygon_vertices[1:insert_index]
            post_points = polygon_vertices[insert_index+1:n_vertices]

            polygon_vertices = vcat(pre_points, [(point[1], point[2])], post_points)
            n_vertices = length(polygon_vertices)
        end
        updated_vertices[poly_id] = polygon_vertices
    end

    for i in keys(polygon_points_map)
        exclusions.geometry[i] = AG.createpolygon(updated_vertices[i])
    end

    return exclusions
end
