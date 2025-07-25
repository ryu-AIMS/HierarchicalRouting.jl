
using Glob

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
    depot::Point{2,Float64}
    targets::Targets
    mothership::Vessel
    tenders::Vessel
end

"""
    load_problem(
        target_path::AbstractString,
        subset_path::AbstractString,
        env_data_path::AbstractString,
        env_disturbance_path::AbstractString,
        depot::NTuple{2,Float64},
        draft_ms::Real,
        draft_t::Real,
        weight_ms::Real,
        weight_t::Real,
        n_tenders::Integer,
        t_cap::Integer;
        target_subset_path::AbstractString="",
        output_dir::AbstractString="outputs/",
        debug_mode::Bool=false,
    )::Problem

Load the problem data to create a `Problem` object from given parameters.

# Arguments
- `target_path`: The path of the target scenario to load, in GeoJSON format.
- `subset_path` : The path to the subset GeoPackage file.
- `env_data_path`: The path to the environmental data directory.
- `env_disturbance_path`: The path to the environmental disturbance GeoDataFrame data.
- `depot`: A Point representing the depot location, in the form (lon, lat) which is the used
    as the coincident start and finish point for the deployment plan.
- `draft_ms`: The draft of the mothership (in metres), given as a negative value, indicating
    depth of the vessel below the surface.
- `draft_t`: The draft of the tenders (in metres), given as a negative value, indicating
    depth of the vessel below the surface.
- `weight_ms`: The weighting factor for the mothership, as a multiplier of distance travelled to
    quantify cost.
- `weight_t`: The weighting factor for the tenders, as a multiplier of distance travelled to
    quantify cost.
- `n_tenders`: The number of tenders available to use for deployment.
- `t_cap`: The capacity of the tenders, determining the maximum number of targets that can
    be visited/deployed by a tender in a single trip.
- `target_subset_path`: The path to save the target subset raster. Default is "".
- `output_dir`: The directory to save output files. Default is "outputs/".
- `debug_mode`: A boolean flag to enable debug mode. Default = false.
     - If `debug_mode=true`, the function will not write to the output directory.

# Returns
The problem object.
"""
function load_problem(
    target_path::AbstractString,
    subset_path::AbstractString,
    env_data_path::AbstractString,
    env_disturbance_path::AbstractString,
    depot::NTuple{2,Float64},
    draft_ms::Real,
    draft_t::Real,
    weight_ms::Real,
    weight_t::Real,
    n_tenders::Integer,
    t_cap::Integer;
    target_subset_path::AbstractString="",
    output_dir::AbstractString="outputs/",
    debug_mode::Bool=false,
)::Problem
    depot = Point{2,Float64}(depot)
    draft_ms = Float64(draft_ms)
    draft_t = Float64(draft_t)
    weight_ms = Float16(weight_ms)
    weight_t = Float16(weight_t)
    n_tenders = Int8(n_tenders)
    t_cap = Int16(t_cap)
    subset = GDF.read(subset_path)
    subset_bbox = get_bbox_bounds_from_df(subset)
    target_gdf_subset = filter_within_bbox(GDF.read(target_path), subset_bbox)

    target_raster = process_targets(target_gdf_subset.geometry, subset, target_subset_path)
    targets_gdf = raster_to_gdf(target_raster)

    env_disturbance = GDF.read(env_disturbance_path)
    disturbance_data_subset = filter_within_bbox(env_disturbance, subset_bbox)
    suitable_targets_subset = process_geometry_targets(target_gdf_subset.geometry)
    indices::Vector{CartesianIndex{2}} = findall(
        x -> x != suitable_targets_subset.missingval,
        suitable_targets_subset
    )
    coords::Vector{Point{2,Float64}} = [
        Point(
            suitable_targets_subset.dims[1][idx[1]],
            suitable_targets_subset.dims[2][idx[2]]
        )
        for idx in indices
    ]
    disturbance_df = create_disturbance_data_dataframe(coords, disturbance_data_subset)
    targets = Targets(targets_gdf, target_path, disturbance_df)

    # Process exclusions
    if !(debug_mode)
        if !isempty(output_dir)
            mkpath(output_dir)
        end
        ms_exclusion_zones_df = read_and_polygonize_exclusions(
            env_data_path,
            draft_ms,
            subset,
            "ms_exclusion",
            output_dir,
        )
        t_exclusions = read_and_polygonize_exclusions(
            env_data_path,
            draft_t,
            subset,
            "t_exclusion",
            output_dir,
        )
    else
        ms_exclusion_zones_df = read_and_polygonize_exclusions(
            env_data_path,
            draft_ms,
            subset,
        )
        t_exclusions = read_and_polygonize_exclusions(
            env_data_path,
            draft_t,
            subset,
        )
    end
    ms_exclusions::DataFrame = ms_exclusion_zones_df |> filter_and_simplify_exclusions! |>
                               buffer_exclusions! |> unionize_overlaps! |> filter_and_simplify_exclusions!

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
        exclusion=ms_exclusions,
        weighting=weight_ms
    )
    tenders = Vessel(
        exclusion=t_exclusions,
        capacity=t_cap,
        number=n_tenders,
        weighting=weight_t
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
    get_bbox_bounds_from_df(df::DataFrame)::NTuple{4,Float64}

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
function get_bbox_bounds_from_df(df::DataFrame)::NTuple{4,Float64}
    min_x = minimum(getfield.(AG.envelope.(df.geom), 1))
    max_x = maximum(getfield.(AG.envelope.(df.geom), 2))
    min_y = minimum(getfield.(AG.envelope.(df.geom), 3))
    max_y = maximum(getfield.(AG.envelope.(df.geom), 4))

    return (min_x, max_x, min_y, max_y)
end

"""
    filter_within_bbox(
        gdf::DataFrame,
        bbox_bounds::NTuple{4,Float64}
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
    bbox_bounds::NTuple{4,Float64}
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
    points::Vector{Point{2,Float64}},
    exclusions::DataFrame
)::DataFrame
    exclusion_geometries = exclusions.geometry
    poly_adjustment_mask = point_in_exclusion.(points, Ref(exclusion_geometries))
    contained_points = points[poly_adjustment_mask]
    containing_poly_ids = containing_exclusion.(contained_points, Ref(exclusion_geometries))

    # Dictionary to map exclusion polygons to contained points
    polygon_points_map = Dict{Int,Vector{Point{2,Float64}}}()
    for (point, poly_id) in zip(contained_points, containing_poly_ids)
        if !haskey(polygon_points_map, poly_id)
            polygon_points_map[poly_id] = Vector{Point{2,Float64}}()
        end
        push!(polygon_points_map[poly_id], point)
    end

    updated_vertices = Dict{Int,Vector{Tuple{Float64,Float64}}}()
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
