const TARGET_SCENARIO_PATH = "data/test_target_locations.geojson"
const SUBSET_PATH = "data/test_target_area.gpkg"
const BATHY_PATH = "data/test_bathy.tif"
const WAVE_DISTURBANCE_PATH = "data/test_env_disturbances.geojson"

"""
Function returns problem instance and initial solution
"""
function initialize_instances(
    depot=(146.175, -16.84),
    draft_ms=Float64(-10.0),
    draft_t=Float64(-5.0),
    weight_ms=Float16(5.0),
    weight_t=Float16(2.0),
    n_tenders=Int8(3),
    t_cap=Int16(2),
)
    test_problem = HierarchicalRouting.load_problem(
        TARGET_SCENARIO_PATH,
        SUBSET_PATH,
        BATHY_PATH,
        WAVE_DISTURBANCE_PATH,
        depot,
        draft_ms,
        draft_t,
        weight_ms,
        weight_t,
        n_tenders,
        t_cap;
        debug_mode=true
    )
    return test_problem, initial_solution(test_problem; k=8)
end
