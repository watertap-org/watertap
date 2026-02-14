from parameter_sweep.loop_tool.loop_tool import loopTool, get_working_dir
from watertap.flowsheets.lsrro import lsrro
import multiprocessing

from watertap.core.solvers import get_solver
from watertap.tools.sensitivity_manager import SensitivityManager

##############################################################
# For loop tool/parameter sweeep interface we need in general, 3 functions:
# 1) build_function - builds the model and configures it for initialization
# 2) initialize_function - initializes the model and configures it for optimization
# 3) optimize_function - optimizes/solves the model
# 4) (optional) probe_function - probes the model for feasibility or other metrics (e.g. ensures we don't try to solve when recovery is too high)

# Your workflow should go as follows:
# 1) Build flowsheet/model and define all functions
# 2) ensure model actually solves using defined functions
# 3) set up loop tool with yaml file and run analysis
# 4) use psPlotKit to analyze results from h5 file - unlike in multi_sweep example, we don't need to
#    define expression or outputs manually, all data on flowsheet will be saved and we can do all the
#    analysis in psPlotKit
##############################################################


def build_function(
    number_of_stages=2,
    has_CP=True,
    add_sensitivity=False,
    target_tds=125,
    target_recovery=0.5,
):
    """Build LSRRO model with optinal arguments, this should only build the model and configure it
    for initialization
    Args:
        number_of_stages (int): number of stages in the LSRRO system
        A_value (float): membrane water permeability in m/(s*Pa)
        permeate_quality_limit (float): permeate quality limit in kg solute/kg solution
        has_CP (bool): whether to include concentration polarization calculations
    Returns:
        m (ConcreteModel): LSRRO model"""

    m = lsrro.build(
        number_of_stages=number_of_stages,
        has_NaCl_solubility_limit=True,
        has_calculated_concentration_polarization=has_CP,
        has_calculated_ro_pressure_drop=True,
    )
    lsrro.set_operating_conditions(m)

    # For use in case configuration analysis, we can store target tds and recovery on model
    m.target_tds = target_tds
    m.target_recovery = target_recovery
    if add_sensitivity:
        add_sensitivity_analysis(m)
    return m


def initialize_model(
    m,
    quick_start,
    A_value,
    permeate_quality_limit,
    feed_tds=35,
    recovery=0.45,
    **kwargs,
):
    """Function to initialize LSRRO model, this function should initialize the model
     and set it up for optimization or use in the sweep. (e.g The next function we run will solve the model!)

    Args:
        m (ConcreteModel): LSRRO model to be initialized
        quick_start (bool): whether to use quick start initialization
        A_value (float): membrane water permeability in m/(s*Pa)
        permeate_quality_limit (float): permeate quality limit in kg solute/kg solution
        **kwargs: additional keyword arguments (These are required for use with paramter sweep and loop tool!)
    Returns:
        None
    """
    if not quick_start:
        lsrro.initialize(m)
    lsrro.solve(m)
    m.fs.feed.flow_mass_phase_comp.unfix()
    m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].fix()
    m.fs.feed.properties[0].flow_vol_phase["Liq"].fix()
    lsrro.optimize_set_up(
        m,
        set_default_bounds_on_module_dimensions=True,
        A_value=A_value,
        permeate_quality_limit=permeate_quality_limit,
    )

    # if we have sense manager on model, activate all sense constraints.
    if m.fs.find_component("sense_manager"):
        # aggregate these pressure constraints that are built during optimize_set_up call, which should
        # really be done during "build call" so we can let sense-manager manage them for
        # our sensitivity analysis
        existing_pressure_constraints = []
        for pump in m.fs.LSRRO_Stages:
            existing_pressure_constraints.append(
                m.fs.PrimaryPumps[pump].max_lsrro_pressure_con
            )
            existing_pressure_constraints.append(
                m.fs.BoosterPumps[pump].max_ro_pressure_con
            )
        # add these to our sense manager
        m.fs.sense_manager.sensitivities["pump_pressure"].existing_constraints = (
            existing_pressure_constraints
        )
        m.fs.sense_manager.activate_sensitivities()

    if quick_start:
        m.fs.water_recovery.fix(0.45)
        m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].fix(35)
        m.fs.feed.properties[0].flow_vol_phase["Liq"].fix()
        lsrro.solve(m, tee=True)

        print("Quick start initialization complete.")
    if feed_tds != 35 or recovery != 0.45:
        m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].fix(feed_tds)
        m.fs.water_recovery.fix(recovery)
        lsrro.solve(m, tee=True)
        print(
            "solved to target tds and recovery of {} g/L and {} %".format(
                feed_tds, recovery
            )
        )
    if m.target_tds != 35 or m.target_recovery != 0.45:
        m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].fix(m.target_tds)
        m.fs.water_recovery.fix(m.target_recovery)
        lsrro.solve(m, tee=True)
        print(
            "solved to target tds and recovery of {} g/L and {} %".format(
                m.target_tds, m.target_recovery
            )
        )


def add_sensitivity_analysis(m):
    """Setup sensitivity manager process block and
    register sensitivity variables for LSRRO model"""
    m.fs.sense_manager = SensitivityManager()

    # aggregate all pressures and efficiencies to register for sensitivity analysis
    pump_pressure = []
    pump_eff = []
    pump_eff.append(m.fs.PrimaryPumps[1].efficiency_pump[0.0])

    for pump in m.fs.LSRRO_Stages:
        pump_pressure.append(m.fs.PrimaryPumps[pump].outlet.pressure[0])
        pump_eff.append(m.fs.PrimaryPumps[pump].efficiency_pump[0.0])

    for pump in m.fs.LSRRO_Stages:
        pump_pressure.append(m.fs.BoosterPumps[pump].outlet.pressure[0])
        pump_eff.append(m.fs.BoosterPumps[pump].efficiency_pump[0.0])

    m.fs.sense_manager.register_sensitivity(
        sensitivity_name="pump_pressure",
        model_variables=pump_pressure,
        sensitivity_type="upper_bound",
        remove_bounds=True,
    )

    # aggreagete all LSRRO A values
    a_values = []
    for ro in m.fs.LSRRO_Stages:
        a_values.append(m.fs.ROUnits[ro].A_comp[0, "H2O"])

    # register each one
    m.fs.sense_manager.register_sensitivity(
        sensitivity_name="a_value",
        model_variables=a_values,
        sensitivity_type="equality",
        remove_bounds=True,
    )
    m.fs.sense_manager.register_sensitivity(
        sensitivity_name="lsrro_membrane_cost",
        model_variables=m.fs.costing.reverse_osmosis.high_pressure_membrane_cost,
        sensitivity_type="equality",
        remove_bounds=True,
    )
    m.fs.sense_manager.register_sensitivity(
        sensitivity_name="pump_cost",
        model_variables=m.fs.costing.high_pressure_pump.cost,
        sensitivity_type="equality",
        remove_bounds=True,
    )
    m.fs.sense_manager.register_sensitivity(
        sensitivity_name="pump_efficiency",
        model_variables=pump_eff,
        sensitivity_type="equality",
        remove_bounds=True,
    )

    # fix and scale them
    m.fs.sense_manager.fix_and_scale()

    # generate templates for our use
    m.fs.sense_manager.generate_multiplier_yaml_template(
        "lsrro_sensitivity_template.yaml",
        multiplier_lb=0.8,
        multiplier_ub=1.2,
    )
    m.fs.sense_manager.generate_absolute_svoi_template(
        "svoi_sweep_template_absolute.yaml",
    )
    # display current sensitivity values
    m.fs.sense_manager.display_sensitivities()


def solve_model(m, **kwargs):
    """Function to optimize LSRRO model, this function should solve the model after it has been initialized
    Args:
        m (ConcreteModel): LSRRO model to be optimized
        **kwargs: additional keyword arguments (These are required for use with paramter sweep and loop tool!)
    Returns:
        result (SolverResults): result of the solver
    """
    result = lsrro.solve(m, tee=False)
    return result


def feasibility_test_function(m, **kwargs):
    """Function to test feasibility of LSRRO model, this function should return True if the model is feasible
    Args:
        m (ConcreteModel): LSRRO model to be tested
        permeate_quality_limit (float): permeate quality limit in kg solute/kg solution
        **kwargs: additional keyword arguments (These are required for use with paramter sweep and loop tool!)
    Returns:
        feasible (bool): True if the model is feasible, False otherwise
    """
    # Example feasibility test: check if permeate concentration is below the limit

    solver = get_solver()
    solver.solve(m.fs.feed)
    salt = m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].value
    recovery = m.fs.water_recovery.value
    water = m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].value
    salt = salt
    mass_fraction = salt / (water * (1 - recovery) + salt)
    print(
        "Testing recovery of ",
        recovery,
        " with feed NaCl of ",
        m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value,
    )
    if mass_fraction > 0.2648:
        print(
            "Would you like some LSRRO with your salt?",
            "Recovery too high for feed NaCl {} and  recovery is {}".format(
                m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value,
                recovery,
            ),
        )
        return False
    else:
        print(
            "Treatment possible for feed NaCl {} and  recovery is {}".format(
                m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].value,
                recovery,
            )
        )
        return True


def test_main_functions():
    # for i in range(5):
    i = 5
    m = build_function(number_of_stages=i + 1, add_sensitivity=True)
    initialize_model(
        m,
        quick_start=True,
        A_value=1.38e-11,
        permeate_quality_limit=1000.0e-6,
    )
    m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].fix(35)
    m.fs.water_recovery.fix(0.45)
    feasible = feasibility_test_function(
        m,
    )
    if feasible:
        m = solve_model(m)
        print(f"Model with {i+1} stages solved successfully.")
    else:
        print(f"Model with {i+1} stages is not feasible.")


def stage_sweep_analysis():
    """Main function to run loop tool analysis on LSRRO model, this will use
    the lssrro_stage_sweep.yaml file to define the parameter sweep and run loop tool in parallel mode,
    where each stage will be run on its own thread, allowing to solve all stages in
    parallel and speed up the analysis.
    """
    cpu_count = multiprocessing.cpu_count() - 2
    if cpu_count > 10:
        cpu_count = 10
    print(f"Working in {get_working_dir()}")
    loopTool(
        get_working_dir() + "/lsrro_stage_sweep.yaml",
        build_function=build_function,
        initialize_function=initialize_model,
        optimize_function=solve_model,
        save_name="lsrro_stage_sweep",
        probe_function=feasibility_test_function,
        saving_dir=get_working_dir(),
        number_of_subprocesses=1,
        num_loop_workers=cpu_count,
    )


def sensitivity_sweep_analysis():
    """Main function to run loop tool analysis on LSRRO model, this will use
    the lsrro_sensitivity_sweep.yaml file to define the parameter sweep and run loop tool in parallel mode,
    where each stage will be run on its own thread, allowing to solve all stages in
    parallel and speed up the analysis.
    """
    cpu_count = multiprocessing.cpu_count() - 2
    if cpu_count > 10:
        cpu_count = 10
    print(f"Working in {get_working_dir()}")
    loopTool(
        get_working_dir() + "/lsrro_sensitivity_sweep.yaml",
        build_function=build_function,
        initialize_function=initialize_model,
        optimize_function=solve_model,
        save_name="lsrro_sensitivity_sweep",
        probe_function=feasibility_test_function,
        saving_dir=get_working_dir(),
        number_of_subprocesses=1,
        num_loop_workers=cpu_count,
    )


if __name__ == "__main__":
    if False:
        ### Use this to verify all our configured functions work properly
        test_main_functions()
    elif False:
        ### Use this to run the full analysis with loop tool
        stage_sweep_analysis()
    elif True:
        sensitivity_sweep_analysis()
