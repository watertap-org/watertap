###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

from pyomo.environ import (
    ConcreteModel,
    Var,
    value,
    Constraint,
    Objective,
    Expression,
    TransformationFactory,
    assert_optimal_termination,
    units as pyunits,
    NonNegativeReals,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import (
    propagate_state,
)
import idaes.core.util.model_statistics as mstat
from pyomo.core.expr.current import identify_variables
from idaes.models.unit_models import Feed, Product, Separator, Mixer
from watertap.unit_models.pressure_changer import Pump
from idaes.models.unit_models.mixer import MixingType, MomentumMixingType
import pandas as pd
import numpy as np
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslogger
from pytest import approx
from watertap.core.util.initialization import check_dof
from watertap.unit_models.electrodialysis_1D import (
    ElectricalOperationMode,
    PressureDropMethod,
    FrictionFactorMethod,
    HydraulicDiameterMethod,
    LimitingCurrentDensityMethod,
)
from watertap.unit_models.electrodialysis_1D import Electrodialysis1D
from watertap.costing.watertap_costing_package import WaterTAPCosting
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
import idaes.core.util.model_diagnostics as m_diag

__author__ = "Xiangyu Bi"

_log = idaeslogger.getIdaesLogger(__name__)


def main():
    m = build()
    m.obj = Objective(expr=0)
    deci_var_dict = {
        "voltage_applied[0]": 10,
        "recovery_vol_H2O": 0.65,
        "channel_height": 1e-3,
    }
    initarg = make_initarg_list([9],flow_rate_vol=1.04e-3)
    initialize_dof0_system(
        m=m, initargs=initarg[0], solve_after_init=True, **deci_var_dict
    )
    print("===INITIALIZE OUTCOME")
    display_model_metrics(m)
    opt_var_dict = {
        "voltage_applied[0]": (10, 0, value(m.fs.voltage_lim)),
        "cell_pair_num": (100, 1, 10000),
        "cell_length": (1.68, 0.01, 10),
        "cell_width": (0.197, 0.01, 10),
        #"channel_height": (7e-4, 1e-4, 1e-2),
    }
    if not m.find_component("obj") is None:
        m.del_component("obj")
        print("obj==0 deleted")
    opt_res = optimize_LCOW_fixed_prod_salinity(m, prod_sal=0.5, **opt_var_dict)
    print("===OPTIMIZE OUTCOME")
    display_model_metrics(m)
    print(f"Final opt solver condition is {opt_res.solver.termination_condition}.")
    dh = m_diag.DegeneracyHunter(m)
    dh.check_residuals(tol=1e-8)


def make_initarg_list(conc_mass_list, mw=0.0585, flow_rate_vol=5.2e-4):
    conc_mass_in = pd.DataFrame(data=conc_mass_list, columns=["C0"])  # g/L
    conc_mol_in = conc_mass_in / mw  # mol m-3
    initarg = []
    for k in conc_mol_in["C0"]:
        initarg.append(
            {
                ("flow_vol_phase", ("Liq")): flow_rate_vol,
                ("conc_mol_phase_comp", ("Liq", "Na_+")): k,
                ("conc_mol_phase_comp", ("Liq", "Cl_-")): k,
            }
        )
    return initarg


def initialize_dof0_system(
    m=None,
    solve_after_init=False,
    terminate_nonoptimal_sol=False,
    report_bad_scaling=False,
    initargs={},
    **deciargs,
):
    """
    To initialize a system at zero dof and give the option to set up the
    decision variables to be controled.

    Keyword Arguments:
        m : model to be initalized.
        solve_after_init: Bool arg to control whether the model is to be solved after setting
                            all intial values.
        terminate_nonoptimal_sol: Bool arg to control whether exception is triggered upon a
                                    non-optimal solution
        initargs : list of dicts of initial states of the feed solution; the dict
                    should have the same length to the state vars.
        report_bad_scaling: control over whether checking badly scaled vars is performed.
        **deciargs: decision vars to be fixed, as a {name: val} dict.

    Returns: a list of the fixed decision vars and the solver results if solve_after_init
                == True.
    """
    if m == None:
        m = build()
    m.fs.feed.properties.calculate_state(
        initargs,
        hold_state=True,
    )
    _condition_base(m)
    deci_comp_list = []
    for deci_var in deciargs:
        comp = m.fs.find_component(str(deci_var))
        try:
            assert not comp is None
        except:
            comp = m.fs.EDstack.find_component(str(deci_var))
            try:
                assert not comp is None
            except:
                raise TypeError(
                    "Var {} in the provided decision_var dict are not found in the model's component".format(
                        deci_var
                    )
                )
        comp.fix(deciargs[deci_var])
        deci_comp_list.append(comp)

    assert mstat.degrees_of_freedom(m) == 0
    try:
        initialize_system(m, report_bad_scaling=report_bad_scaling)
    except Exception as experr:
        _log.warning("Initialization Fails in initialize_dof0_system:{}".format(experr))
        pass
    if solve_after_init:
        res = solve(m, check_termination=terminate_nonoptimal_sol)
        return [deci_comp_list, res]
    else:
        return deci_comp_list


def optimize_LCOW_fixed_prod_salinity(m, tee=False, prod_sal=0.1, mw=0.0585, **optargs):
    """
    Function to optimize LCOW of a defined ED system with a definied product water salinity
    target.  The system should first be initialized to have zero dof before running this method.
    Keyword Arguments:
        m : model to be optimized.
        prod_sal: salinity of product water fixed for optimization
        **optargs: dict argument containing vars to be optimized; The dict member should
                    be in form of {var: (val, lb, ub)}, where (val, lb, ub) passes the initial
                    value, lower bound and upper bound of the var to be optimized. The user is
                    allowed to skip lb and ub, in which case the var will be optimized with no
                    bounds. However, providing lb and ub is encouraged.

    Returns: solver results of the optimization.
    """
    if not mstat.degrees_of_freedom(m) == 0:
        raise TypeError(
            "A model is expected to be fully defined at zero dof before being optimized."
        )
    conc_mol = prod_sal / mw
    m.fs.prod.properties[0].conc_mol_phase_comp["Liq", "Na_+"].fix(conc_mol)
    for opt_var in optargs:
        var = m.fs.find_component(str(opt_var))
        try:
            assert not var is None
        except:
            var = m.fs.EDstack.find_component(str(opt_var))
            try:
                assert not var is None
            except:
                raise TypeError(
                    "Var {} in the provided opt_var dict are not found in the model's component".format(
                        opt_var
                    )
                )
        var.unfix()
        if len(optargs[opt_var]) == 1:
            var.set_value(optargs[opt_var][0])
        elif len(optargs[opt_var]) == 3:
            var.set_value(optargs[opt_var][0])
            var.setlb(optargs[opt_var][1])
            var.setub(optargs[opt_var][2])
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)
    assert mstat.degrees_of_freedom(m) == len(optargs) - 1
    result = solve(m, tee=tee)
    if not m.fs.EDstack.cell_pair_num.value.is_integer():
        _log.warning(
            "The ED cell pair number is not a integer after optimization "
            "and the model is to be re-solved at its closest integer value."
        )
        print("===Outcome of first optimization===")
        display_model_metrics(m)
        m.fs.EDstack.cell_pair_num.fix(round(m.fs.EDstack.cell_pair_num.value))
    result = solve(m, tee=tee)
    if result.solver.termination_condition == "maxIterations":
        while result.solver.termination_condition == "maxIterations":
            result = solve(m, tee=tee)
    return result


def build():
    # ---building model---
    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)
    ion_dict = {
        "solute_list": ["Na_+", "Cl_-"],
        "mw_data": {"H2O": 18e-3, "Na_+": 23e-3, "Cl_-": 35.5e-3},
        "elec_mobility_data": {("Liq", "Na_+"): 5.19e-8, ("Liq", "Cl_-"): 7.92e-8},
        "diffusivity_data": {("Liq", "Na_+"): 1.33e-9, ("Liq", "Cl_-"): 2.03e-9},
        "charge": {"Na_+": 1, "Cl_-": -1},
    }
    m.fs.properties = MCASParameterBlock(**ion_dict)
    m.fs.costing = WaterTAPCosting()
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.sepa0 = Separator(
        property_package=m.fs.properties,
        outlet_list=["to_dil_in", "to_conc_in0"],
    )
    m.fs.mix0 = Mixer(
        property_package=m.fs.properties,
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.none,
        inlet_list=["from_feed", "from_conc_out"],
    )

    m.fs.pump0 = Pump(property_package=m.fs.properties)
    m.fs.pump0.del_component("ratioP")
    m.fs.pump0.del_component("ratioP_calculation")
    m.fs.pump1 = Pump(property_package=m.fs.properties)
    m.fs.pump1.del_component("ratioP")
    m.fs.pump1.del_component("ratioP_calculation")

    # Add electrodialysis (ED) stacks
    m.fs.EDstack = Electrodialysis1D(
        property_package=m.fs.properties,
        operation_mode=ElectricalOperationMode.Constant_Voltage,
        finite_elements=20,
        has_pressure_change=True,
        has_nonohmic_potential_membrane=True,
        has_Nernst_diffusion_layer=True,
        limiting_current_density_method=LimitingCurrentDensityMethod.Theoretical,
        pressure_drop_method=PressureDropMethod.Darcy_Weisbach,
        hydraulic_diameter_method=HydraulicDiameterMethod.spacer_specific_area_known,
        friction_factor_method=FrictionFactorMethod.Gurreri,
    )
    m.fs.sepa1 = Separator(
        property_package=m.fs.properties,
        outlet_list=["to_disp", "to_conc_in1"],
    )
    m.fs.prod = Product(property_package=m.fs.properties)
    m.fs.disp = Product(property_package=m.fs.properties)

    # Touching needed variables for initialization and displaying results
    m.fs.prod.properties[0].flow_vol_phase[...]
    m.fs.disp.properties[0].flow_vol_phase[...]
    m.fs.EDstack.diluate.properties[...].flow_vol_phase[...]
    m.fs.EDstack.concentrate.properties[...].flow_vol_phase[...]

    # costing
    m.fs.EDstack.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.pump0.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.pump1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(
        m.fs.prod.properties[0].flow_vol_phase["Liq"]
    )
    m.fs.costing.add_LCOW(m.fs.prod.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(
        m.fs.prod.properties[0].flow_vol_phase["Liq"]
    )

    # add extra variables and constraints
    m.fs.recovery_vol_H2O = Var(
        initialize=0.5,
        bounds=(0, 1),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="flowsheet level water recovery calculated by volumeric flow rate",
    )

    m.fs.eq_recovery_vol_H2O = Constraint(
        expr=m.fs.recovery_vol_H2O * m.fs.feed.properties[0].flow_vol_phase["Liq"]
        == m.fs.prod.properties[0].flow_vol_phase["Liq"]
    )

    m.fs.eq_electrodialysis_equal_flow = Constraint(
        expr=m.fs.EDstack.diluate.properties[0, 0].flow_vol_phase["Liq"]
        == m.fs.EDstack.concentrate.properties[0, 0].flow_vol_phase["Liq"]
    )
    m.fs.prod_salinity = Expression(
        expr=sum(
            m.fs.prod.properties[0].conc_mass_phase_comp["Liq", j]
            for j in m.fs.properties.solute_set
        )
    )
    m.fs.disp_salinity = Expression(
        expr=sum(
            m.fs.disp.properties[0].conc_mass_phase_comp["Liq", j]
            for j in m.fs.properties.solute_set
        )
    )

    m.fs.mem_area = Expression(
        expr=m.fs.EDstack.cell_width
        * m.fs.EDstack.cell_length
        * m.fs.EDstack.cell_pair_num
    )

    m.fs.voltage_lim = Expression(
        expr=m.fs.EDstack.voltage_x[0, 0].value
        / m.fs.EDstack.current_density_x[0, 0].value
        * m.fs.EDstack.current_dens_lim_x[0, 0].value
    )
    m.fs.voltage_per_cp = Expression(
        expr=m.fs.EDstack.voltage_applied[0].value / m.fs.EDstack.cell_pair_num
    )

    # Add Arcs
    m.fs.arc0 = Arc(source=m.fs.feed.outlet, destination=m.fs.sepa0.inlet)
    m.fs.arc1b = Arc(source=m.fs.sepa0.to_dil_in, destination=m.fs.pump1.inlet)
    m.fs.arc1f = Arc(source=m.fs.pump1.outlet, destination=m.fs.EDstack.inlet_diluate)
    m.fs.arc2 = Arc(
        source=m.fs.sepa0.to_conc_in0,
        destination=m.fs.mix0.from_feed,
    )
    m.fs.arc3b = Arc(source=m.fs.mix0.outlet, destination=m.fs.pump0.inlet)
    m.fs.arc3f = Arc(
        source=m.fs.pump0.outlet, destination=m.fs.EDstack.inlet_concentrate
    )
    m.fs.arc4 = Arc(source=m.fs.EDstack.outlet_diluate, destination=m.fs.prod.inlet)
    m.fs.arc5 = Arc(
        source=m.fs.EDstack.outlet_concentrate, destination=m.fs.sepa1.inlet
    )
    m.fs.arc6 = Arc(source=m.fs.sepa1.to_disp, destination=m.fs.disp.inlet)
    m.fs.arc7 = Arc(source=m.fs.sepa1.to_conc_in1, destination=m.fs.mix0.from_conc_out)
    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def _condition_base(m):
    """
    Internal function to set up a base condition of the ED system. This gives a zero dof
    condition at reasonable values of all variables. Users can set up different intial
    conditions by the "initialize_dof0_system" function.
    """
    # ---specifications---
    # Here is simulated a scenario of a defined EDstack and
    # specific water recovery and product salinity.
    m.fs.feed.properties[0].pressure.fix(101325)
    m.fs.feed.properties[0].temperature.fix(298.15)
    m.fs.pump1.control_volume.properties_in[0].pressure.fix(101325)
    m.fs.pump1.efficiency_pump.fix(0.8)
    m.fs.pump0.efficiency_pump.fix(0.8)

    m.fs.prod.properties[0].pressure.fix(101325)  #
    m.fs.disp.properties[0].pressure.fix(101325)
    m.fs.disp.properties[0].temperature.fix(298.15)

    # Set ED unit vars
    # Operational Properties
    m.fs.EDstack.voltage_applied[0].fix(10)
    m.fs.recovery_vol_H2O.fix(0.7)
    # membrane properties
    m.fs.EDstack.water_trans_number_membrane["cem"].fix(5.8)
    m.fs.EDstack.water_trans_number_membrane["aem"].fix(4.3)
    m.fs.EDstack.water_permeability_membrane["cem"].fix(2.16e-14)
    m.fs.EDstack.water_permeability_membrane["aem"].fix(1.75e-14)
    m.fs.EDstack.membrane_areal_resistance["cem"].fix(1.89e-4)
    m.fs.EDstack.membrane_areal_resistance["aem"].fix(1.77e-4)
    m.fs.EDstack.solute_diffusivity_membrane["cem", "Na_+"].fix(3.28e-11)
    m.fs.EDstack.solute_diffusivity_membrane["aem", "Na_+"].fix(3.28e-11)
    m.fs.EDstack.solute_diffusivity_membrane["cem", "Cl_-"].fix(3.28e-11)
    m.fs.EDstack.solute_diffusivity_membrane["aem", "Cl_-"].fix(3.28e-11)
    m.fs.EDstack.ion_trans_number_membrane["cem", "Na_+"].fix(1)
    m.fs.EDstack.ion_trans_number_membrane["aem", "Na_+"].fix(0)
    m.fs.EDstack.ion_trans_number_membrane["cem", "Cl_-"].fix(0)
    m.fs.EDstack.ion_trans_number_membrane["aem", "Cl_-"].fix(1)
    m.fs.EDstack.membrane_thickness["aem"].fix(1.3e-4)
    m.fs.EDstack.membrane_thickness["cem"].fix(1.3e-4)

    # Stack properties
    m.fs.EDstack.cell_pair_num.fix(100)
    m.fs.EDstack.channel_height.fix(7.1e-4)
    m.fs.EDstack.cell_width.fix(0.197)
    m.fs.EDstack.cell_length.fix(1.68)

    # Spacer properties
    m.fs.EDstack.spacer_porosity.fix(0.83)
    m.fs.EDstack.spacer_specific_area.fix(10400)

    # Electrochemical properties
    m.fs.EDstack.electrodes_resistance.fix(0)
    m.fs.EDstack.current_utilization.fix(1)
    m.fs.EDstack.diffus_mass.fix(1.6e-9)

    assert mstat.degrees_of_freedom(m) == 0


def solve(blk, solver=None, tee=True, check_termination=False):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    if check_termination:
        assert_optimal_termination(results)
    return results


def initialize_system(m, solver=None, report_bad_scaling=False):

    # set up solver
    if solver is None:
        solver = get_solver()
    optarg = solver.options
    # set scaling factors for state vars and call the 'calculate_scaling_factors' function
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 0.001, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e2, index=("Liq", "Na_+")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e2, index=("Liq", "Cl_-")
    )
    iscale.set_scaling_factor(m.fs.EDstack.cell_width, 5)
    iscale.set_scaling_factor(m.fs.EDstack.cell_length, 1)
    iscale.set_scaling_factor(m.fs.EDstack.cell_pair_num, 0.001)
    iscale.set_scaling_factor(m.fs.EDstack.voltage_applied, 0.1)
    iscale.set_scaling_factor(m.fs.pump0.control_volume.work, 1e1)
    iscale.set_scaling_factor(m.fs.pump1.control_volume.work, 1e1)

    iscale.calculate_scaling_factors(m)
    iscale.constraint_scaling_transform(
        m.fs.eq_recovery_vol_H2O,
        10 * iscale.get_scaling_factor(m.fs.feed.properties[0].flow_vol_phase["Liq"]),
    )
    iscale.constraint_scaling_transform(
        m.fs.eq_electrodialysis_equal_flow,
        10 * iscale.get_scaling_factor(m.fs.feed.properties[0].flow_vol_phase["Liq"]),
    )
    # populate intitial properties throughout the system
    m.fs.feed.initialize(optarg=optarg)
    propagate_state(m.fs.arc0)
    m.fs.sepa0.initialize(optarg=optarg)
    propagate_state(m.fs.arc1b)
    m.fs.pump1.deltaP[0].fix(2e5)
    m.fs.pump1.initialize()
    m.fs.pump1.deltaP[0].unfix()
    propagate_state(destination=m.fs.pump0.inlet, source=m.fs.pump1.inlet)
    m.fs.pump0.deltaP[0].fix(2e5)
    m.fs.pump0.initialize()
    m.fs.pump0.deltaP[0].unfix()
    propagate_state(m.fs.arc1f)
    propagate_state(m.fs.arc3f)
    m.fs.EDstack.initialize(optarg=optarg)
    m.fs.EDstack.report()
    propagate_state(m.fs.arc4)
    m.fs.prod.initialize(optarg=optarg)
    propagate_state(m.fs.arc5)
    m.fs.prod.initialize(optarg=optarg)
    init_sepa1_frac = max(2 - value(m.fs.recovery_vol_H2O) ** -1, 1e-8)
    m.fs.sepa1.split_fraction[0, "to_conc_in1"].fix(init_sepa1_frac)
    m.fs.sepa1.initialize(optarg=optarg)
    m.fs.sepa1.split_fraction[0, "to_conc_in1"].unfix()
    propagate_state(m.fs.arc6)
    m.fs.disp.initialize()
    propagate_state(m.fs.arc7)
    propagate_state(m.fs.arc3b, direction="backward")
    propagate_state(m.fs.arc2)
    m.fs.mix0.initialize(optarg=optarg)
    m.fs.costing.initialize()

    # Update essential scaling factors
    iscale.set_scaling_factor(
        m.fs.sepa1.split_fraction[0, "to_conc_in1"], init_sepa1_frac**-1
    )

    iscale.calculate_scaling_factors(m)
    if report_bad_scaling:
        print("BADLY SCALED VARS & CONSTRAINS")
        badly_scaled_var_values = {
            var.name: val for (var, val) in iscale.badly_scaled_var_generator(m)
        }
        for j, k in badly_scaled_var_values.items():
            print(j, ":", k)


def display_model_metrics(m):

    print("---Flow properties in feed, product and disposal---")

    fp = {
        "Feed": [
            value(m.fs.feed.properties[0].flow_vol_phase["Liq"]),
            value(
                sum(
                    m.fs.feed.properties[0].conc_mass_phase_comp["Liq", j]
                    for j in m.fs.properties.ion_set
                )
            ),
        ],
        "ED_diluate_in": [
            value(m.fs.EDstack.diluate.properties[0, 0].flow_vol_phase["Liq"]),
            "",
        ],
        "ED_concentrate_in": [
            value(m.fs.EDstack.concentrate.properties[0, 0].flow_vol_phase["Liq"]),
            "",
        ],
        "Product": [
            value(m.fs.prod.properties[0].flow_vol_phase["Liq"]),
            value(m.fs.prod_salinity),
        ],
        "Disposal": [
            value(m.fs.disp.properties[0].flow_vol_phase["Liq"]),
            value(m.fs.disp_salinity),
        ],
    }
    fp_table = pd.DataFrame(
        data=fp,
        index=["Volume Flow Rate (m3/s)", "Total Ion Mass Concentration (kg/m3)"],
    )
    print(fp_table)

    print("---Performance Metrics---")

    pm_table = pd.DataFrame(
        data=[
            value(m.fs.recovery_vol_H2O),
            value(m.fs.mem_area),
            value(m.fs.EDstack.cell_pair_num),
            value(m.fs.EDstack.channel_height),
            value(m.fs.EDstack.cell_length),
            value(m.fs.EDstack.cell_width),
            value(m.fs.EDstack.voltage_applied[0]),
            value(m.fs.voltage_per_cp),
            value(m.fs.costing.specific_energy_consumption),
            value(m.fs.EDstack.specific_power_electrical[0]),
            value(m.fs.costing.LCOW),
        ],
        columns=["value"],
        index=[
            "Water recovery by volume",
            "Total membrane area (aem or cem), m2",
            "ED cell pair number",
            "ED channel height, m",
            "ED cell flow path length",
            "ED cell width",
            "Operation voltage, V",
            "Cell-pair voltage, V",
            "Specific energy consumption, kWh/m3",
            "Specific energy consumption on ED stack, kWh/m3",
            "Levelized cost of water, $/m3",
        ],
    )
    print(pm_table)
    print("---Pressure and Temperature point checking---")

    pt_dict = {
        "Feed": (
            value(m.fs.feed.outlet.pressure[0]),
            value(m.fs.feed.outlet.temperature[0]),
        ),
        "ED_in_dil": (
            value(m.fs.EDstack.inlet_diluate.pressure[0]),
            value(m.fs.EDstack.inlet_diluate.temperature[0]),
        ),
        "ED_in_conc": (
            value(m.fs.EDstack.inlet_concentrate.pressure[0]),
            value(m.fs.EDstack.inlet_concentrate.temperature[0]),
        ),
        "ED_out_dil": (
            value(m.fs.EDstack.outlet_diluate.pressure[0]),
            value(m.fs.EDstack.outlet_diluate.temperature[0]),
        ),
        "ED_out_conc": (
            value(m.fs.EDstack.outlet_concentrate.pressure[0]),
            value(m.fs.EDstack.outlet_concentrate.temperature[0]),
        ),
        "Prod": (
            value(m.fs.prod.inlet.pressure[0]),
            value(m.fs.prod.inlet.temperature[0]),
        ),
        "Disp": (
            value(m.fs.disp.inlet.pressure[0]),
            value(m.fs.disp.inlet.temperature[0]),
        ),
    }
    pt_table = pd.DataFrame(data=pt_dict, index=["Pressure (Pa)", "Temperature (K)"])
    print(pt_table)


if __name__ == "__main__":
    main()
