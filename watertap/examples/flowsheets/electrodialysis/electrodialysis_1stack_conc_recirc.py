#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

from pyomo.environ import (
    ConcreteModel,
    Var,
    value,
    Constraint,
    Expression,
    Objective,
    TransformationFactory,
    assert_optimal_termination,
    units as pyunits,
    NonNegativeReals,
)
from pyomo.network import Arc

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state
import idaes.core.util.model_statistics as mstat
from idaes.models.unit_models import Feed, Product, Separator, Mixer
from watertap.unit_models.pressure_changer import Pump
from idaes.models.unit_models.mixer import MixingType, MomentumMixingType
import pandas as pd
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslogger
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

__author__ = "Xiangyu Bi"
_log = idaeslogger.getLogger(__name__)

"""
This flowsheet implements a single-stack electrodialysis desaliantion system in a 
feed-and-bleed mode, which recirculate a part of concentrate solution to achieve a
desired water reovery > 50%. 
"""


def main():
    m = build()
    solver = get_solver()

    ## Simulate a fully defined operation
    # A fully-defined system simulation, dof=0
    init_arg = {
        ("flow_vol_phase", ("Liq")): 5.2e-4,
        ("conc_mol_phase_comp", ("Liq", "Na_+")): 34.188,
        ("conc_mol_phase_comp", ("Liq", "Cl_-")): 34.188,
    }  # Corresponding to C_feed = 2g/L
    m.fs.feed.properties.calculate_state(
        init_arg,
        hold_state=True,
    )
    m.fs.EDstack.voltage_applied[0].fix(10)
    m.fs.recovery_vol_H2O.fix(0.7)
    _condition_base(m)

    # Initialize and solve the model
    initialize_system(m, solver=solver)
    solve(m, solver=solver)
    print("\n***---Fully-defined simulation results, Fixed inlet and voltage---***")
    display_model_metrics(m)

    # Perform an optimization over cell_length, cell_pair_mum, and voltage_applied

    ed = m.fs.EDstack
    ulim = (
        ed.voltage_x[0, 0].value
        / ed.current_density_x[0, 0].value
        * ed.current_dens_lim_x[0, 0].value
        * 1.5
    )
    m.fs.EDstack.voltage_applied[0].unfix()
    m.fs.EDstack.voltage_applied[0].setlb(0.1)
    m.fs.EDstack.voltage_applied[0].setub(ulim)
    m.fs.EDstack.cell_pair_num.unfix()
    m.fs.EDstack.cell_pair_num.setlb(10)
    m.fs.EDstack.cell_pair_num.setub(1000)
    m.fs.EDstack.cell_length.unfix()
    m.fs.prod.properties[0].conc_mol_phase_comp["Liq", "Na_+"].fix(
        1.7094
    )  # Corresponding to C_product = 100 ppm
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)
    solve(m, solver=solver, tee=True)
    m.fs.EDstack.cell_pair_num.fix(round(m.fs.EDstack.cell_pair_num.value))
    solve(m, solver=solver, tee=True)
    print("\n***---Optimization results, Product conc = 100 ppb---***")
    display_model_metrics(m)

    return m


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
    m.fs.pump0.ratioP_calculation.deactivate()
    m.fs.pump1 = Pump(property_package=m.fs.properties)
    m.fs.pump0.ratioP_calculation.deactivate()

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

    # touch necessary variables to ensure they are constructed and solved
    m.fs.prod.properties[0].flow_vol_phase[...]
    m.fs.disp.properties[0].flow_vol_phase[...]

    # costing
    m.fs.EDstack.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method_arguments={"cost_electricity_flow": True, "has_rectifier": True},
    )
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

    m.fs.feed_salinity = Var(
        initialize=1,
        domain=NonNegativeReals,
        units=pyunits.kg * pyunits.m**-3,
        doc="Salinity of feed solution. Tested only with NaCl.",
    )

    m.fs.eq_recovery_vol_H2O = Constraint(
        expr=m.fs.recovery_vol_H2O * m.fs.feed.properties[0].flow_vol_phase["Liq"]
        == m.fs.prod.properties[0].flow_vol_phase["Liq"]
    )

    m.fs.eq_electrodialysis_equal_flow = Constraint(
        expr=m.fs.EDstack.diluate.properties[0, 0].flow_vol_phase["Liq"]
        == m.fs.EDstack.concentrate.properties[0, 0].flow_vol_phase["Liq"]
    )
    m.fs.eq_feed_salinity = Constraint(
        expr=m.fs.feed_salinity
        == sum(
            m.fs.feed.properties[0].conc_mass_phase_comp["Liq", j]
            for j in m.fs.properties.ion_set
        )
    )
    m.fs.product_salinity = Expression(
        expr=sum(
            m.fs.prod.properties[0].conc_mass_phase_comp["Liq", j]
            for j in m.fs.properties.solute_set
        )
    )
    m.fs.disposal_salinity = Expression(
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

    m.fs.voltage_per_cp = Expression(
        expr=m.fs.EDstack.voltage_applied[0] / m.fs.EDstack.cell_pair_num
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
    # ---specifications---
    # Here is simulated a scenario of a defined EDstack and
    # specific water recovery and product salinity.
    m.fs.feed.properties[0].pressure.fix(101325)
    m.fs.feed.properties[0].temperature.fix(298.15)
    m.fs.pump0.control_volume.properties_in[0].pressure.fix(101325)
    m.fs.pump1.efficiency_pump.fix(0.8)
    m.fs.pump0.efficiency_pump.fix(0.8)

    m.fs.prod.properties[0].pressure.fix(101325)
    m.fs.disp.properties[0].pressure.fix(101325)
    m.fs.disp.properties[0].temperature.fix(298.15)

    # Set ED unit vars
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
    m.fs.EDstack.cell_pair_num.fix(56)
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


def solve(m, solver=None, tee=True, check_termination=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(m, tee=tee)
    if check_termination:
        assert_optimal_termination(results)
    return results


def initialize_system(m, solver=None):
    # set up solver
    if solver is None:
        solver = get_solver()
    optarg = solver.options
    # set scaling factors for state vars and call the 'calculate_scaling_factors' function
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 0.1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e2, index=("Liq", "Na_+")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e2, index=("Liq", "Cl_-")
    )
    iscale.set_scaling_factor(m.fs.EDstack.cell_width, 5)
    iscale.set_scaling_factor(m.fs.EDstack.cell_length, 1)
    iscale.set_scaling_factor(m.fs.EDstack.cell_pair_num, 0.1)
    iscale.set_scaling_factor(m.fs.EDstack.voltage_applied, 1)
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
    iscale.constraint_scaling_transform(
        m.fs.eq_feed_salinity,
        iscale.get_scaling_factor(
            m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "Na_+"]
        ),
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
    propagate_state(m.fs.arc4)
    m.fs.prod.initialize(optarg=optarg)
    propagate_state(m.fs.arc5)
    m.fs.prod.initialize(optarg=optarg)
    m.fs.sepa1.split_fraction[0, "to_conc_in1"].fix(
        max(2 - value(m.fs.recovery_vol_H2O) ** -1, 1e-8)
    )
    m.fs.sepa1.initialize(optarg=optarg)
    m.fs.sepa1.split_fraction[0, "to_conc_in1"].unfix()
    propagate_state(m.fs.arc6)
    m.fs.disp.initialize()
    propagate_state(m.fs.arc7)
    propagate_state(m.fs.arc3b, direction="backward")
    propagate_state(m.fs.arc2)
    m.fs.mix0.initialize(optarg=optarg)
    m.fs.costing.initialize()

    iscale.calculate_scaling_factors(m)


def display_model_metrics(m):

    print("---Report the ED unit---")

    m.fs.EDstack.report()

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
        "Product": [
            value(m.fs.prod.properties[0].flow_vol_phase["Liq"]),
            value(m.fs.product_salinity),
        ],
        "Disposal": [
            value(m.fs.disp.properties[0].flow_vol_phase["Liq"]),
            value(m.fs.disposal_salinity),
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
            value(m.fs.EDstack.cell_length),
            value(m.fs.EDstack.voltage_applied[0]),
            value(m.fs.costing.specific_energy_consumption),
            value(m.fs.costing.LCOW),
        ],
        columns=["value"],
        index=[
            "Water recovery by volume",
            "Total membrane area (aem or cem), m2",
            "ED cell pair number",
            "ED cell flow path length",
            "Operation Voltage, V",
            "Specific energy consumption, kWh/m3",
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
    m = main()
