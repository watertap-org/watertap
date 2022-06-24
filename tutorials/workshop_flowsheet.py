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
    value,
    Constraint,
    Expression,
    Objective,
    Param,
    TransformationFactory,
    units,
    assert_optimal_termination,
)
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import solve_indexed_blocks, propagate_state
from idaes.models.unit_models import Mixer, Separator, Product, Feed
from idaes.core import UnitModelCostingBlock
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

import watertap.property_models.seawater_prop_pack as props
from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice
from watertap.core.util.initialization import assert_degrees_of_freedom
from watertap.costing import WaterTAPCosting

# for UI:
from watertap.ui.api import export_variables, FlowsheetInterface, WorkflowActions


def flowsheet_interface() -> FlowsheetInterface:
    """Define the interface to the flowsheet for the UI layer.

    Example usage::

        from tutorials import workshop_flowsheet
        from watertap.ui import api
        import json

        fsi = workshop_flowsheet.flowsheet_interface()
        print(f"Using flowsheet: {fsi.name}")

        build, solve = api.WorkflowActions.build, api.WorkflowActions.solve

        fsi.run_action(build) # build the flowsheet

        print(json.dumps(fsi.dict(), indent=2))  # show exported variables

        results = fsi.run_action(solve)  # solve the flowsheet
        print(results) # print results obj
    """
    fsi = FlowsheetInterface(
        {
            "display_name": "Example RO Flowsheet",
            "description": "Example RO flowsheet for " "workshop tutorial",
        }
    )
    fsi.set_action(WorkflowActions.build, ui_build)
    fsi.set_action(WorkflowActions.solve, ui_solve)
    return fsi


def ui_build(ui=None, **kwargs):
    model = build()
    set_operating_conditions(model)
    initialize_system(model)

    # initial solution of square problem
    solve(model)

    ui.set_block(model.fs)
    export_ui_variables(model.fs)


def ui_solve(block=None, **kwargs):
    fs, m = block, block.parent_block()

    # optimize
    optimize_set_up(m)
    optimize(m)

    return display_ui_output(m)


def main():
    # build, set, and initialize
    m = build()
    set_operating_conditions(m)
    initialize_system(m)

    # simulate and display
    solve(m)
    print("\n***---Simulation results---***")
    display_system(m)
    display_design(m)
    display_state(m)
    print(display_ui_output(m))

    # optimize and display
    optimize_set_up(m)
    assert_units_consistent(m)
    optimize(m)
    print("\n***---Optimization results---***")
    display_system(m)
    display_design(m)
    display_state(m)

    # change one parameter and see effect
    m.fs.costing.reverse_osmosis_membrane_cost.fix(60)
    optimize(m)
    print("\n***---Parameter change results---***")
    display_system(m)
    display_design(m)
    display_state(m)
    print(display_ui_input(m))
    print(display_ui_output(m))


def build():
    # flowsheet set up
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = props.SeawaterParameterBlock()
    m.fs.costing = WaterTAPCosting()

    # unit models
    m.fs.feed = Feed(default={"property_package": m.fs.properties})
    m.fs.pump = Pump(default={"property_package": m.fs.properties})
    m.fs.RO = ReverseOsmosis0D(
        default={
            "property_package": m.fs.properties,
            "has_pressure_change": True,
            "pressure_change_type": PressureChangeType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.calculated,
            "concentration_polarization_type": ConcentrationPolarizationType.calculated,
        }
    )
    m.fs.erd = EnergyRecoveryDevice(
        default={
            "property_package": m.fs.properties,
        },
    )
    m.fs.product = Product(default={"property_package": m.fs.properties})
    m.fs.disposal = Product(default={"property_package": m.fs.properties})
    # costing
    m.fs.pump.costing = UnitModelCostingBlock(
        default={"flowsheet_costing_block": m.fs.costing}
    )
    m.fs.RO.costing = UnitModelCostingBlock(
        default={"flowsheet_costing_block": m.fs.costing}
    )
    m.fs.erd.costing = UnitModelCostingBlock(
        default={
            "flowsheet_costing_block": m.fs.costing,
            "costing_method_arguments": {
                "energy_recovery_device_type": "pressure_exchanger"
            },
        }
    )
    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(m.fs.product.properties[0].flow_vol)

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.pump.inlet)
    m.fs.s02 = Arc(source=m.fs.pump.outlet, destination=m.fs.RO.inlet)
    m.fs.s03 = Arc(source=m.fs.RO.permeate, destination=m.fs.product.inlet)
    m.fs.s04 = Arc(source=m.fs.RO.retentate, destination=m.fs.erd.inlet)
    m.fs.s05 = Arc(source=m.fs.erd.outlet, destination=m.fs.disposal.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    # set default property values
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "TDS")
    )
    # set unit model values
    iscale.set_scaling_factor(m.fs.pump.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.erd.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.RO.area, 1e-2)
    # touch properties used in specifying the model
    m.fs.feed.properties[0].flow_vol_phase["Liq"]
    m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"]
    m.fs.disposal.properties[0].flow_vol_phase["Liq"]
    m.fs.disposal.properties[0].mass_frac_phase_comp["Liq", "TDS"]
    # calculate and propagate scaling factors
    iscale.calculate_scaling_factors(m)

    return m


def set_operating_conditions(m):

    # ---specifications---
    # feed
    # state variables
    m.fs.feed.properties[0].pressure.fix(101325)  # feed pressure [Pa]
    m.fs.feed.properties[0].temperature.fix(273.15 + 25)  # feed temperature [K]
    # properties (cannot be fixed for initialization routines, must calculate the state variables)
    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 1e-3,  # feed volumetric flow rate [m3/s]
            ("mass_frac_phase_comp", ("Liq", "TDS")): 0.035,
        },  # feed TDS mass fraction [-]
        hold_state=True,  # fixes the calculated component mass flow rates
    )

    # high pressure pump, 2 degrees of freedom (efficiency and outlet pressure)
    m.fs.pump.efficiency_pump.fix(0.80)  # pump efficiency [-]
    m.fs.pump.control_volume.properties_out[0].pressure.fix(
        75e5
    )  # pump outlet pressure [Pa]

    # RO unit
    m.fs.RO.A_comp.fix(4.2e-12)  # membrane water permeability coefficient [m/s-Pa]
    m.fs.RO.B_comp.fix(3.5e-8)  # membrane salt permeability coefficient [m/s]
    m.fs.RO.channel_height.fix(1e-3)  # channel height in membrane stage [m]
    m.fs.RO.spacer_porosity.fix(0.97)  # spacer porosity in membrane stage [-]
    m.fs.RO.permeate.pressure[0].fix(101325)  # atmospheric pressure [Pa]
    # m.fs.RO.width.fix(5)  # stage width [m]
    # m.fs.RO.N_Re[0, 0].fix(500)
    m.fs.RO.velocity[0, 0].fix(0.15)
    m.fs.RO.recovery_vol_phase[0, "Liq"].fix(0.5)

    # energy recovery device, 2 degrees of freedom (efficiency and outlet pressure)
    m.fs.erd.efficiency_pump.fix(0.80)  # erd efficiency [-]
    m.fs.erd.control_volume.properties_out[0].pressure.fix(
        101325
    )  # atmospheric outlet pressure [Pa]

    # check degrees of freedom
    if degrees_of_freedom(m) != 0:
        raise RuntimeError(
            "The set_operating_conditions function resulted in {} "
            "degrees of freedom rather than 0. This error suggests "
            "that too many or not enough variables are fixed for a "
            "simulation.".format(degrees_of_freedom(m))
        )

    # for UI:
    export_variables(
        m.fs.pump, variables={"efficiency_pump": {"display_name": "pump efficiency"}}
    )
    export_variables(
        m.fs.RO,
        variables={
            "area": {"display_name": "membrane area"},
            "A_comp": {
                "display_name": "water perm coeff",
                "description": "membrane water permeability coefficient",
            },
            "B_comp": {
                "display_name": "salt perm coeff",
                "description": "membrane salt permeability coefficient",
            },
            "channel_height": {
                "display_name": "membrane channel height",
                "description": "channel height in membrane stage",
            },
            "spacer_porosity": {
                "display_name": "membrane spacer porosity",
                "description": "spacer porosity in membrane stage",
            },
            "width": {"display_name": "stage width", "description": "RO stage width"},
        },
    )
    export_variables(
        m.fs.RO.permeate,
        variables={
            "pressure": {
                "display_name": "atm pressure",
                "description": "atmospheric pressure",
            },
        },
    )
    export_variables(m.fs.RO.inlet, variables=["pressure"])


def solve(blk, solver=None, tee=False, check_termination=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    if check_termination:
        assert_optimal_termination(results)
    return results


def initialize_system(m):

    m.fs.feed.initialize()
    propagate_state(m.fs.s01)
    m.fs.pump.initialize()
    propagate_state(m.fs.s02)
    m.fs.RO.initialize()
    propagate_state(m.fs.s03)
    propagate_state(m.fs.s04)
    m.fs.erd.initialize()
    propagate_state(m.fs.s05)

    m.fs.costing.initialize()


def optimize_set_up(m):
    # objective
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)

    # unfix decision variables and add bounds
    # pump 1 and pump 2
    m.fs.pump.control_volume.properties_out[0].pressure.unfix()
    m.fs.pump.control_volume.properties_out[0].pressure.setlb(10e5)
    m.fs.pump.control_volume.properties_out[0].pressure.setub(85e5)
    m.fs.pump.deltaP.setlb(0)

    # RO
    # m.fs.RO.N_Re[0, 0].unfix()
    # m.fs.RO.N_Re.setlb(1)
    # m.fs.RO.N_Re.setub(1000)
    m.fs.RO.velocity[0, 0].unfix()
    m.fs.RO.velocity.setlb(0.01)
    m.fs.RO.velocity.setub(1)
    m.fs.RO.area.setlb(1)
    m.fs.RO.area.setub(150)

    # additional specifications
    m.fs.max_product_salinity = Param(
        initialize=500e-6, mutable=True, units=units.dimensionless
    )  # product TDS mass fraction [-]
    m.fs.max_pressure = Param(
        initialize=85e5, mutable=True, units=units.Pa
    )  # product TDS mass fraction [-]
    m.fs.minimum_water_flux = Param(
        initialize=1.0 / 3600.0, mutable=True, units=units.kg / units.m**2 / units.s
    )  # minimum water flux [kg/m2-s]

    # additional constraints
    m.fs.eq_product_quality = Constraint(
        expr=m.fs.product.properties[0].mass_frac_phase_comp["Liq", "TDS"]
        <= m.fs.max_product_salinity
    )
    iscale.constraint_scaling_transform(
        m.fs.eq_product_quality, 1e3
    )  # scaling constraint
    m.fs.eq_max_pressure = Constraint(
        expr=m.fs.RO.feed_side.properties[0, 0].pressure <= m.fs.max_pressure
    )
    iscale.constraint_scaling_transform(
        m.fs.eq_max_pressure, 1e-6
    )  # scaling constraint
    m.fs.eq_minimum_water_flux = Constraint(
        expr=m.fs.RO.flux_mass_phase_comp[0, 1, "Liq", "H2O"] >= m.fs.minimum_water_flux
    )

    # ---checking model---
    assert_degrees_of_freedom(m, 2)


def optimize(m, check_termination=True):
    # --solve---
    return solve(m, check_termination=check_termination)


def display_system(m):
    print("---system metrics---")
    feed_flow_mass = sum(
        m.fs.feed.flow_mass_phase_comp[0, "Liq", j].value for j in ["H2O", "TDS"]
    )
    feed_mass_frac_TDS = (
        m.fs.feed.flow_mass_phase_comp[0, "Liq", "TDS"].value / feed_flow_mass
    )
    print("Feed: %.2f kg/s, %.0f ppm" % (feed_flow_mass, feed_mass_frac_TDS * 1e6))

    prod_flow_mass = sum(
        m.fs.product.flow_mass_phase_comp[0, "Liq", j].value for j in ["H2O", "TDS"]
    )
    prod_mass_frac_TDS = (
        m.fs.product.flow_mass_phase_comp[0, "Liq", "TDS"].value / prod_flow_mass
    )
    print("Product: %.3f kg/s, %.0f ppm" % (prod_flow_mass, prod_mass_frac_TDS * 1e6))

    print(
        "Volumetric recovery: %.1f%%"
        % (value(m.fs.RO.recovery_vol_phase[0, "Liq"]) * 100)
    )
    print(
        "Water recovery: %.1f%%"
        % (value(m.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"]) * 100)
    )
    print(
        "Energy Consumption: %.1f kWh/m3"
        % value(m.fs.costing.specific_energy_consumption)
    )
    print("Levelized cost of water: %.2f $/m3" % value(m.fs.costing.LCOW))
    # for UI, return a result dict
    return {
        "Product": "%.3f kg/s, %.0f ppm" % (prod_flow_mass, prod_mass_frac_TDS * 1e6),
        "Volumetric recovery": "%.1f%%"
        % (value(m.fs.RO.recovery_vol_phase[0, "Liq"]) * 100),
        "Water recovery": "%.1f%%"
        % (value(m.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"]) * 100),
        "Energy Consumption": "%.1f kWh/m3"
        % value(m.fs.costing.specific_energy_consumption),
        "Levelized cost of water": "%.2f $/m3" % value(m.fs.costing.LCOW),
    }


def display_design(m):
    print("---decision variables---")
    print("Operating pressure %.1f bar" % (m.fs.RO.inlet.pressure[0].value / 1e5))
    print(
        "Membrane\narea %.1f m2\ninlet Reynolds %.1f, inlet velocity %.1f cm/s"
        % (
            m.fs.RO.area.value,
            m.fs.RO.N_Re[0, 0].value,
            m.fs.RO.velocity[0, 0].value * 100,
        )
    )

    print("---system variables---")
    print(
        "Pump\noutlet pressure: %.1f bar\npower %.2f kW"
        % (
            m.fs.pump.outlet.pressure[0].value / 1e5,
            m.fs.pump.work_mechanical[0].value / 1e3,
        )
    )
    print(
        "Membrane"
        "\naverage flux: %.1f LMH"
        "\npressure drop: %.1f bar"
        "\nmax interfacial conc %.1f ppm"
        % (
            value(m.fs.RO.flux_mass_phase_comp_avg[0, "Liq", "H2O"]) * 3600,
            m.fs.RO.deltaP[0].value / 1e5,
            m.fs.RO.feed_side.properties_interface[0, 1]
            .mass_frac_phase_comp["Liq", "TDS"]
            .value
            * 1e6,
        )
    )


def display_state(m):
    print("---state---")

    def print_state(s, b):
        flow_mass = sum(
            b.flow_mass_phase_comp[0, "Liq", j].value for j in ["H2O", "TDS"]
        )
        mass_frac_ppm = b.flow_mass_phase_comp[0, "Liq", "TDS"].value / flow_mass * 1e6
        pressure_bar = b.pressure[0].value / 1e5
        print(
            s
            + ": %.3f kg/s, %.0f ppm, %.1f bar"
            % (flow_mass, mass_frac_ppm, pressure_bar)
        )

    print_state("Feed      ", m.fs.feed.outlet)
    print_state("Pump out  ", m.fs.pump.outlet)
    print_state("RO perm   ", m.fs.RO.permeate)
    print_state("RO reten  ", m.fs.RO.retentate)


def export_ui_variables(fs):
    class Category:
        """Names for categories"""

        feed = "Feed"
        ts = "Treatment specification"
        perf = "Performance parameters"
        cost = "Cost parameters"

    export_variables(
        fs.feed.properties[0],
        variables={
            "flow_vol_phase['Liq']": {
                "display_name": "Volumetric flowrate",
                "to_units": units.m**3 / units.hr,
                "category": Category.feed,
            },
            "mass_frac_phase_comp['Liq', 'TDS']": {
                "display_name": "Salinity",
                "display_units": "ppm",
                "scale_units": 1e6,
                "category": Category.feed,
            },
            "temperature": {
                "display_name": "Temperature",
                "to_units": units.K,
                "category": Category.feed,
            },
            "pressure": {
                "display_name": "Pressure",
                "to_units": units.bar,
                "category": Category.feed,
            },
        },
    )
    export_variables(
        fs.RO,
        variables={
            "recovery_vol_phase[0, 'Liq']": {
                "display_units": "%",
                "scale_units": 100,
                "category": Category.ts,
            }
        },
    )
    export_variables(
        fs,
        variables={
            "max_product_salinity": {
                "display_name": "Maximum product salinity",
                "scale_units": 1e6,
                "display_units": "ppm",
                "category": Category.ts,
            },
            "max_pressure": {
                "display_name": "Maximum allowable pressure",
                "to_units": units.bar,
                "category": Category.ts,
            },
        },
    )


def display_ui_input(m):
    return {
        "Feed": {
            "Volumetric flowrate": (
                round(
                    value(
                        units.convert(
                            m.fs.feed.properties[0].flow_vol_phase["Liq"],
                            to_units=units.m**3 / units.hr,
                        )
                    ),
                    2,
                ),
                "m3/h",
            ),
            "Salinity": (
                round(
                    value(
                        m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"] * 1e6
                    ),
                    0,
                ),
                "ppm",
            ),
            "Temperature": (
                round(
                    value(m.fs.feed.properties[0].temperature),
                    0,
                ),
                "K",
            ),
            "Pressure": (
                round(
                    value(
                        units.convert(
                            m.fs.feed.properties[0].pressure, to_units=units.bar
                        )
                    ),
                    1,
                ),
                "bar",
            ),
        },
        "Treatment specification": {
            "Recovery": (
                round(
                    value(m.fs.RO.recovery_vol_phase[0, "Liq"]) * 100,
                    1,
                ),
                "%",
            ),
            "Maximum product salinity": (
                round(
                    value(m.fs.max_product_salinity) * 1e6,
                    0,
                ),
                "ppm",
            ),
            "Maximum allowable pressure": (
                round(
                    value(units.convert(m.fs.max_pressure, to_units=units.bar)),
                    1,
                ),
                "bar",
            ),
        },
        "Performance parameters": {
            "Pump efficiency": (
                round(
                    value(m.fs.pump.efficiency_pump[0]) * 100,
                    1,
                ),
                "%",
            ),
            "ERD efficiency": (
                round(
                    value(m.fs.erd.efficiency_pump[0]) * 100,
                    1,
                ),
                "%",
            ),
            "Water permeability coeff": (
                round(
                    value(
                        units.convert(
                            m.fs.RO.A_comp[0, "H2O"],
                            to_units=units.mm / units.hr / units.bar,
                        )
                    ),
                    2,
                ),
                "L/(m2-h-bar)",
            ),
            "Salt permeability coeff": (
                round(
                    value(
                        units.convert(
                            m.fs.RO.B_comp[0, "TDS"], to_units=units.mm / units.hr
                        )
                    ),
                    2,
                ),
                "L/(m2-h)",
            ),
            "RO channel height": (
                round(
                    value(units.convert(m.fs.RO.channel_height, to_units=units.mm)),
                    1,
                ),
                "mm",
            ),
            "RO spacer porosity": (
                round(
                    value(m.fs.RO.spacer_porosity) * 100,
                    1,
                ),
                "%",
            ),
        },
        "Cost parameters": {
            "Electricity cost": (
                round(
                    value(m.fs.costing.electricity_base_cost),
                    3,
                ),
                "$/kWh",
            ),
            "Membrane cost": (
                round(
                    value(m.fs.costing.reverse_osmosis_membrane_cost),
                    1,
                ),
                "$/m2",
            ),
            "Pump cost": (
                round(
                    value(
                        units.convert(
                            m.fs.costing.high_pressure_pump_cost,
                            to_units=units.USD_2018 / units.kW,
                        )
                    ),
                    0,
                ),
                "$/kW",
            ),
            "ERD cost": (
                round(
                    value(
                        units.convert(
                            m.fs.costing.erd_pressure_exchanger_cost,
                            to_units=units.USD_2018 / (units.m**3 / units.hr),
                        )
                    ),
                    0,
                ),
                "$/(m3/h)",
            ),
            "Load factor": (
                round(
                    value(m.fs.costing.load_factor) * 100,
                    1,
                ),
                "%",
            ),
            "Capital annualization factor": (
                round(
                    value(m.fs.costing.factor_capital_annualization) * 100,
                    1,
                ),
                "%/year",
            ),
            "Membrane replacement factor": (
                round(
                    value(m.fs.costing.factor_membrane_replacement) * 100,
                    1,
                ),
                "%/year",
            ),
        },
    }


def display_ui_output(m):
    return {
        "System metrics": {
            "Recovery": (
                round(
                    value(m.fs.RO.recovery_vol_phase[0, "Liq"]) * 100,
                    1,
                ),
                "%",
            ),
            "Specific energy consumption": (
                round(
                    value(m.fs.costing.specific_energy_consumption),
                    2,
                ),
                "kWh/m3",
            ),
            "Levelized cost of water": (
                round(
                    value(m.fs.costing.LCOW),
                    2,
                ),
                "$/m3",
            ),
        },
        "Feed": {
            "Volumetric flowrate": (
                round(
                    value(
                        units.convert(
                            m.fs.feed.properties[0].flow_vol_phase["Liq"],
                            to_units=units.m**3 / units.hr,
                        )
                    ),
                    2,
                ),
                "m3/h",
            ),
            "Salinity": (
                round(
                    value(
                        m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "TDS"] * 1e6
                    ),
                    0,
                ),
                "ppm",
            ),
            "Temperature": (
                round(
                    value(m.fs.feed.properties[0].temperature),
                    0,
                ),
                "K",
            ),
            "Pressure": (
                round(
                    value(
                        units.convert(
                            m.fs.feed.properties[0].pressure, to_units=units.bar
                        )
                    ),
                    1,
                ),
                "bar",
            ),
        },
        "Product": {
            "Volumetric flowrate": (
                round(
                    value(
                        units.convert(
                            m.fs.product.properties[0].flow_vol_phase["Liq"],
                            to_units=units.m**3 / units.hr,
                        )
                    ),
                    2,
                ),
                "m3/h",
            ),
            "Salinity": (
                round(
                    value(
                        m.fs.product.properties[0].mass_frac_phase_comp["Liq", "TDS"]
                        * 1e6
                    ),
                    0,
                ),
                "ppm",
            ),
            "Temperature": (
                round(
                    value(m.fs.product.properties[0].temperature),
                    0,
                ),
                "K",
            ),
            "Pressure": (
                round(
                    value(
                        units.convert(
                            m.fs.product.properties[0].pressure, to_units=units.bar
                        )
                    ),
                    1,
                ),
                "bar",
            ),
        },
        "Disposal": {
            "Volumetric flowrate": (
                round(
                    value(
                        units.convert(
                            m.fs.disposal.properties[0].flow_vol_phase["Liq"],
                            to_units=units.m**3 / units.hr,
                        )
                    ),
                    2,
                ),
                "m3/h",
            ),
            "Salinity": (
                round(
                    value(
                        m.fs.disposal.properties[0].mass_frac_phase_comp["Liq", "TDS"]
                        * 1e6
                    ),
                    0,
                ),
                "ppm",
            ),
            "Temperature": (
                round(
                    value(m.fs.disposal.properties[0].temperature),
                    0,
                ),
                "K",
            ),
            "Pressure": (
                round(
                    value(
                        units.convert(
                            m.fs.disposal.properties[0].pressure, to_units=units.bar
                        )
                    ),
                    1,
                ),
                "bar",
            ),
        },
        "Decision variables": {
            "Operating pressure": (
                round(
                    value(
                        units.convert(
                            m.fs.pump.control_volume.properties_out[0].pressure,
                            to_units=units.bar,
                        )
                    ),
                    1,
                ),
                "bar",
            ),
            "Membrane area": (
                round(
                    value(m.fs.RO.area),
                    1,
                ),
                "m2",
            ),
            "Inlet crossflow velocity": (
                round(
                    value(
                        units.convert(
                            m.fs.RO.velocity[0, 0], to_units=units.cm / units.s
                        )
                    ),
                    1,
                ),
                "cm/s",
            ),
        },
        "System variables": {
            "Pump power": (
                round(
                    value(
                        units.convert(
                            m.fs.pump.work_mechanical[0],
                            to_units=units.kW,
                        )
                    ),
                    1,
                ),
                "kW",
            ),
            "ERD power": (
                round(
                    value(
                        units.convert(
                            m.fs.erd.work_mechanical[0],
                            to_units=units.kW,
                        )
                        * -1
                    ),
                    1,
                ),
                "kW",
            ),
            "Average water flux": (
                round(
                    value(
                        units.convert(
                            m.fs.RO.flux_mass_phase_comp_avg[0, "Liq", "H2O"]
                            / (1000 * units.kg / units.m**3),
                            to_units=units.mm / units.hr,
                        )
                    ),
                    1,
                ),
                "L/(m2-h)",
            ),
            "Pressure drop": (
                round(
                    value(
                        units.convert(
                            m.fs.RO.deltaP[0],
                            to_units=units.bar,
                        )
                        * -1
                    ),
                    2,
                ),
                "bar",
            ),
            "Max interfacial salinity": (
                round(
                    value(
                        m.fs.RO.feed_side.properties_interface[
                            0, 1
                        ].mass_frac_phase_comp["Liq", "TDS"]
                    )
                    * 1e6,
                    0,
                ),
                "ppm",
            ),
        },
    }


if __name__ == "__main__":
    main()
