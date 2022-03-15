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
from pyomo.environ import (ConcreteModel,
                           value,
                           Constraint,
                           Expression,
                           Objective,
                           Param,
                           TransformationFactory,
                           units as pyunits,
                           assert_optimal_termination)
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import (solve_indexed_blocks,
                                            propagate_state)
from idaes.generic_models.unit_models import Product, Feed, Mixer
from idaes.generic_models.unit_models.mixer import MomentumMixingType
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

import watertap.property_models.NaCl_prop_pack as props
from watertap.unit_models.reverse_osmosis_0D import (ReverseOsmosis0D,
                                                       ConcentrationPolarizationType,
                                                       MassTransferCoefficient,
                                                       PressureChangeType)
from watertap.unit_models.pump_isothermal import Pump
from watertap.core.util.initialization import assert_degrees_of_freedom
import watertap.examples.flowsheets.high_pressure_RO.financials as financials


def main():
    # set up solver
    solver = get_solver()

    # build, set, and initialize
    m = build()
    specify_model(m, solver=solver)
    initialize_model(m, solver=solver)

    # simulate and display
    solve(m, solver=solver)
    print('\n***---Simulation results---***')
    display_results(m)

    # optimize and display
    set_up_optimization(m)
    optimize(m, solver=solver)
    print('\n***---Optimization results---***')
    display_results(m)

def build():
    # flowsheet set up
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={'dynamic': False})
    m.fs.properties = props.NaClParameterBlock()

    # unit models
    m.fs.feed = Feed(default={'property_package': m.fs.properties})
    m.fs.P1 = Pump(default={'property_package': m.fs.properties})
    m.fs.P2 = Pump(default={'property_package': m.fs.properties})
    m.fs.RO1 = ReverseOsmosis0D(default={
        "property_package": m.fs.properties,
        "has_pressure_change": True,
        "pressure_change_type": PressureChangeType.calculated,
        "mass_transfer_coefficient": MassTransferCoefficient.calculated,
        "concentration_polarization_type": ConcentrationPolarizationType.calculated,
    })
    m.fs.RO2 = ReverseOsmosis0D(default={
        "property_package": m.fs.properties,
        "has_pressure_change": True,
        "pressure_change_type": PressureChangeType.calculated,
        "mass_transfer_coefficient": MassTransferCoefficient.calculated,
        "concentration_polarization_type": ConcentrationPolarizationType.calculated,
    })
    m.fs.M1 = Mixer(default={
        "property_package": m.fs.properties,
        "momentum_mixing_type": MomentumMixingType.equality,
        "inlet_list": ['RO1', 'RO2']})
    m.fs.product = Product(default={'property_package': m.fs.properties})
    m.fs.disposal = Product(default={'property_package': m.fs.properties})

    # build old costing
    financials.add_costing_param_block(m.fs)

    m.fs.P1.get_costing(module=financials, pump_type="High pressure")
    m.fs.P2.get_costing(module=financials, pump_type="High pressure")
    m.fs.RO1.get_costing(module=financials)
    m.fs.RO2.get_costing(module=financials)

    product_flow_vol_total = m.fs.product.properties[0].flow_vol
    m.fs.annual_water_production = Expression(
        expr=pyunits.convert(product_flow_vol_total, to_units=pyunits.m ** 3 / pyunits.year)
             * m.fs.costing_param.load_factor)
    pump_power_total = m.fs.P1.work_mechanical[0] + m.fs.P2.work_mechanical[0]
    m.fs.specific_energy_consumption = Expression(
        expr=pyunits.convert(pump_power_total, to_units=pyunits.kW)
             / pyunits.convert(product_flow_vol_total, to_units=pyunits.m ** 3 / pyunits.hr))

    financials.get_system_costing(m.fs)

    iscale.set_scaling_factor(m.fs.P1.costing.purchase_cost, 1)
    iscale.set_scaling_factor(m.fs.P2.costing.purchase_cost, 1)

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.P1.inlet)
    m.fs.s02 = Arc(source=m.fs.P1.outlet, destination=m.fs.RO1.inlet)
    m.fs.s03 = Arc(source=m.fs.RO1.permeate, destination=m.fs.M1.RO1)
    m.fs.s04 = Arc(source=m.fs.RO1.retentate, destination=m.fs.P2.inlet)
    m.fs.s05 = Arc(source=m.fs.P2.outlet, destination=m.fs.RO2.inlet)
    m.fs.s06 = Arc(source=m.fs.RO2.permeate, destination=m.fs.M1.RO2)
    m.fs.s07 = Arc(source=m.fs.RO2.retentate, destination=m.fs.disposal.inlet)
    m.fs.s08 = Arc(source=m.fs.M1.outlet, destination=m.fs.product.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    # set default property values
    m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
    m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'NaCl'))
    # set unit model values
    iscale.set_scaling_factor(m.fs.P1.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.P2.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.RO1.area, 1)
    iscale.set_scaling_factor(m.fs.RO2.area, 1)
    # touch properties used in specifying and initializing the model
    m.fs.feed.properties[0].flow_vol_phase['Liq']
    m.fs.feed.properties[0].mass_frac_phase_comp['Liq', 'NaCl']
    # calculate and propagate scaling factors
    iscale.calculate_scaling_factors(m)

    return m

def specify_model(m, solver=None):
    if solver is None:
        solver = get_solver()

    # ---specifications---
    # feed
    # state variables
    m.fs.feed.properties[0].pressure.fix(101325)  # feed pressure [Pa]
    m.fs.feed.properties[0].temperature.fix(273.15 + 25)  # feed temperature [K]
    # properties (cannot be fixed for initialization routines, must calculate the state variables)
    m.fs.feed.properties.calculate_state(
        var_args={('flow_vol_phase', 'Liq'): 1e-3,  # feed volumetric flow rate [m3/s]
                  ('mass_frac_phase_comp', ('Liq', 'NaCl')): 0.035},  # feed NaCl mass fraction [-]
        hold_state=True,  # fixes the calculated component mass flow rates
    )

    # pump 1, 2 degrees of freedom (efficiency and outlet pressure)
    m.fs.P1.efficiency_pump.fix(0.80)  # pump efficiency [-]
    m.fs.P1.control_volume.properties_out[0].pressure.fix(65e5)

    # pump 2, 2 degrees of freedom (efficiency and outlet pressure)
    m.fs.P2.efficiency_pump.fix(0.80)  # pump efficiency [-]
    m.fs.P2.control_volume.properties_out[0].pressure.fix(120e5)

    # RO 1, 7 degrees of freedom
    m.fs.RO1.A_comp.fix(4.2e-12)  # membrane water permeability coefficient [m/s-Pa]
    m.fs.RO1.B_comp.fix(3.5e-8)  # membrane salt permeability coefficient [m/s]
    m.fs.RO1.channel_height.fix(1e-3)  # channel height in membrane stage [m]
    m.fs.RO1.spacer_porosity.fix(0.97)  # spacer porosity in membrane stage [-]
    m.fs.RO1.permeate.pressure[0].fix(101325)  # atmospheric pressure [Pa]
    m.fs.RO1.area.fix(100)
    m.fs.RO1.width.fix(20)  # stage width [m]

    # RO 2, 7 degrees of freedom
    m.fs.RO2.A_comp.fix(4.2e-12)  # membrane water permeability coefficient [m/s-Pa]
    m.fs.RO2.B_comp.fix(3.5e-8)  # membrane salt permeability coefficient [m/s]
    m.fs.RO2.channel_height.fix(1e-3)  # channel height in membrane stage [m]
    m.fs.RO2.spacer_porosity.fix(0.97)  # spacer porosity in membrane stage [-]
    # m.fs.RO2.permeate.pressure[0].fix(101325)  # atmospheric pressure [Pa]  # mixer has equality constraint
    m.fs.RO2.area.fix(50)
    m.fs.RO2.width.fix(10)  # stage width [m]

    # check degrees of freedom
    if degrees_of_freedom(m) != 0:
        raise RuntimeError("The specify_model function resulted in {} "
                           "degrees of freedom rather than 0. This error suggests "
                           "that too many or not enough variables are fixed for a "
                           "simulation.".format(degrees_of_freedom(m)))


def solve(blk, solver=None, tee=False):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    assert_optimal_termination(results)


def initialize_model(m, solver=None):
    if solver is None:
        solver = get_solver()
    optarg = solver.options

    m.fs.feed.initialize(optarg=optarg)
    propagate_state(m.fs.s01)
    m.fs.P1.initialize(optarg=optarg)
    propagate_state(m.fs.s02)
    m.fs.RO1.initialize(optarg=optarg)
    propagate_state(m.fs.s03)
    propagate_state(m.fs.s04)
    m.fs.P2.initialize(optarg=optarg)
    propagate_state(m.fs.s05)
    m.fs.RO2.permeate.pressure[0].fix(101325)
    m.fs.RO2.initialize(optarg=optarg)
    m.fs.RO2.permeate.pressure[0].unfix()
    propagate_state(m.fs.s06)
    propagate_state(m.fs.s07)
    m.fs.M1.initialize(optarg=optarg)
    propagate_state(m.fs.s08)
    m.fs.product.initialize(optarg=optarg)
    m.fs.disposal.initialize(optarg=optarg)


def set_up_optimization(m):
    # objective
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)

    # unfix decision variables and add bounds
    # pump 1 and pump 2
    m.fs.P1.control_volume.properties_out[0].pressure.unfix()
    m.fs.P1.control_volume.properties_out[0].pressure.setlb(10e5)
    m.fs.P1.control_volume.properties_out[0].pressure.setub(80e5)
    m.fs.P1.deltaP.setlb(0)
    m.fs.P2.control_volume.properties_out[0].pressure.unfix()
    m.fs.P2.control_volume.properties_out[0].pressure.setlb(50e5)
    m.fs.P2.control_volume.properties_out[0].pressure.setub(150e5)
    m.fs.P2.deltaP.setlb(0)

    # RO
    m.fs.RO1.area.unfix()
    m.fs.RO1.area.setlb(1)
    m.fs.RO1.area.setub(150)
    m.fs.RO2.area.unfix()
    m.fs.RO2.area.setlb(1)
    m.fs.RO2.area.setub(150)

    # additional specifications
    m.fs.product_salinity = Param(initialize=1000e-6, mutable=True)  # product NaCl mass fraction [-]
    m.fs.minimum_water_flux = Param(initialize=1./3600., mutable=True)  # minimum water flux [kg/m2-s]
    m.fs.product_recovery = Param(initialize=0.73, mutable=True)

    # additional constraints
    m.fs.eq_product_quality = Constraint(
        expr=m.fs.product.properties[0].mass_frac_phase_comp['Liq', 'NaCl'] <= m.fs.product_salinity)
    iscale.constraint_scaling_transform(m.fs.eq_product_quality, 1e3)  # scaling constraint
    m.fs.eq_minimum_water_flux_1 = Constraint(
        expr=m.fs.RO1.flux_mass_io_phase_comp[0, 'out', 'Liq', 'H2O'] >= m.fs.minimum_water_flux)
    m.fs.eq_minimum_water_flux_2 = Constraint(
        expr=m.fs.RO2.flux_mass_io_phase_comp[0, 'out', 'Liq', 'H2O'] >= m.fs.minimum_water_flux)
    m.fs.eq_product_recovery = Constraint(
        expr=m.fs.product.properties[0].flow_vol ==
             m.fs.product_recovery * m.fs.feed.properties[0].flow_vol)

    # ---checking model---
    assert_degrees_of_freedom(m, 3)

def optimize(m, solver=None):
    # --solve---
    solve(m, solver=solver)

def display_results(m):
    print('---system metrics---')
    feed_flow_mass = sum(m.fs.feed.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
    feed_mass_frac_NaCl = m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / feed_flow_mass
    print('Feed: %.2f kg/s, %.0f ppm' % (feed_flow_mass, feed_mass_frac_NaCl * 1e6))

    prod_flow_mass = sum(m.fs.product.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
    prod_mass_frac_NaCl = m.fs.product.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / prod_flow_mass
    print('Product: %.3f kg/s, %.0f ppm' % (prod_flow_mass, prod_mass_frac_NaCl * 1e6))

    disp_flow_mass = sum(m.fs.disposal.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
    disp_mass_frac_NaCl = m.fs.disposal.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / disp_flow_mass
    print('Disposal: %.3f kg/s, %.0f ppm' % (disp_flow_mass, disp_mass_frac_NaCl * 1e6))

    print('Volumetric recovery: %.1f%%' % (value(m.fs.product.properties[0].flow_vol
                                                 / m.fs.feed.properties[0].flow_vol) * 100))
    print('Energy Consumption: %.1f kWh/m3' % value(m.fs.specific_energy_consumption))
    print('Levelized cost of water: %.2f $/m3' % value(m.fs.costing.LCOW))

    print('---decision variables---')
    print('Operating pressure: %.1f and %.1f bar' %
          (m.fs.RO1.inlet.pressure[0].value / 1e5, m.fs.RO2.inlet.pressure[0].value / 1e5))
    print('Membrane area %.1f and %.1f m2' %
          (m.fs.RO1.area.value, m.fs.RO2.area.value))


if __name__ == "__main__":
    main()
