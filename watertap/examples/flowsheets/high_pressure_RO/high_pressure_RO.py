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
                                            propagate_state,
                                            fix_state_vars,
                                            revert_state_vars)
from idaes.generic_models.unit_models import Product, Feed, Mixer, Translator
from idaes.generic_models.unit_models.mixer import MomentumMixingType
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

import watertap.property_models.seawater_prop_pack as props
from watertap.unit_models.reverse_osmosis_0D import (ReverseOsmosis0D,
                                                       ConcentrationPolarizationType,
                                                       MassTransferCoefficient,
                                                       PressureChangeType)
from watertap.unit_models.pump_isothermal import Pump
from watertap.core.util.initialization import assert_degrees_of_freedom
from watertap.examples.flowsheets.full_treatment_train.util import solve_block
import watertap.examples.flowsheets.high_pressure_RO.financials as financials
import watertap.examples.flowsheets.high_pressure_RO.components.feed_cases as feed_cases
from watertap.core.util.infeasible import print_close_to_bounds, print_infeasible_bounds, print_infeasible_constraints


def main():

    # build, set, and initialize
    m = build()
    specify_model(m)
    initialize_model(m)

    # simulate and display
    solve(m)
    print('\n***---Simulation results---***')
    display_report(m)
    display_results(m)

    # optimize and display
    set_up_optimization(m)
    optimize(m)
    print('\n***---Optimization results---***')
    display_report(m)
    display_results(m)

    # optimize and display
    m.fs.product_recovery = 0.75
    solve(m)
    print('\n***---Optimization results 2---***')
    display_report(m)
    display_results(m)

def build(case='seawater'):
    # flowsheet set up
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={'dynamic': False})
    m.fs.prop_desal = props.SeawaterParameterBlock()

    # feed and disposal based on the case
    feed_cases.build_specified_feed(m, case=case)
    m.fs.disposal = Product(default={'property_package': m.fs.prop_feed})

    # unit models
    m.fs.P1 = Pump(default={'property_package': m.fs.prop_desal})
    m.fs.P2 = Pump(default={'property_package': m.fs.prop_desal})
    m.fs.RO1 = ReverseOsmosis0D(default={
        "property_package": m.fs.prop_desal,
        "has_pressure_change": True,
        "pressure_change_type": PressureChangeType.calculated,
        "mass_transfer_coefficient": MassTransferCoefficient.calculated,
        "concentration_polarization_type": ConcentrationPolarizationType.calculated,
    })
    m.fs.RO2 = ReverseOsmosis0D(default={
        "property_package": m.fs.prop_desal,
        "has_pressure_change": True,
        "pressure_change_type": PressureChangeType.calculated,
        "mass_transfer_coefficient": MassTransferCoefficient.calculated,
        "concentration_polarization_type": ConcentrationPolarizationType.calculated,
    })
    m.fs.M1 = Mixer(default={
        "property_package": m.fs.prop_desal,
        "momentum_mixing_type": MomentumMixingType.equality,
        "inlet_list": ['RO1', 'RO2']})
    m.fs.product = Product(default={'property_package': m.fs.prop_desal})

    # translator blocks
    m.fs.tb_feed_desal = Translator(
        default={"inlet_property_package": m.fs.prop_feed,
                 "outlet_property_package": m.fs.prop_desal})

    @m.fs.tb_feed_desal.Constraint()
    def eq_flow_vol(blk):
        return blk.properties_in[0].flow_vol == blk.properties_out[0].flow_vol_phase['Liq']

    @m.fs.tb_feed_desal.Constraint()
    def eq_flow_mass_TDS(blk):
        return (sum(blk.properties_in[0].conc_mass_comp[j]
                    for j in blk.properties_in[0].params.solute_set)
                == blk.properties_out[0].conc_mass_phase_comp['Liq', 'TDS'])

    @m.fs.tb_feed_desal.Constraint()
    def eq_pressure(blk):
        return blk.properties_in[0].pressure == blk.properties_out[0].pressure

    @m.fs.tb_feed_desal.Constraint()
    def eq_temperature(blk):
        return blk.properties_in[0].temperature == blk.properties_out[0].temperature

    m.fs.tb_desal_disposal = Translator(
        default={"inlet_property_package": m.fs.prop_desal,
                 "outlet_property_package": m.fs.prop_feed})

    @m.fs.tb_desal_disposal.Constraint()
    def eq_flow_vol(blk):
        return blk.properties_in[0].flow_vol_phase['Liq'] == blk.properties_out[0].flow_vol

    @m.fs.tb_desal_disposal.Constraint(m.fs.prop_feed.solute_set)
    def eq_flow_mass_comp(blk, j):
        flow_mass_out = blk.properties_out[0].conc_mass_comp[j] * blk.properties_out[0].flow_vol
        flow_mass_in = m.fs.feed.properties[0].flow_vol * m.fs.feed.properties[0].conc_mass_comp[j]
        if j == 'NAION':
            flow_mass_in -= m.fs.product.properties[0].flow_mass_phase_comp['Liq', 'TDS'] * 22.99/(22.99 + 35.45)
        elif j == 'CLION':
            flow_mass_in -= m.fs.product.properties[0].flow_mass_phase_comp['Liq', 'TDS'] * 35.45/(22.99 + 35.45)
        return flow_mass_in == flow_mass_out

    @m.fs.tb_desal_disposal.Constraint()
    def eq_pressure(blk):
        return blk.properties_in[0].pressure == blk.properties_out[0].pressure

    @m.fs.tb_desal_disposal.Constraint()
    def eq_temperature(blk):
        return blk.properties_in[0].temperature == blk.properties_out[0].temperature

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
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.tb_feed_desal.inlet)
    m.fs.s02 = Arc(source=m.fs.tb_feed_desal.outlet, destination=m.fs.P1.inlet)
    m.fs.s03 = Arc(source=m.fs.P1.outlet, destination=m.fs.RO1.inlet)
    m.fs.s04 = Arc(source=m.fs.RO1.permeate, destination=m.fs.M1.RO1)
    m.fs.s05 = Arc(source=m.fs.RO1.retentate, destination=m.fs.P2.inlet)
    m.fs.s06 = Arc(source=m.fs.P2.outlet, destination=m.fs.RO2.inlet)
    m.fs.s07 = Arc(source=m.fs.RO2.permeate, destination=m.fs.M1.RO2)
    m.fs.s08 = Arc(source=m.fs.RO2.retentate, destination=m.fs.tb_desal_disposal.inlet)
    m.fs.s09 = Arc(source=m.fs.tb_desal_disposal.outlet, destination=m.fs.disposal.inlet)
    m.fs.s10 = Arc(source=m.fs.M1.outlet, destination=m.fs.product.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    # set default property values
    m.fs.prop_desal.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
    m.fs.prop_desal.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'TDS'))
    # set unit model values
    iscale.set_scaling_factor(m.fs.P1.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.P2.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.RO1.area, 1)
    iscale.set_scaling_factor(m.fs.RO2.area, 1)
    # calculate and propagate scaling factors
    iscale.calculate_scaling_factors(m)

    return m

def specify_model(m):

    # ---specifications---

    # pump 1, 2 degrees of freedom (efficiency and outlet pressure)
    m.fs.P1.efficiency_pump.fix(0.80)  # pump efficiency [-]
    m.fs.P1.control_volume.properties_out[0].pressure.fix(65e5)

    # pump 2, 2 degrees of freedom (efficiency and outlet pressure)
    m.fs.P2.efficiency_pump.fix(0.80)  # pump efficiency [-]
    m.fs.P2.control_volume.properties_out[0].pressure.fix(100e5)

    # RO 1, 7 degrees of freedom
    m.fs.RO1.A_comp.fix(4.2e-12)  # membrane water permeability coefficient [m/s-Pa]
    m.fs.RO1.B_comp.fix(3.5e-8)  # membrane salt permeability coefficient [m/s]
    m.fs.RO1.channel_height.fix(1e-3)  # channel height in membrane stage [m]
    m.fs.RO1.spacer_porosity.fix(0.97)  # spacer porosity in membrane stage [-]
    m.fs.RO1.permeate.pressure[0].fix(101325)  # atmospheric pressure [Pa]
    m.fs.RO1.area.fix(100)
    m.fs.RO1.width.fix(10)  # stage width [m]

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


def solve(blk, tee=False):
    solver = get_solver()
    results = solver.solve(blk, tee=tee)
    print_infeasible_constraints(blk)
    print_infeasible_bounds(blk)
    print_close_to_bounds(blk)
    assert_optimal_termination(results)


def initialize_model(m):

    propagate_state(m.fs.s01)
    flags = fix_state_vars(m.fs.tb_feed_desal.properties_in)
    solve(m.fs.tb_feed_desal)
    revert_state_vars(m.fs.tb_feed_desal.properties_in, flags)
    propagate_state(m.fs.s02)
    m.fs.P1.initialize()
    propagate_state(m.fs.s03)
    m.fs.RO1.initialize()
    propagate_state(m.fs.s04)
    propagate_state(m.fs.s05)
    m.fs.P2.initialize()
    propagate_state(m.fs.s06)
    m.fs.RO2.permeate.pressure[0].fix(101325)
    m.fs.RO2.initialize()
    m.fs.RO2.permeate.pressure[0].unfix()
    propagate_state(m.fs.s07)
    m.fs.M1.initialize()
    propagate_state(m.fs.s08)
    # flags = fix_state_vars(m.fs.tb_desal_disposal.properties_in)
    # solve(m.fs.tb_desal_disposal)
    # revert_state_vars(m.fs.tb_desal_disposal.properties_in, flags)
    propagate_state(m.fs.s09)
    propagate_state(m.fs.s10)
    m.fs.product.initialize()

def set_up_optimization(m):
    # objective
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)

    # unfix decision variables and add bounds
    # pump 1 and pump 2
    m.fs.P1.control_volume.properties_out[0].pressure.unfix()
    m.fs.P1.control_volume.properties_out[0].pressure.setlb(10e5)
    m.fs.P1.control_volume.properties_out[0].pressure.setub(85e5)
    m.fs.P1.deltaP.setlb(0)
    m.fs.P2.control_volume.properties_out[0].pressure.unfix()
    m.fs.P2.control_volume.properties_out[0].pressure.setlb(50e5)
    m.fs.P2.control_volume.properties_out[0].pressure.setub(200e5)
    m.fs.P2.deltaP.setlb(0)

    # RO
    m.fs.RO1.area.unfix()
    m.fs.RO1.area.setlb(1)
    m.fs.RO1.area.setub(150)
    m.fs.RO1.width.unfix()
    m.fs.RO1.width.setlb(1)
    m.fs.RO1.width.setub(75)
    m.fs.RO2.area.unfix()
    m.fs.RO2.area.setlb(1)
    m.fs.RO2.area.setub(150)
    m.fs.RO2.width.unfix()
    m.fs.RO2.width.setlb(1)
    m.fs.RO2.width.setub(75)

    # additional specifications
    m.fs.maximum_product_salinity = Param(initialize=1000e-6, mutable=True)  # product TDS mass fraction [-]
    m.fs.minimum_water_flux = Param(initialize=1./3600., mutable=True)  # minimum water flux [kg/m2-s]
    m.fs.product_recovery = Param(initialize=0.73, mutable=True)

    # additional constraints
    m.fs.eq_product_quality = Constraint(
        expr=m.fs.product.properties[0].mass_frac_phase_comp['Liq', 'TDS'] <= m.fs.maximum_product_salinity)
    iscale.constraint_scaling_transform(m.fs.eq_product_quality, 1e3)  # scaling constraint
    m.fs.eq_minimum_water_flux_1 = Constraint(
        expr=m.fs.RO1.flux_mass_io_phase_comp[0, 'out', 'Liq', 'H2O'] >= m.fs.minimum_water_flux)
    m.fs.eq_minimum_water_flux_2 = Constraint(
        expr=m.fs.RO2.flux_mass_io_phase_comp[0, 'out', 'Liq', 'H2O'] >= m.fs.minimum_water_flux)
    m.fs.eq_product_recovery = Constraint(
        expr=m.fs.product.properties[0].flow_vol ==
             m.fs.product_recovery * m.fs.feed.properties[0].flow_vol)

    # ---checking model---
    assert_degrees_of_freedom(m, 5)


def optimize(m, check_termination=True):
    return solve_block(m, tee=False, fail_flag=check_termination)


def display_results(m):
    print('---system metrics---')
    feed_flow_mass = sum(m.fs.tb_feed_desal.properties_out[0].flow_mass_phase_comp['Liq', j].value
                         for j in ['H2O', 'TDS'])
    feed_mass_frac_TDS = (m.fs.tb_feed_desal.properties_out[0].flow_mass_phase_comp['Liq', 'TDS'].value
                          / feed_flow_mass)
    print('Feed: %.2f kg/s, %.0f ppm' % (feed_flow_mass, feed_mass_frac_TDS * 1e6))

    prod_flow_mass = sum(m.fs.product.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'TDS'])
    prod_mass_frac_TDS = m.fs.product.flow_mass_phase_comp[0, 'Liq', 'TDS'].value / prod_flow_mass
    print('Product: %.3f kg/s, %.0f ppm' % (prod_flow_mass, prod_mass_frac_TDS * 1e6))

    disp_flow_mass = sum(m.fs.tb_desal_disposal.properties_in[0].flow_mass_phase_comp['Liq', j].value
                         for j in ['H2O', 'TDS'])
    disp_mass_frac_TDS = (m.fs.tb_desal_disposal.properties_in[0].flow_mass_phase_comp['Liq', 'TDS'].value
                          / disp_flow_mass)
    print('Disposal: %.3f kg/s, %.0f ppm' % (disp_flow_mass, disp_mass_frac_TDS * 1e6))

    print('Volumetric recovery: %.1f%%' % (value(m.fs.product.properties[0].flow_vol
                                                 / m.fs.feed.properties[0].flow_vol) * 100))
    print('Energy Consumption: %.1f kWh/m3' % value(m.fs.specific_energy_consumption))
    print('Levelized cost of water: %.2f $/m3' % value(m.fs.costing.LCOW))

    print('---decision variables---')
    print('Operating pressure: %.1f and %.1f bar' %
          (m.fs.RO1.inlet.pressure[0].value / 1e5, m.fs.RO2.inlet.pressure[0].value / 1e5))
    print('Membrane area %.1f and %.1f m2' %
          (m.fs.RO1.area.value, m.fs.RO2.area.value))


def display_report(m):
    m.fs.feed.report()
    # m.fs.tb_feed_desal.report()
    # m.fs.P1.report()
    # m.fs.RO1.report()
    # m.fs.P2.report()
    # m.fs.RO2.report()
    # m.fs.M1.report()
    # m.fs.product.report()
    # m.fs.tb_desal_disposal.report()
    m.fs.disposal.report()


if __name__ == "__main__":
    main()
