from pyomo.environ import (ConcreteModel,
                           SolverFactory,
                           TerminationCondition,
                           value,
                           Constraint,
                           Expression,
                           Objective,
                           TransformationFactory,
                           units as pyunits)
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import (solve_indexed_blocks,
                                            propagate_state)
from idaes.generic_models.unit_models import Mixer, Separator, Product, Feed
from idaes.generic_models.unit_models.mixer import MomentumMixingType
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

# import proteuslib.property_models.seawater_prop_pack as props
import proteuslib.property_models.NaCl_prop_pack as props
from proteuslib.unit_models.reverse_osmosis_0D import ReverseOsmosis0D
from proteuslib.unit_models.pressure_exchanger import PressureExchanger
from proteuslib.unit_models.pump_isothermal import Pump
import proteuslib.flowsheets.RO_with_energy_recovery.financials as financials

def main():
    RO_system = ReverseOsmosisSystem()

    RO_system.build()

    RO_system.set_operating_conditions()

    RO_system.initialize_system()
    RO_system.solve(RO_system.m)

    print('\n***---Simulation results---****')
    RO_system.display_system(RO_system.m)
    RO_system.display_design(RO_system.m)
    RO_system.display_state(RO_system.m)

    RO_system.optimize()

    print('\n***---Optimization results---****')
    RO_system.display_system(RO_system.m)
    RO_system.display_design(RO_system.m)
    RO_system.display_state(RO_system.m)


class ReverseOsmosisSystem():
    """
    This class contains methods that build, set operating conditions, initialize,
    and optimize an RO system with energy recovery through a pressure exchanger
    """
    def __init__(self):
        # set up pyomo model
        self.m = ConcreteModel()
        # set up solver
        self.solver_str = 'ipopt'
        self.solver_opt = {'nlp_scaling_method': 'user-scaling'}
        self.solver = SolverFactory(self.solver_str)
        self.solver.options = self.solver_opt

    def build(self):
        # flowsheet set up
        self.m.fs = FlowsheetBlock(default={'dynamic': False})
        self.m.fs.properties = props.NaClParameterBlock()
        financials.add_costing_param_block(self.m.fs)

        # unit models
        self.m.fs.feed = Feed(default={'property_package': self.m.fs.properties})
        self.m.fs.S1 = Separator(default={
            "property_package": self.m.fs.properties,
            "outlet_list": ['P1', 'PXR']})
        self.m.fs.P1 = Pump(default={'property_package': self.m.fs.properties})
        self.m.fs.PXR = PressureExchanger(default={'property_package': self.m.fs.properties})
        self.m.fs.P2 = Pump(default={'property_package': self.m.fs.properties})
        self.m.fs.M1 = Mixer(default={
            "property_package": self.m.fs.properties,
            "momentum_mixing_type": MomentumMixingType.equality,  # booster pump will match pressure
            "inlet_list": ['P1', 'P2']})
        self.m.fs.RO = ReverseOsmosis0D(default={
            "property_package": self.m.fs.properties,
            "has_pressure_change": True})
        self.m.fs.product = Product(default={'property_package': self.m.fs.properties})
        self.m.fs.disposal = Product(default={'property_package': self.m.fs.properties})

        # additional variables or expressions
        feed_flow_vol_total = self.m.fs.feed.properties[0].flow_vol
        product_flow_vol_total = self.m.fs.product.properties[0].flow_vol
        self.m.fs.recovery = Expression(
            expr=product_flow_vol_total/feed_flow_vol_total)
        self.m.fs.annual_water_production = Expression(
            expr=pyunits.convert(product_flow_vol_total, to_units=pyunits.m ** 3 / pyunits.year)
                 * self.m.fs.costing_param.load_factor)
        pump_power_total = self.m.fs.P1.work_mechanical[0] + self.m.fs.P2.work_mechanical[0]
        self.m.fs.specific_energy_consumption = Expression(
            expr=pyunits.convert(pump_power_total, to_units=pyunits.kW)
                 / pyunits.convert(product_flow_vol_total, to_units=pyunits.m**3 / pyunits.hr))

        # costing
        self.m.fs.P1.get_costing(module=financials, pump_type="High pressure")
        self.m.fs.P2.get_costing(module=financials, pump_type="High pressure")
        self.m.fs.RO.get_costing(module=financials)
        self.m.fs.PXR.get_costing(module=financials)
        financials.get_system_costing(self.m.fs)

        # connections
        self.m.fs.s01 = Arc(source=self.m.fs.feed.outlet, destination=self.m.fs.S1.inlet)
        self.m.fs.s02 = Arc(source=self.m.fs.S1.P1, destination=self.m.fs.P1.inlet)
        self.m.fs.s03 = Arc(source=self.m.fs.P1.outlet, destination=self.m.fs.M1.P1)
        self.m.fs.s04 = Arc(source=self.m.fs.M1.outlet, destination=self.m.fs.RO.inlet)
        self.m.fs.s05 = Arc(source=self.m.fs.RO.permeate, destination=self.m.fs.product.inlet)
        self.m.fs.s06 = Arc(source=self.m.fs.RO.retentate, destination=self.m.fs.PXR.high_pressure_inlet)
        self.m.fs.s07 = Arc(source=self.m.fs.PXR.high_pressure_outlet, destination=self.m.fs.disposal.inlet)
        self.m.fs.s08 = Arc(source=self.m.fs.S1.PXR, destination=self.m.fs.PXR.low_pressure_inlet)
        self.m.fs.s09 = Arc(source=self.m.fs.PXR.low_pressure_outlet, destination=self.m.fs.P2.inlet)
        self.m.fs.s10 = Arc(source=self.m.fs.P2.outlet, destination=self.m.fs.M1.P2)
        TransformationFactory("network.expand_arcs").apply_to(self.m)

        # scaling
        self.m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
        self.m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'NaCl'))
        iscale.calculate_scaling_factors(self.m)

    def set_operating_conditions(self):
        # ---specifications---
        # feed
        feed_flow_mass = 1  # feed mass flow rate [kg/s]
        feed_mass_frac_NaCl = 0.035  # feed NaCl mass fraction [-]
        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl # feed H20 mass fraction [-]

        self.m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'NaCl'].fix(feed_flow_mass * feed_mass_frac_NaCl)
        self.m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(feed_flow_mass * feed_mass_frac_H2O)
        self.m.fs.feed.pressure.fix(101325)  # atmospheric pressure [Pa]
        self.m.fs.feed.temperature.fix(273.15 + 25)  # room temperature [K]

        # separator, no degrees of freedom (i.e. equal flow rates in PXR determines split fraction)

        # pump 1, high pressure pump, 2 degrees of freedom (efficiency and outlet pressure)
        self.m.fs.P1.efficiency_pump.fix(0.80)  # pump efficiency [-]
        operating_pressure = self.calculate_operating_pressure(
            over_pressure=0.3,
            water_recovery=0.5,
            NaCl_passage=0.01,
            feed_state_block=self.m.fs.feed.properties[0])
        self.m.fs.P1.control_volume.properties_out[0].pressure.fix(operating_pressure)

        # pressure exchanger
        self.m.fs.PXR.efficiency_pressure_exchanger.fix(0.95)  # pressure exchanger efficiency [-]

        # pump 2, booster pump, 1 degree of freedom (efficiency, pressure must match high pressure pump)
        self.m.fs.P2.efficiency_pump.fix(0.80)

        # mixer, no degrees of freedom

        # RO unit
        self.m.fs.RO.deltaP.fix(-3e5)  # pressure drop in membrane stage [Pa]
        self.m.fs.RO.A_comp.fix(4.2e-12)  # membrane water permeability coefficient [m/s-Pa]
        self.m.fs.RO.B_comp.fix(3.5e-8)  # membrane salt permeability coefficient [m/s]
        self.m.fs.RO.permeate.pressure[0].fix(101325)  # atmospheric pressure [Pa]
        # initiate RO feed values to determine area
        self.m.fs.RO.feed_side.properties_in[0].flow_mass_phase_comp['Liq', 'H2O'] = \
            value(self.m.fs.feed.properties[0].flow_mass_phase_comp['Liq', 'H2O'])
        self.m.fs.RO.feed_side.properties_in[0].flow_mass_phase_comp['Liq', 'NaCl'] = \
            value(self.m.fs.feed.properties[0].flow_mass_phase_comp['Liq', 'NaCl'])
        self.m.fs.RO.feed_side.properties_in[0].temperature = \
            value(self.m.fs.feed.properties[0].temperature)
        self.m.fs.RO.feed_side.properties_in[0].pressure = \
            value(self.m.fs.P1.control_volume.properties_out[0].pressure)
        RO_area = self.calculate_RO_area(unit=self.m.fs.RO, water_recovery=0.5)
        self.m.fs.RO.area.fix(RO_area)

        # check degrees of freedom
        if degrees_of_freedom(self.m) != 0:
            raise RuntimeError("The set_operating_conditions function resulted in {} "
                               "degrees of freedom rather than 0. This error suggests "
                               "that too many or not enough variables are fixed for a "
                               "simulation.".format(degrees_of_freedom(self.m)))

    def calculate_operating_pressure(self, feed_state_block=None,
            over_pressure=0.15, water_recovery=0.5, NaCl_passage=0.01):
        """
        estimate operating pressure for RO unit model given the following arguments:
            over_pressure:  the amount of operating pressure above the brine osmotic pressure
                            represented as a fraction (default=0.15)
            water_recovery: the mass-based fraction of inlet H2O that becomes permeate
                            (default=0.5)
            NaCl_passage:   the mass-based fraction of inlet NaCl that becomes permeate
                            (default=0.01)
            feed_state_block:    the state block of the RO feed that has the non-pressure state
                                 variables initialized to their values (default=None)
        """
        t = ConcreteModel()  # create temporary model
        prop = feed_state_block.config.parameters
        t.brine = prop.build_state_block([0], default={})

        # specify state block
        t.brine[0].flow_mass_phase_comp['Liq', 'H2O'].fix(
            value(feed_state_block.flow_mass_phase_comp['Liq', 'H2O']) * (1 - water_recovery))
        t.brine[0].flow_mass_phase_comp['Liq', 'NaCl'].fix(
            value(feed_state_block.flow_mass_phase_comp['Liq', 'NaCl']) * (1 - NaCl_passage))
        t.brine[0].pressure.fix(101325)  # valid when osmotic pressure is independent of hydraulic pressure
        t.brine[0].temperature.fix(value(feed_state_block.temperature))

        # calculate osmotic pressure
        # since properties are created on demand, we must touch the property to create it
        t.brine[0].pressure_osm
        # solve state block
        results = solve_indexed_blocks(self.solver, [t.brine])
        self.check_solve(results)

        return value(t.brine[0].pressure_osm) * (1 + over_pressure)

    def calculate_RO_area(self, unit=None, water_recovery=0.5):
        """
        determine RO membrane area required to achieve the specified water recovery:
            unit:  the RO unit model, e.g. m.fs.RO, it should have its inlet feed state block
                   initiated to the correct values (default=None)
            water_recovery: the mass-based fraction of inlet H2O that becomes permeate
                            (default=0.5)
        """
        # # fix inlet conditions
        unit.feed_side.properties_in[0].flow_mass_phase_comp.fix()
        unit.feed_side.properties_in[0].temperature.fix()
        unit.feed_side.properties_in[0].pressure.fix()
        # fix unit water recovery
        unit.feed_side.properties_out[0].flow_mass_phase_comp['Liq', 'H2O'].fix(
            unit.feed_side.properties_in[0].flow_mass_phase_comp['Liq', 'H2O'].value * (1 - water_recovery))
        # solve for unit area
        self.check_dof(unit)
        results = self.solver.solve(unit)
        self.check_solve(results)
        # unfix variables
        unit.feed_side.properties_in[0].flow_mass_phase_comp.unfix()
        unit.feed_side.properties_in[0].temperature.unfix()
        unit.feed_side.properties_in[0].pressure.unfix()
        unit.feed_side.properties_out[0].flow_mass_phase_comp['Liq', 'H2O'].unfix()
        return unit.area.value

    def initialize_system(self):

        # ---initialize feed block---
        self.m.fs.feed.initialize(solver=self.solver_str, optarg=self.solver_opt)

        # ---initialize splitter and pressure exchanger---
        # pressure exchanger high pressure inlet
        propagate_state(self.m.fs.s06)  # propagate to PXR high pressure inlet from RO retentate
        self.m.fs.PXR.high_pressure_side.properties_in.initialize(
            solver=self.solver_str, optarg=self.solver_opt)

        # splitter inlet
        propagate_state(self.m.fs.s01)  # propagate to splitter inlet from feed
        self.m.fs.S1.mixed_state[0].mass_frac_phase_comp  # touch property, so that it is built and can be solved for
        self.m.fs.S1.mixed_state.initialize(solver=self.solver_str, optarg=self.solver_opt)

        # splitter outlet to PXR, enforce same flow_vol as PXR high pressure inlet
        self.m.fs.S1.PXR_state[0].pressure.fix(value(self.m.fs.S1.mixed_state[0].pressure))
        self.m.fs.S1.PXR_state[0].temperature.fix(value(self.m.fs.S1.mixed_state[0].temperature))
        self.m.fs.S1.PXR_state[0].flow_vol_phase['Liq'].fix(
            value(self.m.fs.PXR.high_pressure_side.properties_in[0].flow_vol_phase['Liq']))
        self.m.fs.S1.PXR_state[0].mass_frac_phase_comp['Liq', 'NaCl'].fix(
            value(self.m.fs.S1.mixed_state[0].mass_frac_phase_comp['Liq', 'NaCl']))

        self.check_dof(self.m.fs.S1.PXR_state[0])
        results = solve_indexed_blocks(self.solver, [self.m.fs.S1.PXR_state])
        self.check_solve(results)

        # unfix PXR_state state variables and properties
        self.m.fs.S1.PXR_state[0].pressure.unfix()
        self.m.fs.S1.PXR_state[0].temperature.unfix()
        self.m.fs.S1.PXR_state[0].flow_vol_phase['Liq'].unfix()
        self.m.fs.S1.PXR_state[0].mass_frac_phase_comp['Liq', 'NaCl'].unfix()
        self.m.fs.S1.PXR_state[0].flow_mass_phase_comp['Liq', 'NaCl'].fix()

        # splitter initialization
        self.m.fs.S1.initialize(solver=self.solver_str, optarg=self.solver_opt)
        self.m.fs.S1.PXR_state[0].flow_mass_phase_comp['Liq', 'NaCl'].unfix()

        # pressure exchanger low pressure inlet
        propagate_state(self.m.fs.s08)

        # pressure exchanger initialization
        self.m.fs.PXR.initialize(solver=self.solver_str, optarg=self.solver_opt)

        # ---initialize pump 1---
        propagate_state(self.m.fs.s02)
        self.m.fs.P1.initialize(solver=self.solver_str, optarg=self.solver_opt)

        # ---initialize pump 2---
        propagate_state(self.m.fs.s09)
        self.m.fs.P2.control_volume.properties_out[0].pressure.fix(
            value(self.m.fs.P2.control_volume.properties_out[0].pressure))
        self.m.fs.P2.initialize(solver=self.solver_str, optarg=self.solver_opt)
        self.m.fs.P2.control_volume.properties_out[0].pressure.unfix()

        # ---initialize mixer---
        propagate_state(self.m.fs.s03)
        propagate_state(self.m.fs.s10)
        self.m.fs.M1.initialize(solver=self.solver_str, optarg=self.solver_opt, outlvl=idaeslog.INFO)


    def solve(self, blk, tee=False):
        results = self.solver.solve(blk, tee=tee)
        self.check_solve(results)

    @staticmethod
    def check_dof(blk, dof_expected=0):
        if degrees_of_freedom(blk) != dof_expected:
            raise RuntimeError("The degrees of freedom on {blk} were {dof} but {dof_e} "
                               "were expected, check the fixed variables on that block".format(
                blk=blk, dof=degrees_of_freedom(blk), dof_e=dof_expected))

    @staticmethod
    def check_solve(results):
        if results.solver.termination_condition != TerminationCondition.optimal:
            raise RuntimeError("The solver failed to converge to an optimal solution. "
                               "This suggests that the user provided infeasible inputs "
                               "or that the model is poorly scaled.")


    def optimize(self):
        # objective
        self.m.fs.objective = Objective(expr=self.m.fs.costing.LCOW)

        # unfix decision variables and add bounds
        # pump 1 and pump 2
        self.m.fs.P1.control_volume.properties_out[0].pressure.unfix()
        self.m.fs.P1.control_volume.properties_out[0].pressure.setlb(10e5)
        self.m.fs.P1.control_volume.properties_out[0].pressure.setub(80e5)
        self.m.fs.P1.deltaP.setlb(0)
        self.m.fs.P2.control_volume.properties_out[0].pressure.setlb(10e5)
        self.m.fs.P2.control_volume.properties_out[0].pressure.setub(80e5)
        self.m.fs.P2.deltaP.setlb(0)

        # RO
        self.m.fs.RO.area.unfix()  # area in membrane stage [m2]
        self.m.fs.RO.area.setlb(1)
        self.m.fs.RO.area.setub(100)

        # additional specifications
        product_recovery = 0.5  # product mass flow rate fraction of feed [-]
        product_salinity = 500e-6  # product NaCl mass fraction [-]

        # additional constraints
        self.m.fs.eq_recovery = Constraint(expr=product_recovery == self.m.fs.recovery)
        self.m.fs.eq_product_quality = Constraint(
            expr=self.m.fs.product.properties[0].mass_frac_phase_comp['Liq', 'NaCl'] <= product_salinity)
        iscale.constraint_scaling_transform(self.m.fs.eq_product_quality, 1e3)  # scaling constraint

        # ---checking model---
        self.check_dof(self.m, dof_expected=1)

        # --solve---
        self.solve(self.m)

    @staticmethod
    def display_system(m):
        print('---system metrics---')
        feed_flow_mass = sum(m.fs.feed.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
        feed_mass_frac_NaCl = m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / feed_flow_mass
        print('Feed: %.2f kg/s, %.0f ppm' % (feed_flow_mass, feed_mass_frac_NaCl * 1e6))

        prod_flow_mass = sum(m.fs.product.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
        prod_mass_frac_NaCl = m.fs.product.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / prod_flow_mass
        print('Product: %.3f kg/s, %.0f ppm' % (prod_flow_mass, prod_mass_frac_NaCl * 1e6))

        print('Recovery: %.1f%%' % (value(m.fs.recovery) * 100))
        print('Energy Consumption: %.1f kWh/m3' % value(m.fs.specific_energy_consumption))
        print('Levelized cost of water: %.2f $/m3' % value(m.fs.costing.LCOW))

    @staticmethod
    def display_design(m):
        print('---decision variables---')
        print('Operating pressure %.1f bar' % (m.fs.RO.inlet.pressure[0].value/1e5))
        print('Membrane area %.1f m2' % (m.fs.RO.area.value))

        print('---design variables---')
        print('Separator')
        print('Split fraction %.2f' % (m.fs.S1.split_fraction[0, 'PXR'].value*100))
        print('Pump 1 \noutlet pressure: %.1f bar \npower %.2f kW'
              % (m.fs.P1.outlet.pressure[0].value / 1e5, m.fs.P1.work_mechanical[0].value / 1e3))
        print('Pump 2 \noutlet pressure: %.1f bar \npower %.2f kW'
              % (m.fs.P2.outlet.pressure[0].value / 1e5, m.fs.P2.work_mechanical[0].value / 1e3))

    @staticmethod
    def display_state(m):
        print('---state---')
        def print_state(s, b):
            flow_mass = sum(b.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
            mass_frac_ppm = b.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / flow_mass * 1e6
            pressure_bar = b.pressure[0].value / 1e5
            print(s + ': %.3f kg/s \t%.0f ppm \t%.1f bar' % (flow_mass, mass_frac_ppm, pressure_bar))

        print_state('Feed      ', m.fs.feed.outlet)
        print_state('Split 1   ', m.fs.S1.P1)
        print_state('P1 out    ', m.fs.P1.outlet)
        print_state('Split 2   ', m.fs.S1.PXR)
        print_state('PXR LP out', m.fs.PXR.low_pressure_outlet)
        print_state('P2 out    ', m.fs.P2.outlet)
        print_state('Mix out   ', m.fs.M1.outlet)
        print_state('RO perm   ', m.fs.RO.permeate)
        print_state('RO reten  ', m.fs.RO.retentate)
        print_state('PXR HP out', m.fs.PXR.high_pressure_outlet)


if __name__ == "__main__":
    main()
