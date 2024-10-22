import os
import json
import pandas as pd
from collections import defaultdict
from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Expression,
    Block,
    Param,
    value,
    Var,
    NonNegativeReals,
    assert_optimal_termination,
    Objective,
    units as pyunits,
)
from pyomo.environ import TransformationFactory
from pyomo.network import Arc

import idaes.core.util.scaling as iscale
from idaes.core.util import DiagnosticsToolbox
from idaes.core.util.initialization import propagate_state
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.scaling import calculate_scaling_factors, set_scaling_factor, constraint_scaling_transform
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import Product, Feed
from idaes.models.unit_models.mixer import (
    Mixer,
    MomentumMixingType,
    MaterialBalanceType,
)
from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel

from watertap.costing import (
    WaterTAPCosting,
    PumpType,
    MixerType,
    ROType,
)
from watertap.unit_models.pressure_changer import Pump
from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock
from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
)
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *
from watertap.core.solvers import get_solver

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

atmospheric_pressure = 101325 * pyunits.Pa

class CCRO:

    """
    CCRO multiperiod model for validation and comparison against dynamic version.

    All input data must be in the units provided in data:
        volumetric flow: L/min
        concentration: g/L
            converted from conductivity with:
            (g/L) = 0.2674 * (mS/cm) ** 1.2404
        pressure: psi

    """

    def __init__(
        self,
        feed_conc=None,  # feed concentration is constant
        feed_flow=None,  # feed flow rate is constant
        reject_conc_start=None, # L/min
        reject_flow=None, # L/min
        water_recovery=None,
        n_time_points=3,
        temperature_start=None, # degC
        p1_pressure_start=None, # psi
        p2_pressure_start=None, # psi
        p1_eff=0.8,
        p2_eff=0.8,
        rho=1000,
        A_comp=4.2e-12,
        B_comp=3.5e-8,
        membrane_area=7.2,  # m2
        membrane_length=0.9626,  # m
        channel_height=1e-3,
        spacer_porosity=0.97,
        add_costing=False,
    ):
        self.rho = rho * pyunits.kg / pyunits.m**3

        self.feed_conc = feed_conc * pyunits.gram / pyunits.liter
        self.feed_flow = feed_flow * pyunits.liter / pyunits.minute
        self.reject_conc_start = reject_conc_start * pyunits.gram / pyunits.liter
        self.reject_flow = reject_flow * pyunits.liter / pyunits.minute
        self.water_recovery = water_recovery

        self.temperature_start = (temperature_start + 273.15) * pyunits.degK
        self.p1_pressure_start = p1_pressure_start * pyunits.psi

        self.n_time_points = n_time_points

        if p2_pressure_start is not None:
            self.p2_pressure_start = p2_pressure_start * pyunits.psi
        else:
            self.p2_pressure_start = self.p1_pressure_start

        self.p1_eff = p1_eff * pyunits.dimensionless
        self.p2_eff = p2_eff * pyunits.dimensionless

        self.feed_flow_mass_water = pyunits.convert(
            self.rho * self.feed_flow, to_units=pyunits.kg / pyunits.s
        )
        self.feed_flow_mass_salt = pyunits.convert(
            self.feed_conc * self.feed_flow, to_units=pyunits.kg / pyunits.s
        )
        self.reject_flow_mass_water = pyunits.convert(
            self.rho * self.reject_flow, to_units=pyunits.kg / pyunits.s
        )
        self.reject_flow_mass_salt = pyunits.convert(
            self.reject_conc_start * self.reject_flow,
            to_units=pyunits.kg / pyunits.s,
        )

        self.inlet_flow = self.feed_flow + self.reject_flow
        self.inlet_flow_mass_water = (
            self.feed_flow_mass_water + self.reject_flow_mass_water
        )
        self.inlet_conc_start = pyunits.convert(
            (self.feed_flow_mass_water + self.reject_flow_mass_water) / self.inlet_flow,
            to_units=pyunits.gram / pyunits.liter,
        )

        self.A_comp = A_comp
        self.B_comp = B_comp
        self.membrane_area = membrane_area * pyunits.m**2
        self.membrane_length = membrane_length * pyunits.m
        self.channel_height = channel_height * pyunits.m
        self.spacer_porosity = spacer_porosity * pyunits.dimensionless

        self.m = self.build_system()
        self.set_operating_conditions()
        # self.scale_system()
        # self.initialize()

    def check_jac(self, m, print_extreme_jacobian_values=True):
        jac, jac_scaled, nlp = iscale.constraint_autoscale_large_jac(m, min_scale=1e-8)
        try:
            cond_number = iscale.jacobian_cond(m, jac=jac_scaled) / 1e10
            print("--------------------------")
            print("COND NUMBER:", cond_number)
        except:
            print("Cond number failed")
            cond_number = None
        if print_extreme_jacobian_values:
            print("--------------------------")
            print("Extreme Jacobian entries:")
            extreme_entries = iscale.extreme_jacobian_entries(
                m, jac=jac_scaled, nlp=nlp, zero=1e-20, large=100
            )
            for val, var, con in extreme_entries:
                print(val, var.name, con.name)
            print("--------------------------")
            print("Extreme Jacobian columns:")
            extreme_cols = iscale.extreme_jacobian_columns(
                m, jac=jac_scaled, nlp=nlp, small=1e-3
            )
            for val, var in extreme_cols:
                print(val, var.name)
            print("------------------------")
            print("Extreme Jacobian rows:")
            extreme_rows = iscale.extreme_jacobian_rows(
                m, jac=jac_scaled, nlp=nlp, small=1e-3
            )
            for val, con in extreme_rows:
                print(val, con.name)
        return cond_number

    def build_system(self, time_blk=None):
        """
        Build steady-state model
        """
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = NaClParameterBlock()

        m.fs.feed = Feed(property_package=m.fs.properties)
        m.fs.product = Product(property_package=m.fs.properties)

        m.fs.P1 = Pump(property_package=m.fs.properties)
        m.fs.P2 = Pump(property_package=m.fs.properties)

        m.fs.M1 = Mixer(
            property_package=m.fs.properties,
            has_holdup=False,
            num_inlets=2,
            momentum_mixing_type=MomentumMixingType.none,
        )

        m.fs.RO = ReverseOsmosis0D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            pressure_change_type=PressureChangeType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            module_type="spiral_wound",
            has_full_reporting=True,
        )

        # connect unit models
        m.fs.feed_to_P1 = Arc(source=m.fs.feed.outlet, destination=m.fs.P1.inlet)
        m.fs.P1_to_M1 = Arc(source=m.fs.P1.outlet, destination=m.fs.M1.inlet_1)
        m.fs.P2_to_M1 = Arc(source=m.fs.P2.outlet, destination=m.fs.M1.inlet_2)
        m.fs.M1_to_RO = Arc(source=m.fs.M1.outlet, destination=m.fs.RO.inlet)

        m.fs.RO_permeate_to_product = Arc(
            source=m.fs.RO.permeate, destination=m.fs.product.inlet
        )
        m.fs.RO_retentate_to_P2 = Arc(
            source=m.fs.RO.retentate, destination=m.fs.P2.inlet
        )

        TransformationFactory("network.expand_arcs").apply_to(m)

        m.fs.water_recovery = Var(
            initialize=0.5,
            bounds=(0, 0.99),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="System Water Recovery",
        )

        m.fs.feed_salinity = Var(
            initialize=self.feed_conc,
            bounds=(0, 2000),
            domain=NonNegativeReals,
            units=pyunits.kg / pyunits.m**3,
            doc="Feed salinity",
        )

        m.fs.feed_flow_mass_water = Var(
            initialize=self.feed_flow_mass_water,
            bounds=(0.00001, 1e6),
            domain=NonNegativeReals,
            units=pyunits.kg / pyunits.s,
            doc="Mass flow water",
        )

        m.fs.feed_flow_vol_water = Var(
            initialize=self.feed_flow,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.liter / pyunits.min,
            doc="Feed tank, volumetric flow water",
        )

        m.fs.inlet_flow_vol_water = Expression(
            expr=pyunits.convert(
                m.fs.M1.mixed_state[0].flow_vol_phase["Liq"],
                to_units=pyunits.liter / pyunits.minute,
            )
        )

        m.fs.feed.properties[0].flow_vol_phase
        m.fs.feed.properties[0].conc_mass_phase_comp
        m.fs.M1.inlet_1_state[0].flow_vol_phase
        m.fs.M1.inlet_1_state[0].conc_mass_phase_comp
        m.fs.M1.inlet_2_state[0].flow_vol_phase
        m.fs.M1.inlet_2_state[0].conc_mass_phase_comp
        m.fs.M1.mixed_state[0].flow_vol_phase
        m.fs.M1.mixed_state[0].conc_mass_phase_comp

        self.scale_system(m=m)

        return m

    def scale_system(self, m=None):
        """
        Scale steady-state model.
        """
        if m is None:
            m = self.m

        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
        )

        set_scaling_factor(m.fs.P1.control_volume.work, 1e-2)
        set_scaling_factor(m.fs.P2.control_volume.work, 1e-2)
        set_scaling_factor(m.fs.RO.area, 1)

        set_scaling_factor(m.fs.water_recovery, 1e2)
        set_scaling_factor(m.fs.feed_flow_mass_water, 1e1)
        set_scaling_factor(m.fs.feed_salinity, 1)

        constraint_scaling_transform(m.fs.RO.feed_side.eq_K[0.0,0.0,'NaCl'], 1e7)
        constraint_scaling_transform(m.fs.RO.feed_side.eq_K[0.0,1.0,'NaCl'], 1e7)
        constraint_scaling_transform(m.fs.RO.eq_flux_mass[0.0,0.0,'Liq','NaCl'], 1e7)
        constraint_scaling_transform(m.fs.RO.eq_flux_mass[0.0,1.0,'Liq','NaCl'], 1e7)
        constraint_scaling_transform(m.fs.RO.eq_flux_mass[0.0,0.0,'Liq','H2O'], 1e4)
        constraint_scaling_transform(m.fs.RO.eq_flux_mass[0.0,1.0,'Liq','H2O'], 1e4)
        constraint_scaling_transform(m.fs.RO.eq_recovery_mass_phase_comp[0.0,'NaCl'], 1e4)
        constraint_scaling_transform(m.fs.RO.eq_mass_frac_permeate[0.0,0.0,'NaCl'], 1e5)
        constraint_scaling_transform(m.fs.RO.eq_mass_frac_permeate[0.0,1.0,'NaCl'], 1e5)
        constraint_scaling_transform(m.fs.RO.eq_permeate_production[0.0,'Liq','NaCl'], 1e4)

        calculate_scaling_factors(m)

    def set_operating_conditions(self, m=None):
        """
        Set operating conditions as initial conditions.
        """
        if m is None:
            m = self.m

        m.fs.feed_flow_mass_water.fix(self.feed_flow_mass_water)
        m.fs.feed_flow_vol_water.fix(self.feed_flow)
        m.fs.feed_salinity.fix(self.feed_conc)

        """
        Feed block operating conditions
        """

        m.fs.feed.properties[0].pressure.fix(atmospheric_pressure)
        m.fs.feed.properties[0].temperature.fix(self.temperature_start)

        m.fs.feed_flow_mass_water_constraint = Constraint(
            expr=m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
            # == m.fs.feed_flow_vol_water * self.rho
            == m.fs.feed_flow_mass_water
        )
        m.fs.feed_flow_salt_constraint = Constraint(
            expr=m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"] * self.rho
            == m.fs.feed_flow_mass_water * self.feed_conc
        )

        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].set_value(
            self.feed_flow_mass_water
        )
        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].set_value(
            self.feed_flow_mass_salt
        )

        """
        Pump 1 operating conditions
        """

        m.fs.P1.efficiency_pump.fix(self.p1_eff)
        m.fs.P1.control_volume.properties_out[0].pressure.fix(self.p1_pressure_start)
        
        """
        Pump 2 operating conditions
        """

        m.fs.P2.efficiency_pump.fix(self.p2_eff)
        m.fs.P2.control_volume.properties_out[0].pressure.fix(self.p1_pressure_start)
        m.fs.P2.control_volume.properties_out[0].pressure.setub(6e6)

        """
        Mixer operating conditions
        """

        # m.fs.M1_constraint_1 = Constraint(
        #     expr=m.fs.M1.inlet_2_state[0].pressure
        #     == m.fs.M1.inlet_1_state[0].pressure
        # )
        m.fs.M1_constraint_2 = Constraint(
            expr=m.fs.M1.inlet_1_state[0].pressure
            == m.fs.M1.mixed_state[0].pressure
        )

        # m.fs.M1.mixed_state[0].conc_mass_phase_comp["Liq", "NaCl"].setub(5)

        """
        RO operating conditions
        """

        m.fs.RO.permeate.pressure[0].fix(atmospheric_pressure)
        m.fs.RO.A_comp.fix(self.A_comp)
        m.fs.RO.B_comp.fix(self.B_comp)
        m.fs.RO.area.fix(self.membrane_area)
        m.fs.RO.length.fix(self.membrane_length)
        m.fs.RO.feed_side.channel_height.fix(self.channel_height)
        m.fs.RO.feed_side.spacer_porosity.fix(self.spacer_porosity)
        m.fs.RO.RO_constraint_1 = Constraint(expr=m.fs.RO.mixed_permeate[0].conc_mass_phase_comp["Liq", "NaCl"] <= 0.035) # Activating this reduces condition number by 8 OoM

        print("DOF =", degrees_of_freedom(m))
        print("DOF FEED =", degrees_of_freedom(m.fs.feed))
        print("DOF PUMP 1 =", degrees_of_freedom(m.fs.P1))
        print("DOF PUMP 2 =", degrees_of_freedom(m.fs.P2))
        print("DOF MIXER =", degrees_of_freedom(m.fs.M1))
        print("DOF RO =", degrees_of_freedom(m.fs.RO))
        assert_no_degrees_of_freedom(m)

    def initialize_mixer(self, m=None, guess=False):
        m.fs.M1.inlet_2_state[0].flow_mass_phase_comp["Liq", "H2O"].set_value(0.75465)
        m.fs.M1.inlet_2_state[0].flow_mass_phase_comp["Liq", "NaCl"].set_value(0.0029261)
        m.fs.M1.mixed_state[0].flow_mass_phase_comp["Liq", "H2O"].set_value(0.80332)
        m.fs.M1.mixed_state[0].flow_mass_phase_comp["Liq", "NaCl"].set_value(0.0030916)
        m.fs.M1.mixed_state[0].conc_mass_phase_comp["Liq", "NaCl"].set_value(4.0)
        m.fs.M1.inlet_2_state[0].pressure.set_value(2.1098e+06)
        m.fs.M1.initialize()

    def do_forward_initialization_pass(self, m=None, pass_num=1):
        """
        Initialize steady-state model
        """
        if m is None:
            m = self.m
        
        m.fs.feed.initialize()

        propagate_state(m.fs.feed_to_P1)
        m.fs.P1.initialize()

        propagate_state(m.fs.P1_to_M1)
       
        print('M1 first initialization')
        m.fs.M1.report()
        if pass_num > 0:
            m.fs.M1.initialize()
        else:
            self.initialize_mixer(m)
        m.fs.M1.report()

        propagate_state(m.fs.M1_to_RO)
        
        print('RO first initialization')
        m.fs.RO.feed_side.properties_in[0].pressure_osm_phase
        m.fs.RO.feed_side.properties_in[0].temperature = value(
            m.fs.feed.properties[0].temperature
        )
        m.fs.RO.feed_side.properties_in[0].pressure = value(
            m.fs.P1.control_volume.properties_out[0].pressure
        )
        m.fs.RO.feed_side.properties_out[0].flow_mass_phase_comp["Liq", "H2O"].set_value(0.75253)
        m.fs.RO.feed_side.properties_out[0].flow_mass_phase_comp["Liq", "NaCl"].set_value(0.0030899)
        m.fs.RO.feed_side.properties_out[0].pressure.set_value(1976757.5687263503)
        m.fs.RO.report()
        m.fs.RO.initialize()
        m.fs.RO.report()

        propagate_state(m.fs.RO_permeate_to_product)
        m.fs.product.initialize()

        propagate_state(m.fs.RO_retentate_to_P2)
        m.fs.P2.report()
        m.fs.P2.initialize()
        m.fs.P2.report()

        propagate_state(m.fs.P2_to_M1)
    
    def initialize(self, m=None):
        # initialize unit by unit
        dt = DiagnosticsToolbox(self.m)
        dt.report_structural_issues()
        
        # dt.display_constraints_with_extreme_jacobians()
        for idx, init_pass in enumerate(range(1)):
            print(f"\n\nINITIALIZATION PASS {idx+1}\n\n")
            self.do_forward_initialization_pass(m, pass_num=idx)
            self.print_results()
        self.check_jac(self.m)
        dt.report_numerical_issues()
        dt.display_constraints_with_large_residuals()
        # dt.display_variables_with_extreme_jacobians()
        # self.m.fs.P1.control_volume.properties_out[0].pressure.unfix()
        # self.m.fs.M1_constraint_1.deactivate()

    def fix_dof_and_initialize(self, m=None):
        """
        Fix DOF for MP model and initialize steady-state models.
        """
        if m is None:
            m = self.m
        self.set_operating_conditions(m=m)
        self.initialize_system(m=m)

    def solve(self, model=None, solver=None, tee=False, raise_on_failure=False):
        # ---solving---
        if solver is None:
            solver = get_solver()

        if model is None:
            model = self.m

        print("\n--------- SOLVING ---------\n")
        print(f"Degrees of Freedom: {degrees_of_freedom(model)}")
        optarg = {"tol": 1e-5, "linear_solver": "ma27", "max_iter": 1000}
        solver.options = optarg
        results = solver.solve(model, tee=tee)

        if check_optimal_termination(results):
            print("\n--------- OPTIMAL SOLVE!!! ---------\n")
            return results
        msg = "The current configuration is infeasible. Please adjust the decision variables."
        if raise_on_failure:
            print("\n--------- INFEASIBLE SOLVE!!! ---------\n")

            print("\n--------- CLOSE TO BOUNDS ---------\n")
            print_close_to_bounds(model)

            print("\n--------- INFEASIBLE BOUNDS ---------\n")
            print_infeasible_bounds(model)

            print("\n--------- INFEASIBLE CONSTRAINTS ---------\n")
            print_infeasible_constraints(model)

            raise RuntimeError(msg)
        else:
            print(msg)
            # debug(model, solver=solver, automate_rescale=False, resolve=False)
            # check_jac(model)
            # assert False
    
    def print_results(self, m=None):
        print("\n\n")
        print(
            f'MIXER INLET 1: {value(self.m.fs.M1.inlet_1_state[0].flow_mass_phase_comp["Liq", "H2O"]):<5.2f}'
        )
        print(
            f'MIXER INLET 2: {value(self.m.fs.M1.inlet_2_state[0].flow_mass_phase_comp["Liq", "H2O"]):<5.2f}'
        )
        print(
            f'MIXER OUTLET: {value(self.m.fs.M1.mixed_state[0].flow_mass_phase_comp["Liq", "H2O"]):<5.2f}'
        )
        print(
            f'MIXER CONC: {value(pyunits.convert(self.m.fs.M1.mixed_state[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.g/pyunits.L)):<5.2f} {pyunits.get_units(pyunits.convert(self.m.fs.M1.mixed_state[0].conc_mass_phase_comp["Liq", "NaCl"], to_units=pyunits.g/pyunits.L))}'
        )
        print("\n")
        print(
            f'PUMP 1 INLET: {value(self.m.fs.P1.control_volume.properties_in[0.0].flow_mass_phase_comp["Liq", "H2O"]):<5.2f}'
        )
        print(
            f'PUMP 1 OUTLET: {value(self.m.fs.P1.control_volume.properties_out[0.0].flow_mass_phase_comp["Liq", "H2O"]):<5.2f}'
        )
        print(
            f"PUMP 1 PRESSURE: {value(pyunits.convert(self.m.fs.P1.control_volume.properties_out[0.0].pressure, to_units=pyunits.bar)):<5.2f}"
        )
        print("\n")
        print(
            f'PUMP 2 INLET: {value(self.m.fs.P2.control_volume.properties_in[0.0].flow_mass_phase_comp["Liq", "H2O"]):<5.2f}'
        )
        print(
            f'PUMP 2 OUTLET: {value(self.m.fs.P2.control_volume.properties_out[0.0].flow_mass_phase_comp["Liq", "H2O"]):<5.2f}'
        )
        print(
            f"PUMP 2 PRESSURE: {value(pyunits.convert(self.m.fs.P2.control_volume.properties_out[0.0].pressure, to_units=pyunits.bar)):<5.2f}"
        )
        print("\n")
        print(f'RO FEED: {value(self.m.fs.RO.inlet.flow_mass_phase_comp[0,"Liq", "H2O"]):<5.2f}')
        print(
            f'RO PRODUCT: {value(self.m.fs.RO.permeate.flow_mass_phase_comp[0,"Liq", "H2O"]):<5.2f}'
        )
        print(
            f'RO BRINE: {value(self.m.fs.RO.retentate.flow_mass_phase_comp[0,"Liq", "H2O"]):<5.2f}'
        )
        print("\n\n")
        print(self.m.fs.M1.report())
        print(self.m.fs.RO.report())

if __name__ == "__main__":

    initial_conditions = {
        "feed_flow": 2.92,
        "feed_conc": 3.4,
        "reject_flow": 45.25,
        "reject_conc_start": 3.9,
        "temperature_start": 25,  # degC
        "p1_pressure_start": 306,  # psi
        "A_comp": 4.422e-12,
        "B_comp": 5.613e-8,
        "channel_height": 0.0008636,
        "spacer_porosity": 0.7081,
        "water_recovery": 0.063,
        "n_time_points": 5,
    }

    ccro = CCRO(**initial_conditions)
    ccro.initialize()
    ccro.solve(tee=False)
    dt = DiagnosticsToolbox(ccro.m)
    dt.compute_infeasibility_explanation()
    ccro.mp_df