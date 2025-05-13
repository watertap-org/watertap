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

from pyomo.util.calc_var_value import calculate_variable_from_constraint

import idaes.core.util.scaling as iscale
from idaes.core.util.initialization import propagate_state
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.scaling import calculate_scaling_factors, set_scaling_factor
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
    # ReverseOsmosis1D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.unit_models.reverse_osmosis_0D import (
    # ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
)
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.core.util.initialization import *
from watertap.core.solvers import get_solver
from watertap.flowsheets.ccro.multiperiod.model_state_tool import ModelState


from watertap.unit_models.pseudo_steady_state import CCRO1D, DeadVolume0D

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

atmospheric_pressure = 101325 * pyunits.Pa

convert_dict = json.load(open(f"{__location__}/column_convert.json"))
mp_results_cols = json.load(open(f"{__location__}/mp_results_cols.json"))


def get_variable_pairs(t1, t2):

    return [
        (
            t2.fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"],
            t1.fs.dead_volume.dead_volume.properties_out[0].mass_frac_phase_comp[
                "Liq", "NaCl"
            ],
        ),
        (
            t2.fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"],
            t1.fs.dead_volume.dead_volume.properties_out[0].dens_mass_phase["Liq"],
        ),
    ]


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
        reject_conc_start=None,  # L/min
        reject_flow=None,  # L/min
        water_recovery=None,
        n_time_points=3,
        temperature_start=None,  # degC
        p1_pressure_start=None,  # psi
        p2_pressure_start=None,  # psi
        p1_eff=0.8,
        p2_eff=0.8,
        rho=1000,
        A_comp=4.2e-12,
        B_comp=3.5e-8,
        membrane_area=7.2,  # m2
        membrane_length=0.9626,  # m
        channel_height=1e-3,
        spacer_porosity=0.97,
        include_costing=False,
        dead_volume=0.1,
        accumulation_time=300,
    ):
        self.rho = rho * pyunits.kg / pyunits.m**3
        self.dead_volume = dead_volume * pyunits.m**3

        self.feed_conc = feed_conc * pyunits.gram / pyunits.liter
        self.feed_flow = feed_flow * pyunits.liter / pyunits.minute
        self.reject_conc_start = reject_conc_start * pyunits.gram / pyunits.liter
        self.reject_flow = reject_flow * pyunits.liter / pyunits.minute
        self.water_recovery = water_recovery

        self.temperature_start = (temperature_start + 273.15) * pyunits.degK
        self.p1_pressure_start = p1_pressure_start * pyunits.psi

        self.n_time_points = n_time_points

        self.accumulation_time = accumulation_time * pyunits.s
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
        self.feed_as_reject_flow_mass_water = pyunits.convert(
            self.rho * self.reject_flow, to_units=pyunits.kg / pyunits.s
        )
        self.feed_as_reject_flow_mass_salt = pyunits.convert(
            self.feed_conc * self.reject_flow,
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

        self.include_costing = include_costing

        watertap_solver = get_solver()

        self.mp = MultiPeriodModel(
            n_time_points=self.n_time_points,
            process_model_func=self.build_system,
            linking_variable_func=get_variable_pairs,
            initialization_func=self.fix_dof_and_initialize,
            unfix_dof_func=self.unfix_dof,
            solver=watertap_solver,
            outlvl=logging.WARNING,
        )

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
            momentum_mixing_type=MomentumMixingType.equality,
        )

        m.fs.RO = CCRO1D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            pressure_change_type=PressureChangeType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            transformation_scheme="BACKWARD",
            transformation_method="dae.finite_difference",
            finite_elements=10,
            module_type="spiral_wound",
            has_full_reporting=True,
            # feed_pump=m.fs.P1,
        )

        m.fs.dead_volume = DeadVolume0D(property_package=m.fs.properties)

        # connect unit models
        m.fs.feed_to_P1 = Arc(source=m.fs.feed.outlet, destination=m.fs.P1.inlet)
        m.fs.P1_to_M1 = Arc(source=m.fs.P1.outlet, destination=m.fs.M1.inlet_1)
        m.fs.P2_to_M1 = Arc(source=m.fs.P2.outlet, destination=m.fs.M1.inlet_2)
        m.fs.M1_to_RO = Arc(source=m.fs.M1.outlet, destination=m.fs.RO.inlet)

        m.fs.RO_permeate_to_product = Arc(
            source=m.fs.RO.permeate, destination=m.fs.product.inlet
        )
        m.fs.RO_retentate_to_dead_volume = Arc(
            source=m.fs.RO.retentate, destination=m.fs.dead_volume.inlet
        )

        m.fs.dead_volume_to_P2 = Arc(
            source=m.fs.dead_volume.outlet, destination=m.fs.P2.inlet
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

        dv = m.fs.dead_volume.dead_volume
        dv.properties_in[0].flow_vol_phase["Liq"]
        dv.properties_in[0].conc_mass_phase_comp["Liq", "NaCl"]
        dv.properties_out[0].flow_vol_phase["Liq"]
        dv.properties_out[0].conc_mass_phase_comp["Liq", "NaCl"]

        ### MUST BUILD THIS ON THE MODEL BEFORE PASSING INTO MULTIPERIOD MODEL BLOCK
        m.fs.feed_flow_mass_water_constraint = Constraint(
            expr=m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]
            == m.fs.feed_flow_mass_water
        )

        m.fs.feed.properties[0].dens_mass_phase["Liq"].set_value(self.rho)
        m.fs.feed_flow_salt_constraint = Constraint(
            expr=m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"]
            * m.fs.feed.properties[0].dens_mass_phase["Liq"]
            == m.fs.feed_flow_mass_water * self.feed_conc
        )
        return m

    def add_costing(self, m=None):
        """
        Add costing blocks to steady-state model.
        """
        if m is None:
            m = self.m
        m.fs.costing = WaterTAPCosting()
        m.fs.RO.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        # m.fs.P1.costing = UnitModelCostingBlock(
        #     flowsheet_costing_block=m.fs.costing,
        # )
        # m.fs.P2.costing = UnitModelCostingBlock(
        #     flowsheet_costing_block=m.fs.costing,
        # )

        m.fs.costing.cost_process()

        m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol_phase["Liq"])
        m.fs.costing.add_specific_energy_consumption(
            m.fs.product.properties[0].flow_vol_phase["Liq"], name="SEC"
        )

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

        set_scaling_factor(m.fs.P1.control_volume.work, 1e-3)
        set_scaling_factor(m.fs.P2.control_volume.work, 1e-3)
        # set_scaling_factor(m.fs.P1.control_volume.work, 1e-5)
        # set_scaling_factor(m.fs.P2.control_volume.work, 1e-5)
        set_scaling_factor(m.fs.RO.area, 1e-2)

        set_scaling_factor(m.fs.water_recovery, 10)
        set_scaling_factor(m.fs.feed_flow_mass_water, 1)
        set_scaling_factor(m.fs.feed_salinity, 1)

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
        # TODO:  These constraints need to be scaled!
        calculate_variable_from_constraint(
            m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"],
            m.fs.feed_flow_mass_water_constraint,
        )
        calculate_variable_from_constraint(
            m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"],
            m.fs.feed_flow_salt_constraint,
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
        # m.fs.P2.control_volume.properties_out[0].pressure.set_value(self.p2_pressure_start)

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

        m.fs.RO.feed_side.K.setlb(1e-6)
        m.fs.RO.feed_side.friction_factor_darcy.setub(200)
        m.fs.RO.flux_mass_phase_comp.setub(1)
        m.fs.RO.feed_side.cp_modulus.setub(50)
        m.fs.RO.feed_side.cp_modulus.setlb(0.1)
        m.fs.RO.deltaP.setlb(None)

        for e in m.fs.RO.permeate_side:
            if e[-1] != 0:
                m.fs.RO.permeate_side[e].pressure_osm_phase["Liq"].setlb(200)
                m.fs.RO.permeate_side[e].molality_phase_comp["Liq", "NaCl"].setlb(1e-8)
        # assume there is no change in volume
        m.fs.dead_volume.volume.fix(self.dead_volume)
        m.fs.dead_volume.delta_state.volume.fix(self.dead_volume)
        m.fs.dead_volume.accumulation_time.fix(self.accumulation_time)

        # initialize state to assumed initial condition, use concentration of feed
        # and density of feed
        # Solve feed block first

        m.fs.feed.properties[0].pressure_osm_phase["Liq"]
        m.fs.feed.properties[0].flow_vol_phase["Liq"].fix(self.reject_flow)
        m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].fix(
            self.reject_conc_start
        )
        solver = get_solver()
        solver.solve(m.fs.feed)
        m.fs.feed.properties[0].flow_vol_phase["Liq"].unfix()
        m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].unfix()

        # I found fixing mass fraction and density is easiest way to get initial state
        # we will also use these as connection points between current and future state.

        m.fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"].fix(
            m.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"].value
        )
        m.fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"].fix(
            m.fs.feed.properties[0].dens_mass_phase["Liq"].value
        )

        # initialize feed to desired flow and conc
        m.fs.feed.properties[0].flow_vol_phase["Liq"].fix(self.feed_flow)
        m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].fix(self.feed_conc)
        solver = get_solver()
        solver.solve(m.fs.feed)
        m.fs.feed.properties[0].flow_vol_phase["Liq"].unfix()
        m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].unfix()
        print("DOF =", degrees_of_freedom(m))
        print("DOF FEED =", degrees_of_freedom(m.fs.feed))
        print("DOF PUMP 1 =", degrees_of_freedom(m.fs.P1))
        print("DOF PUMP 2 =", degrees_of_freedom(m.fs.P2))
        print("DOF MIXER =", degrees_of_freedom(m.fs.M1))
        print("DOF RO =", degrees_of_freedom(m.fs.RO))
        print("DOF Dead Volume =", degrees_of_freedom(m.fs.dead_volume))
        assert_no_degrees_of_freedom(m)
        self.scale_system(m)

    def copy_inlet_state_for_mixer(self, m):
        for idx, obj in m.fs.M1.inlet_2.flow_mass_phase_comp.items():
            obj.value = m.fs.M1.inlet_1.flow_mass_phase_comp[idx].value * 1

    def initialize_system(self, m=None):
        """
        Initialize steady-state model
        """
        if m is None:
            m = self.m

        # feed is initialized when we setup our fixed operating conditions

        propagate_state(m.fs.feed_to_P1)
        m.fs.P1.outlet.pressure[0].fix(
            m.fs.feed.properties[0].pressure_osm_phase["Liq"].value * 5
        )
        m.fs.P2.outlet.pressure[0].fix(
            m.fs.feed.properties[0].pressure_osm_phase["Liq"].value * 5
        )

        m.fs.P1.initialize()

        propagate_state(m.fs.P1_to_M1)
        self.copy_inlet_state_for_mixer(m)
        m.fs.M1.initialize()
        propagate_state(m.fs.M1_to_RO)

        m.fs.RO.initialize()

        propagate_state(m.fs.RO_permeate_to_product)
        propagate_state(m.fs.RO_retentate_to_dead_volume)
        m.fs.dead_volume.initialize()

        propagate_state(m.fs.dead_volume_to_P2)
        m.fs.P2.initialize()
        m.fs.P2.outlet.pressure[0].unfix()
        propagate_state(m.fs.P2_to_M1)
        m.fs.product.initialize()

    def create_multiperiod(self):
        """
        Create MultiPeriod model
        """

        # watertap_solver = get_solver()

        # self.mp = mp = MultiPeriodModel(
        #     n_time_points=self.n_time_points,
        #     process_model_func=self.build_system,
        #     linking_variable_func=get_variable_pairs,
        #     initialization_func=self.fix_dof_and_initialize,
        #     unfix_dof_func=self.unfix_dof,
        #     solver=watertap_solver,
        #     outlvl=logging.WARNING,
        # )

        self.flowsheet_options = {
            t: {
                "time_blk": t,
            }
            for t in range(self.n_time_points)
        }

        mp = self.mp

        mp.build_multi_period_model(
            model_data_kwargs=self.flowsheet_options,
            # unfix_dof_options={"water_recovery": self.water_recovery, "Q_ro": self.inlet_flow_mass_water},
        )

        tot_perm_vol = 0
        tot_brine_vol = 0
        for t, m in enumerate(mp.get_active_process_blocks()):
            tot_perm_vol += m.fs.RO.mixed_permeate[0].flow_vol_phase["Liq"]
            tot_brine_vol += m.fs.RO.feed_side.properties[0, 1.0].flow_vol_phase["Liq"]
            if t == 0:
                self.m0 = m
                self.fix_dof_and_initialize(m=m)
                self.unfix_dof(
                    m=m, time_idx=t
                )  # ensure we do not unfix dead volume stuff
                old_m = m
                self.add_costing(m=m)
            else:
                self.copy_state_prop_time_period_links(old_m, m)

            results = self.solve(model=m, tee=True)
            assert_optimal_termination(results)
            self.unfix_dof(m=m, time_idx=t)
            old_m = m

        mp.recovery = Expression(
            expr=tot_perm_vol / (tot_perm_vol + tot_brine_vol),
        )

    def copy_state_prop_time_period_links(self, m_old, m_new):
        self.copy_state(m_old, m_new)
        m_new.fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"].fix(
            m_old.fs.dead_volume.dead_volume.properties_out[0]
            .mass_frac_phase_comp["Liq", "NaCl"]
            .value
        )
        m_new.fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"].fix(
            m_old.fs.dead_volume.dead_volume.properties_out[0]
            .dens_mass_phase["Liq"]
            .value
        )

    def copy_state(self, old_model, new_model):
        model_state = ModelState()
        model_state.get_model_state(old_model)
        model_state.set_model_state(new_model)

    def fix_dof_and_initialize(self, m=None):
        """
        Fix DOF for MP model and initialize steady-state models.
        """
        if m is None:
            m = self.m

        self.set_operating_conditions(m=m)
        self.initialize_system(m=m)

    def unfix_dof(self, time_idx=0, m=None):
        """
        Unfix linking variables in MP model
        """
        if m is None:
            m = self.m
        m.fs.RO.recovery_vol_phase[0, "Liq"].fix(self.water_recovery)
        m.fs.P1.control_volume.properties_out[0].pressure.unfix()
        m.fs.P2.control_volume.properties_out[0].pressure.unfix()

        if time_idx > 0:
            m.fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"].unfix()
            m.fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"].unfix()

    def solve(self, model=None, solver=None, tee=False, raise_on_failure=True):
        # ---solving---
        if solver is None:
            solver = get_solver()

        if model is None:
            model = self.mp

        print("\n--------- SOLVING ---------\n")
        print(f"Degrees of Freedom: {degrees_of_freedom(model)}")
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
            assert False

    def extract_multiperiod_data(self, mp=None):
        """
        Extract results from MP model as time series.

        class attribute mp_df will have index corresponding to each period
        with columns found in mp_results_col (loaded as .json at top of file)
        and select columns converted to units found in dataset.

        """
        if mp is None:
            mp = self.mp

        unit_model_blocks = ["feed", "P1", "M1", "P2", "RO", "dead_volume", "product"]
        self.results_dict = results_dict = dict()

        for t, m in enumerate(mp.get_active_process_blocks()):
            results_dict[t] = dict()

            for fs in m.component_objects(
                Block,
                descend_into=False,
            ):
                for b in fs.component_objects(Block, descend_into=False):
                    bname = b.name.replace(f"blocks[{t}].process.fs.", "").split(".")[0]
                    if bname not in unit_model_blocks:
                        continue
                    results_dict[t][bname] = defaultdict(float)

                    for x in b.component_objects(
                        [Var, Expression],
                        descend_into=True,
                    ):
                        if (
                            x.name.replace(f"blocks[{t}].process.fs.{bname}.", "")[0]
                            == "_"
                        ):  # skip property _refs
                            continue
                        if x.is_indexed():
                            for i, xx in x.items():
                                xname = xx.name.replace(
                                    f"blocks[{t}].process.fs.{bname}.", ""
                                ).replace("0.0", "0")
                                results_dict[t][bname][xname] = value(x[i])
                        else:
                            xname = x.name.replace(
                                f"blocks[{t}].process.fs.{bname}.", ""
                            ).replace("0.0", "0")
                            results_dict[t][bname][xname] = value(x)

        self.mp_results_dict = mp_results_dict = defaultdict(list)

        for t, d in results_dict.items():
            for b, p in d.items():
                for k, v in p.items():
                    if "costing" in k:
                        continue
                    if b in mp_results_cols.keys():
                        # if k in mp_results_cols[b]:
                        mp_results_dict[f"{b}.{k}"].append(v)

        self.mp_df = mp_df = pd.DataFrame.from_dict(mp_results_dict)

        for c in mp_df.columns:
            if c in convert_dict.keys():
                from_u = convert_dict[c]["units"][0]
                to_u = convert_dict[c]["units"][1]
                if "/" in from_u:
                    fnum, fden = [*from_u.split("/")]
                    from_units = getattr(pyunits, fnum) / getattr(pyunits, fden)
                    tnum, tden = [*to_u.split("/")]
                    to_units = getattr(pyunits, tnum) / getattr(pyunits, tden)
                else:
                    from_units = getattr(pyunits, from_u)
                    to_units = getattr(pyunits, to_u)

                mp_df[convert_dict[c]["col_name"]] = mp_df[c].apply(
                    lambda x: value(pyunits.convert(x * from_units, to_units=to_units))
                )
        for i in mp_df.index:
            mp_df.at[i, "time"] = i * value(self.accumulation_time)


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
        "n_time_points": 3,
    }

    ccro = CCRO(**initial_conditions)

    ccro.create_multiperiod()
    ccro.solve(tee=True)
    ccro.extract_multiperiod_data()
    ccro.mp_df
