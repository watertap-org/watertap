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
from idaes.models.unit_models import Product, Feed, StateJunction
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
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
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
from idaes.core.scaling import  report_scaling_factors, AutoScaler
import numpy as np
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

atmospheric_pressure = 101325 * pyunits.Pa

   
def build_system():
    """
    Build steady-state model
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True,
                          time_set=list(np.linspace(0, 200, 6)),
                          time_units=pyunits.s)
    m.fs.properties = NaClParameterBlock()

    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)

    m.fs.P1 = Pump(property_package=m.fs.properties)
    m.fs.P2 = Pump(property_package=m.fs.properties)

    m.fs.M1 = Mixer(
        property_package=m.fs.properties,
        # has_holdup=False,
        # num_inlets=2,
        momentum_mixing_type=MomentumMixingType.equality,
    )
    m.fs.M1.pressure_equality_constraints[0,2].deactivate()

    m.fs.RO = ReverseOsmosis0D(
        dynamic=True,
        has_holdup=True,
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        module_type="spiral_wound",
        has_full_reporting=True,
    )
    time_nfe = len(m.fs.time) - 1
    TransformationFactory("dae.finite_difference").apply_to(
        m.fs, nfe=time_nfe, wrt=m.fs.time, scheme="BACKWARD"
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

    # m.fs.water_recovery = Var(
    #     initialize=0.5,
    #     bounds=(0, 0.99),
    #     domain=NonNegativeReals,
    #     units=pyunits.dimensionless,
    #     doc="System Water Recovery",
    # )

    # m.fs.feed_salinity = Var(
    #     initialize=self.feed_conc,
    #     bounds=(0, 2000),
    #     domain=NonNegativeReals,
    #     units=pyunits.kg / pyunits.m**3,
    #     doc="Feed salinity",
    # )

    # m.fs.feed_flow_mass_water = Var(
    #     initialize=self.feed_flow_mass_water,
    #     bounds=(0.00001, 1e6),
    #     domain=NonNegativeReals,
    #     units=pyunits.kg / pyunits.s,
    #     doc="Mass flow water",
    # )

    # m.fs.feed_flow_vol_water = Var(
    #     initialize=self.feed_flow,
    #     bounds=(0, None),
    #     domain=NonNegativeReals,
    #     units=pyunits.liter / pyunits.min,
    #     doc="Feed tank, volumetric flow water",
    # )

    # m.fs.inlet_flow_vol_water = Expression(
    #     expr=pyunits.convert(
    #         m.fs.M1.mixed_state[0].flow_vol_phase["Liq"],
    #         to_units=pyunits.liter / pyunits.minute,
    #     )
    # )

    # m.fs.feed.properties[0].flow_vol_phase
    # m.fs.feed.properties[0].conc_mass_phase_comp
    # m.fs.M1.inlet_1_state[0].flow_vol_phase
    # m.fs.M1.inlet_1_state[0].conc_mass_phase_comp
    # m.fs.M1.inlet_2_state[0].flow_vol_phase
    # m.fs.M1.inlet_2_state[0].conc_mass_phase_comp
    # m.fs.M1.mixed_state[0].flow_vol_phase
    # m.fs.M1.mixed_state[0].conc_mass_phase_comp


    return m

def scale_system(m):
    """
    Scale steady-state model.
    """

    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e4, index=("Liq", "NaCl")
    )

    # set_scaling_factor(m.fs.RO.permeate_side[0,0].flow_mass_phase_comp["Liq", "NaCl"], 1e6)
    # set_scaling_factor(m.fs.RO.permeate_side[0,1].flow_mass_phase_comp["Liq", "NaCl"], 1e4)
    # set_scaling_factor(m.fs.RO.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"], 1e6)
    # set_scaling_factor(m.fs.P1.control_volume.work, 1e-2)
    # set_scaling_factor(m.fs.P2.control_volume.work, 1e-2)
    # set_scaling_factor(m.fs.RO.area, 1)

    # set_scaling_factor(m.fs.water_recovery, 1e2)
    # set_scaling_factor(m.fs.feed_flow_mass_water, 1e1)
    # set_scaling_factor(m.fs.feed_salinity, 1)

    # constraint_scaling_transform(m.fs.RO.feed_side.eq_K[0.0,0.0,'NaCl'], 1e7)
    # constraint_scaling_transform(m.fs.RO.feed_side.eq_K[0.0,1.0,'NaCl'], 1e7)
    # constraint_scaling_transform(m.fs.RO.eq_flux_mass[0.0,0.0,'Liq','NaCl'], 1e7)
    # constraint_scaling_transform(m.fs.RO.eq_flux_mass[0.0,1.0,'Liq','NaCl'], 1e7)
    # constraint_scaling_transform(m.fs.RO.eq_flux_mass[0.0,0.0,'Liq','H2O'], 1e4)
    # constraint_scaling_transform(m.fs.RO.eq_flux_mass[0.0,1.0,'Liq','H2O'], 1e4)
    # constraint_scaling_transform(m.fs.RO.eq_recovery_mass_phase_comp[0.0,'NaCl'], 1e4)
    # constraint_scaling_transform(m.fs.RO.eq_mass_frac_permeate[0.0,0.0,'NaCl'], 1e5)
    # constraint_scaling_transform(m.fs.RO.eq_mass_frac_permeate[0.0,1.0,'NaCl'], 1e5)
    # constraint_scaling_transform(m.fs.RO.eq_permeate_production[0.0,'Liq','NaCl'], 1e4)

    calculate_scaling_factors(m)

def set_operating_conditions(m):
    """
    Set operating conditions as initial conditions.
    """

    # m.fs.feed_flow_mass_water.fix(self.feed_flow_mass_water)
    # m.fs.feed_flow_vol_water.fix(self.feed_flow)
    # m.fs.feed_salinity.fix(self.feed_conc)

    """
    Feed block operating conditions
    """
    feed_temp = 25 +273.15
    conc_mass_tds = 3.4 * pyunits.kg/pyunits.m**3       
    flow_vol = 2.92 * pyunits.L/pyunits.min



    m.fs.feed.properties[:].pressure.fix(atmospheric_pressure)
    m.fs.feed.properties[:].temperature.fix(feed_temp)

    m.fs.feed.properties.calculate_state(
        var_args={
            ("conc_mass_phase_comp", ("Liq", "NaCl")): conc_mass_tds,
            ("flow_vol_phase", "Liq"): flow_vol,
        },
        hold_state=True,
    )
    

    """
    Pump 1 operating conditions
    """
    p1_eff = 0.8    
    p1_pressure_start = 306 * pyunits.psi#306 * pyunits.psi

    m.fs.P1.efficiency_pump.fix(p1_eff)
    m.fs.P1.control_volume.properties_out[0].pressure.fix(p1_pressure_start)
    
    """
    Pump 2 operating conditions
    """
    p2_eff = 0.8
    m.fs.P2.efficiency_pump.fix(p2_eff)
    m.fs.P2.control_volume.properties_out[0].pressure.fix(p1_pressure_start)
    # m.fs.P2.control_volume.properties_out[0].pressure.setub(6e6)

    """
    Mixer operating conditions
    """

    # m.fs.M1_constraint_1 = Constraint(
    #     expr=m.fs.M1.inlet_2_state[0].pressure
    #     == m.fs.M1.inlet_1_state[0].pressure
    # )
    # m.fs.M1_constraint_2 = Constraint(
    #     expr=m.fs.M1.inlet_1_state[0].pressure
    #     == m.fs.M1.mixed_state[0].pressure
    # )

    # m.fs.M1.mixed_state[0].conc_mass_phase_comp["Liq", "NaCl"].setub(5)

    """
    RO operating conditions
    """

    A_comp= 4.422e-12
    B_comp=5.613e-8
    membrane_area=7.2 
    membrane_length=0.9626
    channel_height=0.0008636
    spacer_porosity=0.7081
    m.fs.RO.permeate.pressure[0].fix(atmospheric_pressure)
    m.fs.RO.A_comp.fix(A_comp)
    m.fs.RO.B_comp.fix(B_comp)
    set_scaling_factor(m.fs.RO.B_comp[ 0, "NaCl"], 10/value(m.fs.RO.B_comp[ 0, "NaCl"]))
    m.fs.RO.area.fix(membrane_area)
    m.fs.RO.length.fix(membrane_length)
    m.fs.RO.feed_side.channel_height.fix(channel_height)
    m.fs.RO.feed_side.spacer_porosity.fix(spacer_porosity)
  
    # m.fs.RO.recovery_vol_phase[0, "Liq"].fix(0.06)
    m.fs.unit.feed_side.material_accumulation[:, :, :].value = 0
    m.fs.unit.feed_side.material_accumulation[0, :, :].fix(0)
    
    m.fs.RO.feed_side.K.setlb(1e-9)
    m.fs.RO.feed_side.friction_factor_darcy.setub(None)
    m.fs.RO.permeate_side[...].conc_mass_phase_comp.setlb(None)
    # m.fs.RO.RO_constraint_1 = Constraint(expr=m.fs.RO.mixed_permeate[0].conc_mass_phase_comp["Liq", "NaCl"] <= 0.5 * pyunits.kg/pyunits.m**3) # Activating this reduces condition number by 8 OoM
    # m.fs.RO.deltaP.fix(-1e5)
    print("DOF =", degrees_of_freedom(m))
    print("DOF FEED =", degrees_of_freedom(m.fs.feed))
    print("DOF PUMP 1 =", degrees_of_freedom(m.fs.P1))
    print("DOF PUMP 2 =", degrees_of_freedom(m.fs.P2))
    print("DOF MIXER =", degrees_of_freedom(m.fs.M1))
    print("DOF RO =", degrees_of_freedom(m.fs.RO))
    assert_no_degrees_of_freedom(m)

def initialize_system(m):
    
    m.fs.feed.initialize()

    propagate_state(m.fs.feed_to_P1)
    m.fs.P1.initialize()

    propagate_state(m.fs.P1_to_M1)



    master_initialize_with_recirculation(m, count=3)
    
def master_initialize_with_recirculation(m, count=1):
    solved = 0 
    counter = 0
    while not solved:
        try:
            initialize_with_recirculation(m)
            res= solve(m, tee=True)
        except:
            pass
        counter = counter + 1
        if counter == count:
            break 
        if check_optimal_termination(res):
            solved = 1
        m.fs.report()
        m.fs.RO.report()

def initialize_with_recirculation(m):
    propagate_state(source=m.fs.P1.outlet, destination=m.fs.P2.outlet)
    propagate_state(m.fs.P2_to_M1)
    m.fs.M1.initialize()
    m.fs.M1.pressure_equality_constraints[0,2].deactivate()

    # mixer to RO 
    propagate_state(m.fs.M1_to_RO)

    try:
        m.fs.RO.initialize()
    except:
        pass

    # RO brine to P2
    propagate_state(m.fs.RO_permeate_to_product) 
    propagate_state(m.fs.RO_retentate_to_P2)
    # P2 initialize
    m.fs.P2.initialize()



def solve(blk, solver=None, tee=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    return results


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

    # initial_conditions = {
    #     "feed_flow": 2.92,
    #     "feed_conc": 3.4,
    #     "reject_flow": 45.25,
    #     "reject_conc_start": 3.9,
    #     "temperature_start": 25,  # degC
    #     "p1_pressure_start": 306,  # psi
    #     "A_comp": 4.422e-12,
    #     "B_comp": 5.613e-8,
    #     "channel_height": 0.0008636,
    #     "spacer_porosity": 0.7081,
    #     "water_recovery": 0.063,
    #     "n_time_points": 5,
    # }

    m = build_system()
    dt = DiagnosticsToolbox(m)
    set_operating_conditions(m)
    scale_system(m)
    # autoscaler = AutoScaler()
    # autoscaler.scale_model(m,descend_into=True)
    initialize_system(m)
    m.fs.report()
    m.fs.RO.report()

    m.fs.RO.feed_side.N_Re.setlb(None)
    solver =get_solver()

    solver.options["max_iter"] = 1
    res = solver.solve(m,tee=True)
    
    # assert_optimal_termination(res)
    # dt.compute_infeasibility_explanation()
    # ccro.mp_df