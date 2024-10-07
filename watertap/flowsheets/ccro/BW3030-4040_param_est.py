import matplotlib.pyplot as plt
from pyomo.environ import (
    ConcreteModel,
    check_optimal_termination,
    value,
    assert_optimal_termination,
    units as pyunits,
)

import pyomo.contrib.parmest.parmest as parmest
from pyomo.contrib.parmest.experiment import Experiment

from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import FlowsheetBlock
import idaes.core.util.scaling as iscale

from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D as RO1D,
    PressureChangeType,
    MassTransferCoefficient,
    ConcentrationPolarizationType,
)
from watertap.core.util.model_diagnostics.infeasible import *
from watertap.property_models import NaCl_T_dep_prop_pack as props
from watertap.core.solvers import get_solver

solver = get_solver()

#  Spec sheet, data for BW30-4040: https://www.lenntech.com/Data-sheets/DuPont-FilmTec-BW30-4040-L.pdf

# Permeate flow and salt rejection based on the following test conditions:
# feed pressure = 225 psi
# feed water concentration = 2000 mg/L
# salt rejection = 99.5% 
# water recovery = 15% 
# membrane area = 37.2 m2
# perm flow = 9.1 m3/d

rho = 997.0 * pyunits.kg / pyunits.m**3

feed_conc = 2 * pyunits.g / pyunits.liter
recovery = 0.15
pressure = 225 * pyunits.psi
perm_vol_flow = 9.1 * pyunits.m**3 / pyunits.day
salt_rej = 0.995
temperature = 25 

mem_length = (1.016 - 2 * 0.0267) * pyunits.m
mem_area = 7.2 * pyunits.m**2
pressure_loss = -15 * pyunits.psi

perm_mass_flow = pyunits.convert(rho * perm_vol_flow, to_units=pyunits.kg / pyunits.s)

feed_vol_flow = pyunits.convert(
    perm_vol_flow / recovery, to_units=pyunits.m**3 / pyunits.s
)
feed_mass_flow_water = value(perm_mass_flow / recovery)
feed_mass_flow_salt = value(
    pyunits.convert(feed_vol_flow * feed_conc, to_units=pyunits.kg / pyunits.s)
)

spacer_thickness = 34 # mil
channel_height = spacer_thickness * 2.54e-5 # mil to m

def solve(m):
    # ---solving---
    solver = get_solver()

    print("\n--------- SOLVING ---------\n")
    results = solver.solve(m)

    if check_optimal_termination(results):
        print("\n--------- OPTIMAL SOLVE!!! ---------\n")
        print(f"water_perm = {m.fs.RO.A_comp[0,'H2O']()} m/s/Pa")
        print(f"salt_perm = {m.fs.RO.B_comp[0,'NaCl']()} m/s")
        print(f"porosity = {m.fs.RO.feed_side.spacer_porosity()}")
        print(f"channel_height = {m.fs.RO.feed_side.channel_height()}")
        print("\n")

        return results
    
    assert False


def estimate_params(
    water_perm=4.2e-12,
    salt_perm=3.5e-8,
    porosity=0.95,
):

    global m

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()

    m.fs.RO = RO1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        transformation_scheme="BACKWARD",
        transformation_method="dae.finite_difference",
        module_type="spiral_wound",
        finite_elements=10,
        has_full_reporting=True,
    )

    print("RO DOF:", degrees_of_freedom(m.fs.RO))

    m.fs.RO.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(feed_mass_flow_salt)
    m.fs.RO.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(feed_mass_flow_water)
    m.fs.RO.inlet.pressure[0].fix(pressure)
    m.fs.RO.inlet.temperature[0].fix(temperature + 273.15)

    m.fs.RO.permeate.pressure[0].fix(101325)
    # m.fs.RO.feed_side.channel_height.fix(1e-3)
    m.fs.RO.feed_side.channel_height.fix(channel_height)
    m.fs.RO.length.fix(mem_length)

    m.fs.RO.area.fix(mem_area)
    m.fs.RO.A_comp.fix(water_perm)
    m.fs.RO.B_comp.fix(salt_perm)
    m.fs.RO.feed_side.spacer_porosity.fix(porosity)
    m.fs.RO.feed_side.spacer_porosity.setlb(0.65)

    print("DOF = ", degrees_of_freedom(m))
    print("RO DOF = ", degrees_of_freedom(m.fs.RO))
    # assert_no_degrees_of_freedom(m)
    iscale.set_scaling_factor(m.fs.RO.area, 1e-2)
    iscale.set_scaling_factor(m.fs.RO.feed_side.area, 1e-2)
    iscale.set_scaling_factor(m.fs.RO.width, 1e-2)

    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    iscale.calculate_scaling_factors(m)

    m.fs.RO.initialize()
    results = solve(m)

    # Unfix A variable
    print("unfix A, fix perm flow...")
    m.fs.RO.A_comp.unfix()

    # Fix the permeate flow
    m.fs.RO.mixed_permeate[0.0].flow_mass_phase_comp["Liq", "H2O"].fix(perm_mass_flow)

    print("DOF = ", degrees_of_freedom(m))

    results = solve(m)

    # Unfix B variable
    print("unfix B, fix rejection...")
    m.fs.RO.B_comp.unfix()

    # Fix the salt rejection
    m.fs.RO.rejection_phase_comp[0, "Liq", "NaCl"].fix(salt_rej)

    print("DOF = ", degrees_of_freedom(m))

    results = solve(m)

    # Unfix the permeate flow rate and concentration
    m.fs.RO.mixed_permeate[0.0].flow_mass_phase_comp["Liq", "H2O"].unfix()
    m.fs.RO.rejection_phase_comp[0, "Liq", "NaCl"].unfix()

    # Fix the new A & B values. This will fix them at the solved values from the previous steps
    m.fs.RO.A_comp.fix()
    m.fs.RO.B_comp.fix()

    # Fix the new flow rate
    m.fs.RO.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(feed_mass_flow_water)
    # Fix the pressure drop
    m.fs.RO.deltaP.fix(pressure_loss)

    # Unfix the spacer porosity
    print("unfix spacer porosity...")
    m.fs.RO.feed_side.spacer_porosity.unfix()

    print("DOF = ", degrees_of_freedom(m))
    # m.fs.RO.feed_side.channel_height.unfix()
    # print("\n")
    # print("DOF = ", degrees_of_freedom(m))

    results = solve(m)

    water_perm = m.fs.RO.A_comp[0, "H2O"]()
    salt_perm = m.fs.RO.B_comp[0, "NaCl"]()
    porosity = m.fs.RO.feed_side.spacer_porosity()

    return water_perm, salt_perm, porosity


if __name__ == "__main__":
    wp = list()
    sp = list()
    por = list()
    num = list()

    water_perm = 4.2e-12
    salt_perm = 3.5e-8
    porosity = 0.9
    n = 0

    # starting point
    wp.append(water_perm)
    sp.append(salt_perm)
    por.append(porosity)
    num.append(n)
    n += 1

    num_iter = 5

    while n != num_iter:

        water_perm, salt_perm, porosity = estimate_params(
            water_perm=water_perm, salt_perm=salt_perm, porosity=porosity
        )
        wp.append(water_perm)
        sp.append(salt_perm)
        por.append(porosity)
        num.append(n)
        # print(n, water_perm, salt_perm, porosity)
        n += 1

    _, ax = plt.subplots()
    ax.scatter(num, wp, marker="o", color="r")
    ax.set_title("water perm")
    ax.set_xlabel("iteration")

    _, ax = plt.subplots()
    ax.scatter(num, sp, marker="o", color="g")
    ax.set_title("salt perm")
    ax.set_xlabel("iteration")

    _, ax = plt.subplots()
    ax.scatter(num, por, marker="o", color="k")
    ax.set_title("porosity")
    ax.set_xlabel("iteration")

    plt.show()
