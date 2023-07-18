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

import pytest
import idaes.logger as idaeslog
from pyomo.environ import (
    ConcreteModel,
    value,
    Var,
    Constraint,
    assert_optimal_termination,
    units as pyunits,
    NonNegativeReals,
)
from pyomo.util.check_units import assert_units_consistent
from pyomo.network import Port
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    FlowDirection,
)
from watertap.unit_models.osmotically_assisted_reverse_osmosis_0D import (
    OsmoticallyAssistedReverseOsmosis0D,
)
import watertap.property_models.NaCl_prop_pack as props

from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    badly_scaled_var_generator,
)

from watertap.core import (
    MembraneChannel0DBlock,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
    FrictionFactor,
)

import matplotlib.pyplot as plt
import numpy as np
from debug import *

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


def main():
    # set up solver
    # solver = get_solver()

    # build, set, and initialize
    m = build(water_recovery=0.5)
    results = solve(m)
    # set_operating_conditions(m, number_of_stages=number_of_stages)
    # initialize_system(m, number_of_stages, solver=solver)
    assert_optimal_termination(results)
    display_state(m)
    display_design(m)
    plot(m)

    # print("\n***---Optimization results---***")
    # display_system(m)
    # display_design(m)
    # if erd_type == ERDtype.pump_as_turbine:
    #     display_state(m)
    # else:
    #     pass
    return m

def solve(m):
    solver = get_solver()
    try:
        results = solver.solve(m, tee=True)
        debug(m)
    except:
        debug(m, automate_rescale=True, resolve=True)

    return results

def build(water_recovery=0.5):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = OsmoticallyAssistedReverseOsmosis0D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        pressure_change_type=PressureChangeType.calculated,
    )

    # fully specify system
    feed_flow_mass = 5 / 18
    feed_mass_frac_NaCl = 0.075
    feed_pressure = 65e5
    feed_temperature = 273.15 + 25
    membrane_area = 50
    A = 1e-12
    B = 7.7e-8
    pressure_atmospheric = 101325

    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
    m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.feed_inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.feed_inlet.temperature[0].fix(feed_temperature)

    permeate_flow_mass = 0.33 * feed_flow_mass
    permeate_mass_frac_NaCl = 0.1
    permeate_mass_frac_H2O = 1 - permeate_mass_frac_NaCl
    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        permeate_flow_mass * permeate_mass_frac_H2O
    )
    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        permeate_flow_mass * permeate_mass_frac_NaCl
    )
    m.fs.unit.permeate_inlet.pressure[0].fix(1.5e5)
    m.fs.unit.permeate_inlet.temperature[0].fix(feed_temperature)

    m.fs.unit.area.fix(membrane_area)

    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)

    m.fs.unit.structural_parameter.fix(1200e-6)

    m.fs.unit.permeate_side.channel_height.fix(0.002)
    m.fs.unit.permeate_side.spacer_porosity.fix(0.89)
    m.fs.unit.feed_side.channel_height.fix(0.002)
    m.fs.unit.feed_side.spacer_porosity.fix(0.95)
    # m.fs.unit.feed_side.velocity[0, 0].fix(0.1)
    m.fs.unit.feed_side.N_Re[0, 0].fix(375)

    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    calculate_scaling_factors(m)

    print(f"DOF: {degrees_of_freedom(m)}")

    m.fs.unit.initialize()

    m.fs.mass_water_recovery = Var(
        initialize=0.5,
        bounds=(0, 1),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="System Volumetric Recovery of Water",
    )
    m.fs.eq_mass_water_recovery = Constraint(
        expr= m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "H2O"] * m.fs.mass_water_recovery == 
        (m.fs.unit.permeate_outlet.flow_mass_phase_comp[0, "Liq", "H2O"] - 
        m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"])
    )

    m.fs.unit.permeate_inlet.pressure[0].unfix()

    m.fs.unit.feed_side.velocity[0, 0].unfix()
    m.fs.unit.feed_side.velocity[0, 0].setlb(0)
    m.fs.unit.feed_side.velocity[0, 0].setub(1)

    m.fs.unit.area.unfix()

    m.fs.mass_water_recovery.fix(water_recovery)
    m.fs.unit.permeate_outlet.pressure[0].fix(1e5)
    m.fs.unit.feed_side.N_Re[0, 0].fix(375)

    print(f"DOF: {degrees_of_freedom(m)}")

    return m


def display_design(m):
    print("--decision variables--")
    print(
        "OARO Stage feed side water flux: %.1f L/m2/h"
        % (
            value(m.fs.unit.flux_mass_phase_comp[0, 0, "Liq", "H2O"])
            / 1e3
            * 1000
            * 3600,
        )
    )
    print(
        "OARO permeate side water flux: %.1f L/m2/h"
        % (
            value(m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "H2O"])
            / 1e3
            * 1000
            * 3600,
        )
    )
    print(
        "OARO average water flux: %.1f L/m2/h"
        % (
            value(m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"])
            / 1e3
            * 1000
            * 3600,
        )
    )
    print(
        "OARO average salt flux: %.1f g/m2/h"
        % (value(m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]) * 1000 * 3600,)
    )
    print(
        "OARO feed operating pressure: %.1f bar"
        % (m.fs.unit.feed_inlet.pressure[0].value / 1e5)
    )
    print(
        "OARO feed side pressure drop: %.1f bar"
        % (-m.fs.unit.feed_side.deltaP[0].value / 1e5)
    )
    print(
        "OARO permeate operating pressure: %.1f bar"
        % (m.fs.unit.permeate_inlet.pressure[0].value / 1e5)
    )
    print(
        "OARO permeate side pressure drop: %.1f bar"
        % (-m.fs.unit.permeate_side.deltaP[0].value / 1e5)
    )
    print("OARO membrane area:      %.1f m2" % (m.fs.unit.area.value))
    print("OARO membrane width:      %.1f m" % (m.fs.unit.width.value))
    print("OARO membrane length:      %.1f m" % (m.fs.unit.length.value))
    print(
        "OARO feed side average Reynolds number: %.1f"
        % value(m.fs.unit.feed_side.N_Re_avg[0])
    )
    print(
        "OARO permeate side average Reynolds number: %.1f"
        % value(m.fs.unit.permeate_side.N_Re_avg[0])
    )
    print(
        "OARO feed side average mass transfer coeff.: %.1f mm/h"
        % value(m.fs.unit.feed_side.K_avg[0, "NaCl"] * 1000 * 3600)
    )
    print(
        "OARO permeate side average mass transfer coeff.: %.1f mm/h"
        % value(m.fs.unit.permeate_side.K_avg[0, "NaCl"] * 1000 * 3600)
    )
    print(
        "OARO water perm. coeff.:  %.3f LMH/bar"
        % (m.fs.unit.A_comp[0, "H2O"].value * (3.6e11))
    )
    print(
        "OARO salt perm. coeff.:  %.3f LMH/bar"
        % (m.fs.unit.B_comp[0, "NaCl"].value * (1000.0 * 3600.0))
    )


def display_state(m):
    print("--------state---------")

    def print_state(s, b):
        feed_flow_mass = (
            sum(
                m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", j].value
                for j in ["H2O", "NaCl"]
            )
            * 3600
        )
        flow_mass = (
            sum(b.flow_mass_phase_comp[0, "Liq", j].value for j in ["H2O", "NaCl"])
            * 3600
        )
        normalized_flow_mass = flow_mass / feed_flow_mass * 100
        mass_frac_ppm = (
            b.flow_mass_phase_comp[0, "Liq", "NaCl"].value / (flow_mass / 3600) * 1e3
        )
        pressure_bar = b.pressure[0].value / 1e5
        print(
            s.ljust(20)
            + ": %.2f kg/h,  %.0f, %.3f g/L, %.1f bar"
            % (flow_mass, normalized_flow_mass, mass_frac_ppm, pressure_bar)
        )

    print_state(f"OARO feed inlet:", m.fs.unit.feed_inlet)
    print_state(f"OARO permeate inlet:", m.fs.unit.permeate_inlet)
    print_state(f"OARO feed outlet:", m.fs.unit.feed_outlet)
    print_state(f"OARO permeate outlet:", m.fs.unit.permeate_outlet)


def plot(m):
    feed_conc_in = (
        m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].value
        / sum(
            m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", j].value
            for j in ["H2O", "NaCl"]
        )
        * 1e3
    )
    feed_conc_out = (
        m.fs.unit.feed_outlet.flow_mass_phase_comp[0, "Liq", "NaCl"].value
        / sum(
            m.fs.unit.feed_outlet.flow_mass_phase_comp[0, "Liq", j].value
            for j in ["H2O", "NaCl"]
        )
        * 1e3
    )
    permeate_conc_in = (
        m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].value
        / sum(
            m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", j].value
            for j in ["H2O", "NaCl"]
        )
        * 1e3
    )
    permeate_conc_out = (
        m.fs.unit.permeate_outlet.flow_mass_phase_comp[0, "Liq", "NaCl"].value
        / sum(
            m.fs.unit.permeate_outlet.flow_mass_phase_comp[0, "Liq", j].value
            for j in ["H2O", "NaCl"]
        )
        * 1e3
    )
    feed_flux = (
        value(m.fs.unit.flux_mass_phase_comp[0, 0, "Liq", "H2O"]) / 1e3 * 1000 * 3600
    )
    permeate_flux = (
        value(m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "H2O"]) / 1e3 * 1000 * 3600
    )
    
    salt_flux_in = (
        value(pyunits.convert(m.fs.unit.flux_mass_phase_comp[0, 0, "Liq", "NaCl"], to_units=pyunits.gram / pyunits.m**2 / pyunits.hour))
    )
    salt_flux_out = (
        value(pyunits.convert(m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "NaCl"], to_units=pyunits.gram / pyunits.m**2 / pyunits.hour))
    )
    xpoints = np.array([0, 1])
    ypoints1 = np.array([feed_conc_in, feed_conc_out])
    ypoints2 = np.array([permeate_conc_out, permeate_conc_in])
    ypoints3 = np.array([feed_flux, permeate_flux])
    ypoints4 = np.array([salt_flux_out, salt_flux_in])

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(xpoints, ypoints1, "k")
    ax.plot(xpoints, ypoints2, "k--")
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 175])
    ax.set_ylabel("Concentration (g/L)", fontsize=12)
    ax.legend(
        ["feed side concentration", "permeate side concentration"], loc="upper left"
    )

    ax.set_xlabel("Normalized Membrane Length", fontsize=12)

    ax.tick_params(axis="x", labelsize=12)
    ax.tick_params(axis="y", labelsize=12)
    plt.locator_params(axis="y", nbins=8)

    ax2 = ax.twinx()
    ax2.plot(xpoints, ypoints3)
    ax2.set_ylim([0, 10])
    ax2.set_ylabel("Water flux (LMH)", fontsize=12)
    ax2.legend(["water flux"])
    ax2.tick_params(axis="x", labelsize=12)
    ax2.tick_params(axis="y", labelsize=12)
    ax2.yaxis.label.set_color('#1f77b4')
    ax2.spines["right"].set_color('#1f77b4')
    ax2.tick_params(axis="y", colors='#1f77b4')

    ax3 = ax.twinx()
    ax3.spines.right.set_position(("axes", 1.25))
    ax3.plot(xpoints, ypoints4, color="#c07432")
    ax3.set_ylabel("Salt Flux (kg/m2-hr)", fontsize=12)
    ax3.set_ylim([0, 20])
    ax3.yaxis.label.set_color('#c07432')
    ax3.spines["right"].set_color("#c07432")
    ax3.tick_params(axis="y", colors="#c07432")

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    m = main()

    # NOTE: I think this can better match the paper if the salt flux is increased at the beginning of the permeate side and reduced at the beginning of the feed side
