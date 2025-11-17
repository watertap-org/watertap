#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

__author__ = "Alexander Dudchenko"

from pyomo.environ import (
    ConcreteModel,
    value,
    units as pyunits,
    TransformationFactory,
    assert_optimal_termination,
)

from pyomo.network import Arc

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state
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
import yaml
from idaes.models.unit_models import Feed
import os


def fit_bw30_4040_to_spec_sheet(save_location=None):
    # https://www.dupont.com/content/dam/water/amer/us/en/water/public/documents/en/RO-FilmTec-BW30-PRO-4040-and-BW30-PRO-2540-PDS-45-D03970-en.pdf
    return fit_ro_module_to_spec_sheet(
        membrane_name="BW30 PRO-4040",
        water_production_rate=9.8 * pyunits.m**3 / pyunits.day,
        nacl_rejection=99.7 * pyunits.percent,  # percent
        feed_conc=2000 * pyunits.mg / pyunits.liter,
        recovery=0.15,
        module_length=1 * pyunits.m,
        pressure=225 * pyunits.psi,
        membrane_area=7.9 * pyunits.m**2,
        channel_height=1 * pyunits.mm,
        spacer_porosity=0.95,
        temperature=25,  # degrees C
        save_location=save_location,
    )


def fit_sw30_4040_to_spec_sheet(save_location=None):
    # https://www.dupont.com/content/dam/water/amer/us/en/water/public/documents/en/RO-FilmTec-SW30-Seawater-PDS-45-D01519-en.pdf
    return fit_ro_module_to_spec_sheet(
        membrane_name="SW30-4040",
        water_production_rate=7.4 * pyunits.m**3 / pyunits.day,
        nacl_rejection=99.7 * pyunits.percent,  # percent
        feed_conc=32000 * pyunits.mg / pyunits.liter,
        recovery=0.08,
        module_length=1 * pyunits.m,
        pressure=55 * pyunits.bar,
        membrane_area=7.9 * pyunits.m**2,
        channel_height=1 * pyunits.mm,
        spacer_porosity=0.95,
        temperature=25,  # degrees C
        save_location=save_location,
    )


def fit_ro_module_to_spec_sheet(
    membrane_name="RO_module",
    water_production_rate=38 * pyunits.m**3 / pyunits.day,
    nacl_rejection=99.5 * pyunits.percent,  # percent
    feed_conc=1500 * pyunits.mg / pyunits.liter,
    recovery=0.15,
    module_length=1 * pyunits.m,
    pressure=150 * pyunits.psi,
    membrane_area=40 * pyunits.m**2,
    channel_height=1e-3 * pyunits.m,
    spacer_porosity=0.95,
    temperature=25,  # degrees C
    save_location=None,
    save_to_yaml=True,
):
    """Method for fitting RO spec sheet data to 1D RO model to find A and B membrane parameters. Use pyunits to specify units
    for each input.

    Args:
        membrane_name : str
            Name of the membrane module.
        water_production_rate : float
            Water production rate (volume flow rate) [m3/day].
        nacl_rejection : float
            Salt rejection [%].
        feed_conc : float
            Feed concentration [mg/L].
        recovery : float
            Recovery [%].
        module_length : float
            Length of the RO module [m].
        pressure : float
            Feed pressure [Pa].
        membrane_area : float
            Membrane area [m2].
        channel_height : float
            Channel height [m].
        spacer_porosity : float
            Spacer porosity [-].
        temperature : float
            Temperature [deg C].
        save_location : str
            Directory to save the fitted membrane parameters yaml file.
        save_to_yaml : bool

    Output:
        Dictionary that contains fitted membrane parameters A and B along with other
        membrane design parameters.
    """

    solver = get_solver()

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()
    # seem feed and ro model
    m.fs.feed = Feed(property_package=m.fs.properties)
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

    m.fs.feed_to_ro = Arc(source=m.fs.feed.outlet, destination=m.fs.RO.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    # specify feed conditions
    m.fs.feed.properties[0].flow_vol_phase["Liq"].fix(water_production_rate / recovery)
    m.fs.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].fix(feed_conc)
    m.fs.feed.properties[0].temperature.fix(temperature + 273.15)
    m.fs.feed.properties[0].pressure.fix(pressure)

    print("Degrees of freedom before initialization: ", degrees_of_freedom(m.fs.feed))
    assert degrees_of_freedom(m.fs.feed) == 0
    # Solve the feed to get mass flows and concentrations through out
    result = solver.solve(m.fs.feed)
    assert_optimal_termination(result)
    propagate_state(m.fs.feed_to_ro)

    # configure RO
    m.fs.RO.length.fix(module_length)
    iscale.set_scaling_factor(m.fs.RO.length, value(1 / m.fs.RO.length))
    m.fs.RO.area.fix(membrane_area)
    iscale.set_scaling_factor(m.fs.RO.area, value(1 / m.fs.RO.area))
    m.fs.RO.feed_side.channel_height.fix(channel_height)
    m.fs.RO.feed_side.spacer_porosity.fix(spacer_porosity)
    # initial guess for intialization
    m.fs.RO.A_comp.fix(4.2e-12)
    m.fs.RO.B_comp.fix(3.5e-8)
    m.fs.RO.permeate.pressure[0].fix(101325)  # 1 atm

    # scale mass flow units
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        value(1 / m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]),
        index=("Liq", "H2O"),
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp",
        value(1 / m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"]),
        index=("Liq", "NaCl"),
    )
    iscale.calculate_scaling_factors(m)
    # initialization guess
    m.fs.RO.initialize()

    print("Degrees of freedom before RO box solve: ", degrees_of_freedom(m.fs.feed))
    assert degrees_of_freedom(m) == 0
    results = solver.solve(m, tee=False)
    assert_optimal_termination(results)
    # now fix recovery and rejection to spec sheet values and unfix A and B to solve for them
    m.fs.RO.A_comp.unfix()
    m.fs.RO.B_comp.unfix()
    m.fs.RO.rejection_phase_comp[0, "Liq", "NaCl"].fix(nacl_rejection)
    m.fs.RO.recovery_vol_phase[0.0, "Liq"].fix(recovery)

    print("Degrees of freedom before A/B solve: ", degrees_of_freedom(m.fs.feed))
    assert degrees_of_freedom(m) == 0
    results = solver.solve(m, tee=True)
    assert_optimal_termination(results)
    print('Fit successfully completed for membrane: "', membrane_name, '"')
    membrane_design = {
        "A": {
            "value": m.fs.RO.A_comp[0, "H2O"].value,
            "units": str(m.fs.RO.A_comp[0, "H2O"].get_units()),
        },
        "B": {
            "value": m.fs.RO.B_comp[0, "NaCl"].value,
            "units": str(m.fs.RO.B_comp[0, "NaCl"].get_units()),
        },
        "A (LMH/bar)": {
            "value": value(
                pyunits.convert(
                    m.fs.RO.A_comp[0, "H2O"],
                    to_units=pyunits.L / (pyunits.m**2 * pyunits.hr * pyunits.bar),
                )
            ),
            "units": "LMH/bar",
        },
        "B (LMH)": {
            "value": value(
                pyunits.convert(
                    m.fs.RO.B_comp[0, "NaCl"],
                    to_units=pyunits.L / (pyunits.m**2 * pyunits.hr),
                )
            ),
            "units": "LMH",
        },
        "Area": {"value": m.fs.RO.area.value, "units": str(m.fs.RO.area.get_units())},
        "Length": {
            "value": m.fs.RO.length.value,
            "units": str(m.fs.RO.length.get_units()),
        },
        "Porosity": {
            "value": m.fs.RO.feed_side.spacer_porosity.value,
            "units": str(m.fs.RO.feed_side.spacer_porosity.get_units()),
        },
        "Channel_height": {
            "value": m.fs.RO.feed_side.channel_height.value,
            "units": str(m.fs.RO.feed_side.channel_height.get_units()),
        },
        "Rejection_NaCl": {
            "value": m.fs.RO.rejection_phase_comp[0, "Liq", "NaCl"].value,
            "units": str(m.fs.RO.rejection_phase_comp[0, "Liq", "NaCl"].get_units()),
        },
        "Recovery": {
            "value": m.fs.RO.recovery_vol_phase[0.0, "Liq"].value,
            "units": str(m.fs.RO.recovery_vol_phase[0.0, "Liq"].get_units()),
        },
    }
    for key, val in membrane_design.items():
        print(f"{key}: {val}")
    if save_to_yaml:
        save_dir = save_location if save_location is not None else os.getcwd()
        os.makedirs(save_dir, exist_ok=True)

        print("Saving membrane design to directory: ", save_dir)
        with open(f"{save_dir}/{membrane_name}.yaml", "w") as outfile:
            yaml.dump(
                membrane_design, outfile, default_flow_style=False, sort_keys=False
            )
    return membrane_design


if __name__ == "__main__":
    # fit_bw30_4040_to_spec_sheet()
    fit_sw30_4040_to_spec_sheet()
