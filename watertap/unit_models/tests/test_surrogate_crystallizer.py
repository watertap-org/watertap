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

import os
import pytest
from pyomo.environ import (
    ConcreteModel,
    Var,
    Constraint,
    units as pyunits,
    assert_optimal_termination,
)
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock
from idaes.core.util.exceptions import ConfigurationError
from watertap.property_models.water_prop_pack import WaterParameterBlock
from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
    MaterialFlowBasis,
)
from idaes.core.util.model_statistics import degrees_of_freedom

from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core.surrogate.pysmo_surrogate import (
    PysmoSurrogate,
)
from watertap.unit_models.surrogate_crystallizer import SurrogateCrystallizer

from idaes.core import UnitModelCostingBlock
from watertap.costing import WaterTAPCosting
from watertap.core.solvers import get_solver

__author__ = "Oluwamayowa Amusat, Adam Atia"


# TODO: can consider adding NN model surrogate function and test
def add_crystallizer_rbf_model(
    blk,
    surrogate_inputs_with_bounds,
    surrogate_outputs,
    surrogate_to_flowsheet_basis_ratio,
):
    ##############################################################################################################
    # block loading into IDAES Surrogate Block
    ##############################################################################################################
    crystallizer_inputs = [
        surrogate_inputs_with_bounds[k]["flowsheet_var"]
        for k in surrogate_inputs_with_bounds.keys()
    ]
    crystallizer_outputs = [
        surrogate_outputs[k]["flowsheet_var"] for k in surrogate_outputs.keys()
    ]

    try:
        assert set(list(blk.solids_list)).issubset(surrogate_outputs.keys())
    except AssertionError:
        raise ConfigurationError(
            "List of solids provided in crystallizer model must match surrogate output keys for solids. Please check."
        )

    filename = [
        "Vapor_Pressure",
        "Calcite_g",
        "Anhydrite_g",
        "Glauberite_g",
        "Halite_g",
    ]

    for sm in range(0, len(filename)):
        block_name = "crystallizer_surrogate" + "_" + filename[sm]
        blk.add_component(block_name, SurrogateBlock(concrete=True))
        surrogate_name = f"{filename[sm]}.json"
        surrogate_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..",
            "..",
            "data",
            "surrogate_defaults",
            "surrogate_crystallizer_defaults",
            surrogate_name,
        )
        current_surrogate = PysmoSurrogate.load_from_file(surrogate_path)

        getattr(blk, block_name).build_model(
            current_surrogate,
            input_vars=crystallizer_inputs,
            output_vars=crystallizer_outputs[sm],
        )

    # Add constraints tying provided mass variables with crystallizer unit variables - done to avoid unit challenges
    def eq_mass_conversion_constraint(b, j):
        return (
            b.flow_mass_sol_comp_apparent[j]
            == pyunits.convert(
                surrogate_outputs[j]["flowsheet_var"], to_units=pyunits.kg / pyunits.s
            )
            / surrogate_to_flowsheet_basis_ratio
        )

    blk.mass_conversion_constraints = Constraint(
        blk.solids_list, rule=eq_mass_conversion_constraint
    )

    # Add a constraint to tie pressures together after unit conversion
    blk.pressure_conversion_constraints = Constraint(
        expr=blk.pressure_operating
        == pyunits.convert(
            surrogate_outputs["Vapor Pressure (atm)"]["flowsheet_var"],
            to_units=pyunits.Pa,
        )
    )


@pytest.mark.component
def test_rbf_surrogate():
    m = ConcreteModel()
    m.case = "BGW1"
    m.fs = FlowsheetBlock(dynamic=False)

    input_ions = ["Cl_-", "Na_+", "SO4_2-", "Mg_2+", "Ca_2+", "K_+", "HCO3_-"]
    solids_list = {
        "Calcite_g": {"Ca_2+": 1, "HCO3_-": 1},
        "Anhydrite_g": {"Ca_2+": 1, "SO4_2-": 1},
        "Glauberite_g": {"Ca_2+": 1, "Na_+": 2, "SO4_2-": 2},
        "Halite_g": {"Na_+": 1, "Cl_-": 1},
    }
    m.fs.cryst_prop_feed = MCASParameterBlock(
        solute_list=input_ions,
        material_flow_basis=MaterialFlowBasis.mass,
    )

    m.fs.water_properties_vapor = WaterParameterBlock()

    # Create crystallizer framework
    m.fs.cryst = SurrogateCrystallizer(
        property_package=m.fs.cryst_prop_feed,
        vapor_property_package=m.fs.water_properties_vapor,
        solids_ions_dict=solids_list,
    )

    # Feed composition
    feed = {
        "ion_composition_g_kg": {
            "Na_+": 91.49859876,
            "K_+": 3.185931759,
            "Cl_-": 159.0941278,
            "Ca_2+": 0.788476605,
            "Mg_2+": 10.5806906,
            "HCO3_-": 1.614773395,
            "SO4_2-": 22.20893266,
        },
        "water_content_kg": 12.0764972138554,
    }

    g_to_kg = 1e-3
    for i in m.fs.cryst_prop_feed.component_list:
        if i == "H2O":
            m.fs.cryst.inlet.flow_mass_phase_comp[0.0, "Liq", i].fix(
                feed["water_content_kg"]
            )
        else:
            m.fs.cryst.inlet.flow_mass_phase_comp[0.0, "Liq", i].fix(
                feed["water_content_kg"] * feed["ion_composition_g_kg"][i] * g_to_kg
            )

    m.fs.cryst.inlet.temperature.fix(298.15)
    m.fs.cryst.inlet.pressure.fix(101325)

    assert_units_consistent(m)

    # Create dummy variables for inputs and outputs
    surrogate_inputs = {
        "Temperature": {
            "flowsheet_var": m.fs.cryst.temperature_operating,
            "var_bounds": (303.15, 372.15),
        },
        "Evaporation percent": {
            "flowsheet_var": m.fs.cryst.evaporation_percent,
            "var_bounds": (30.0, 98.0),
        },
    }
    m.fs.o1 = Var(units=pyunits.atm)
    m.fs.o2 = Var(units=pyunits.g / pyunits.s)
    m.fs.o3 = Var(units=pyunits.g / pyunits.s)
    m.fs.o4 = Var(units=pyunits.g / pyunits.s)
    m.fs.o5 = Var(units=pyunits.g / pyunits.s)
    surrogate_outputs = {
        "Vapor Pressure (atm)": {"flowsheet_var": m.fs.o1, "Solid": False},
        "Calcite_g": {"flowsheet_var": m.fs.o2, "Solid": True},
        "Anhydrite_g": {"flowsheet_var": m.fs.o3, "Solid": True},
        "Glauberite_g": {"flowsheet_var": m.fs.o4, "Solid": True},
        "Halite_g": {"flowsheet_var": m.fs.o5, "Solid": True},
    }

    add_crystallizer_rbf_model(
        m.fs.cryst,
        surrogate_inputs_with_bounds=surrogate_inputs,
        surrogate_outputs=surrogate_outputs,
        surrogate_to_flowsheet_basis_ratio=1,
    )

    # Costing
    m.fs.costing = WaterTAPCosting()
    m.fs.cryst.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
    )
    m.fs.costing.cost_process()

    # Edit bound on liquid pressures from property package
    m.fs.cryst.properties_out_liq[0].pressure.setlb(1000)

    # 1. Simulate single case
    m.fs.cryst.temperature_operating.fix(40 + 273.15)
    m.fs.cryst.evaporation_percent.fix(60)

    # Initialize and solve
    assert degrees_of_freedom(m) == 0
    m.fs.cryst.initialize()
    solver = get_solver()
    res = solver.solve(m)
    assert_optimal_termination(res)
    m.fs.cryst.report()
