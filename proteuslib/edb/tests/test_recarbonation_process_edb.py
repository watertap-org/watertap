###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################
"""
This file is to run the example of running a phase change for the dissolution
of CO2 into water.
"""
from pprint import pformat, pprint
import pytest

# Import specific pyomo objects
from pyomo.environ import (
    ConcreteModel,
    SolverStatus,
    TerminationCondition,
    value,
    Suffix,
)

# Import pyomo methods to check the system units
from pyomo.util.check_units import assert_units_consistent
from pyomo.environ import units as pyunits

from idaes.core.util import scaling as iscale

from idaes.core.util import get_solver

# Import idaes methods to check the model during construction
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    fixed_variables_set,
    activated_constraints_set,
)

# Import the idaes objects for Generic Properties and Reactions
from idaes.generic_models.properties.core.generic.generic_property import (
    GenericParameterBlock,
)
from idaes.generic_models.properties.core.generic.generic_reaction import (
    GenericReactionParameterBlock,
    ConcentrationForm,
)


from idaes.generic_models.properties.core.reactions.dh_rxn import constant_dh_rxn
from idaes.generic_models.properties.core.reactions.equilibrium_constant import gibbs_energy
from idaes.generic_models.properties.core.reactions.equilibrium_constant import (
    van_t_hoff,
)

# Import safe log power law equation
from idaes.generic_models.properties.core.reactions.equilibrium_forms import (
    log_power_law_equil,
)

# Import the idaes object for the EquilibriumReactor unit model
from idaes.generic_models.unit_models.equilibrium_reactor import EquilibriumReactor
from idaes.generic_models.properties.core.pure.Perrys import Perrys
from idaes.generic_models.properties.core.pure.NIST import NIST

# Import the core idaes objects for Flowsheets and types of balances
from idaes.core import FlowsheetBlock

# Import statements to be used in the starter config dict
from idaes.core import VaporPhase, AqueousPhase
from idaes.core.components import Solvent, Solute, Cation, Anion
from idaes.core.phases import PhaseType as PT
from idaes.generic_models.properties.core.phase_equil.forms import fugacity
from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.ideal import Ideal
from idaes.generic_models.properties.core.phase_equil import SmoothVLE
from idaes.generic_models.properties.core.phase_equil.bubble_dew import IdealBubbleDew

# Import log10 function from pyomo
from pyomo.environ import log10

from proteuslib.edb import ElectrolyteDB
from . import data as test_data

__author__ = "Srikanth Allu, Dan Gunter"


# Try to connect to DB at import time, so we can use the state of that connection
# to determine whether to skip tests
electrolyte_db = None
try:
    electrolyte_db = ElectrolyteDB(db="edb")
except Exception as err:
    print(f"Failed to connect to ElectrolyteDB: {err}")

if electrolyte_db is None:
    pytest.skip("Could not connect to electrolyte database", allow_module_level=True)


component_names = ["H2O", "CO2", "H +", "OH -", "H2CO3", "HCO3 -", "CO3 2-"]


@pytest.fixture
def thermo_config():
    thermo_base = next(electrolyte_db.get_base("thermo"))
    # Defining phase equilibria
    thermo_base.data.update(
        {
            "phases_in_equilibrium": [("Vap", "Liq")],
            "phase_equilibrium_state": {("Vap", "Liq"): SmoothVLE},
            "bubble_dew_method": IdealBubbleDew,
        }
    )
    thermo_base.data["phases"]["Vap"] = {"type": VaporPhase, "equation_of_state": Ideal}
    thermo_base.data["phases"]["Liq"]["type"] = AqueousPhase
    for comp in electrolyte_db.get_components(component_names):
        thermo_base.add(comp)
    result = thermo_base.idaes_config
    for name, component in result["components"].items():
        if name == "H2CO3":
            component["valid_phase_types"] = PT.aqueousPhase
        elif "valid_phase_types" in component.keys():
            del component["valid_phase_types"]
    assert compare_configs("thermo", result, reference=test_data.recarbonation_thermo_config) is not False
    return result


@pytest.fixture
def reaction_config():
    reaction_base = next(electrolyte_db.get_base("water_reaction"))
    # reaction_names = ("H2O_Kw", "CO2_to_H2CO3", "H2CO3_Ka1", "H2CO3_Ka2")
    reactions = list(electrolyte_db.get_reactions(component_names))
    print(f"Got reactions: {', '.join([r.name for r in reactions])}")
    for react in reactions:
        reaction_base.add(react)
    result = reaction_base.idaes_config
    assert compare_configs("reaction", result, reference=test_data.recarbonation_reaction_config) is not False
    return result


def compare_configs(name, conf, reference=None):
    print(f"Diffing {name} generated and reference")
    result = dict_diff(conf, reference)
    if result:
        print(f"They are different\nHow:")
        for difference in result:
            print(difference)
        return False  # different
    return True


def dict_diff(d1, d2, result=[], pfx=""):
    if isinstance(d1, list) and isinstance(d2, list):
        if len(d1) != len(d2):
            result.append(f"Array length at {pfx} first({len(d1)}) != second({len(d2)})")
        else:
            pass  # good enough
    elif not isinstance(d1, dict) or not isinstance(d2, dict):
        if type(d1) == type(d2):
            same = None
            try:
                same = d1 == d2
            except:   # cannot compare them
                pass  # good enough
            if same is False:
                result.append(f"value at {pfx} first != second")
        else:
            result.append(f"type of value at {pfx} first({type(d1)} != second({type(d2)}")
    else:
        if set(d1.keys()) != set(d2.keys()):
            for k in d1:
                if k not in d2:
                    result.append(f"{pfx}{k} in first, not in second")
            for k in d2:
                if k not in d1:
                    result.append(f"{pfx}{k} in first, not in second")
        for k in d1:
            if k in d2:
                pfx += f"{k}."
                dict_diff(d1[k], d2[k], result=result, pfx=pfx)
    return result


@pytest.fixture
def solver():
    return get_solver()


@pytest.fixture
def equilibrium_config(thermo_config, reaction_config):
    # Create a pyomo model object
    model = ConcreteModel()
    model.fs = FlowsheetBlock(default={"dynamic": False})

    print(f"thermo_config components: {thermo_config['components'].keys()}")

    model.fs.thermo_params = GenericParameterBlock(default=thermo_config)

    model.fs.rxn_params = GenericReactionParameterBlock(
        default={"property_package": model.fs.thermo_params, **reaction_config}
    )

    model.fs.unit = EquilibriumReactor(
        default={
            "property_package": model.fs.thermo_params,
            "reaction_package": model.fs.rxn_params,
            "has_rate_reactions": False,
            "has_equilibrium_reactions": True,
            "has_heat_transfer": False,
            "has_heat_of_reaction": False,
            "has_pressure_change": False,
        }
    )

    model.fs.unit.inlet.flow_mol.fix(10)
    model.fs.unit.inlet.pressure.fix(101325.0)
    model.fs.unit.inlet.temperature.fix(300.0)
    model.fs.unit.inlet.mole_frac_comp[0, "H +"].fix(0.0)
    model.fs.unit.inlet.mole_frac_comp[0, "OH -"].fix(0.0)
    model.fs.unit.inlet.mole_frac_comp[0, "H2CO3"].fix(0.0)
    model.fs.unit.inlet.mole_frac_comp[0, "HCO3 -"].fix(0.0)
    model.fs.unit.inlet.mole_frac_comp[0, "CO3 2-"].fix(0.0)
    model.fs.unit.inlet.mole_frac_comp[0, "CO2"].fix(0.0005)
    model.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix(1 - 0.0005)

    print(
        f"thermo_params components: {model.fs.thermo_params.component_list.value_list}"
    )
    return model


@pytest.mark.skip
@pytest.mark.unit
def test_build_model_equilibrium(equilibrium_config):
    model = equilibrium_config

    assert hasattr(model.fs.thermo_params, "component_list")
    assert len(model.fs.thermo_params.component_list) == 7
    assert "H2O" in model.fs.thermo_params.component_list
    assert isinstance(model.fs.thermo_params.H2O, Solvent)
    assert "H +" in model.fs.thermo_params.component_list
    assert isinstance(model.fs.thermo_params.component("H +"), Cation)
    assert "OH -" in model.fs.thermo_params.component_list
    assert isinstance(model.fs.thermo_params.component("OH -"), Anion)

    assert hasattr(model.fs.thermo_params, "phase_list")
    assert len(model.fs.thermo_params.phase_list) == 2
    assert isinstance(model.fs.thermo_params.Liq, AqueousPhase)
    assert isinstance(model.fs.thermo_params.Vap, VaporPhase)


@pytest.mark.skip
@pytest.mark.unit
def test_units_equilibrium(equilibrium_config):
    model = equilibrium_config
    assert_units_consistent(model)


@pytest.mark.skip
@pytest.mark.unit
def test_dof_equilibrium(equilibrium_config):
    model = equilibrium_config
    assert degrees_of_freedom(model) == 0


@pytest.mark.skip
@pytest.mark.component
def test_scaling_equilibrium(equilibrium_config):
    model = equilibrium_config
    iscale.calculate_scaling_factors(model.fs.unit)

    assert isinstance(model.fs.unit.control_volume.scaling_factor, Suffix)

    assert isinstance(
        model.fs.unit.control_volume.properties_out[0.0].scaling_factor, Suffix
    )

    assert isinstance(
        model.fs.unit.control_volume.properties_in[0.0].scaling_factor, Suffix
    )

    # When using equilibrium reactions, there are another set of scaling factors calculated
    assert isinstance(
        model.fs.unit.control_volume.reactions[0.0].scaling_factor, Suffix
    )


@pytest.mark.skip
@pytest.mark.component
def test_initialize_solver(equilibrium_config, solver):
    model = equilibrium_config
    solver.options["bound_push"] = 1e-10
    solver.options["mu_init"] = 1e-6
    model.fs.unit.initialize(optarg=solver.options)
    assert degrees_of_freedom(model) == 0


@pytest.mark.skip
@pytest.mark.component
def test_solve_equilibrium(equilibrium_config, solver):
    model = equilibrium_config
    solver.options["max_iter"] = 100
    results = solver.solve(model)
    print(results.solver.termination_condition)
    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok


@pytest.mark.skip
@pytest.mark.component
def test_solution_equilibrium(equilibrium_config):
    model = equilibrium_config

    assert pytest.approx(300, rel=1e-5) == value(model.fs.unit.outlet.temperature[0])
    assert pytest.approx(10, rel=1e-5) == value(model.fs.unit.outlet.flow_mol[0])
    assert pytest.approx(101325, rel=1e-5) == value(model.fs.unit.outlet.pressure[0])

    total_molar_density = (
        value(model.fs.unit.control_volume.properties_out[0.0].dens_mol_phase["Liq"])
        / 1000
    )
    assert pytest.approx(55.1847856, rel=1e-5) == total_molar_density
    pH = -value(
        log10(model.fs.unit.outlet.mole_frac_comp[0, "H +"] * total_molar_density)
    )
    pOH = -value(
        log10(model.fs.unit.outlet.mole_frac_comp[0, "OH -"] * total_molar_density)
    )
    assert pytest.approx(5.33698164, rel=1e-5) == pH
    assert pytest.approx(8.598909, rel=1e-5) == pOH

    CO2_sorbed = value(
        model.fs.unit.control_volume.properties_out[0.0].conc_mol_phase_comp[
            ("Liq", "CO2")
        ]
    )
    assert pytest.approx(27.5409958, rel=1e-5) == CO2_sorbed
