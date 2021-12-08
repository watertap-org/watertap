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
import pytest
from watertap.examples.edb.the_basics import run_the_basics_with_mockdb
from watertap.examples.edb.simple_acid import run_simple_acid_with_mockdb
from watertap.examples.edb.vapor_liquid_equilibrium import run_vap_liq_with_mockdb

# Import pyomo methods to check the system units
from pyomo.util.check_units import assert_units_consistent
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
)
from idaes.core import AqueousPhase, VaporPhase, SolidPhase
from idaes.core.components import Solvent, Solute, Cation, Anion

__author__ = "Austin Ladshaw"

@pytest.mark.component
def test_the_basics(edb):
    assert run_the_basics_with_mockdb(edb) == True

@pytest.mark.component
def test_simple_acid(edb):
    model = run_simple_acid_with_mockdb(edb)

    # Check model stats
    assert_units_consistent(model)
    assert degrees_of_freedom(model) == 8
    assert len(model.fs.unit.inlet) == 1
    assert len(model.fs.unit.inlet.flow_mol_phase_comp) == 6

    # Check for correct species and phases
    assert hasattr(model.fs.thermo_params, 'component_list')
    assert len(model.fs.thermo_params.component_list) == 6
    assert 'H2O' in model.fs.thermo_params.component_list
    assert isinstance(model.fs.thermo_params.H2O, Solvent)
    assert 'H_+' in model.fs.thermo_params.component_list
    assert isinstance(model.fs.thermo_params.component('H_+'), Cation)
    assert 'OH_-' in model.fs.thermo_params.component_list
    assert isinstance(model.fs.thermo_params.component('OH_-'), Anion)
    assert 'HCO3_-' in model.fs.thermo_params.component_list
    assert isinstance(model.fs.thermo_params.component('HCO3_-'), Anion)
    assert 'CO3_2-' in model.fs.thermo_params.component_list
    assert isinstance(model.fs.thermo_params.component('CO3_2-'), Anion)
    assert 'H2CO3' in model.fs.thermo_params.component_list
    assert isinstance(model.fs.thermo_params.H2CO3, Solute)

    assert hasattr(model.fs.thermo_params, 'phase_list')
    assert len(model.fs.thermo_params.phase_list) == 1
    assert isinstance(model.fs.thermo_params.Liq, AqueousPhase)


@pytest.mark.component
def test_vapor_liquid_equilibrium(edb):
    model = run_vap_liq_with_mockdb(edb)

    # Check model stats
    assert_units_consistent(model)
    assert degrees_of_freedom(model) == 11

    # Check for correct species and phases
    assert hasattr(model.fs.thermo_params, 'component_list')
    assert len(model.fs.thermo_params.component_list) == 7
    assert 'H2O' in model.fs.thermo_params.component_list
    assert isinstance(model.fs.thermo_params.H2O, Solvent)
    assert 'CO2' in model.fs.thermo_params.component_list
    assert isinstance(model.fs.thermo_params.CO2, Solute)
    assert 'H_+' in model.fs.thermo_params.component_list
    assert isinstance(model.fs.thermo_params.component('H_+'), Cation)
    assert 'OH_-' in model.fs.thermo_params.component_list
    assert isinstance(model.fs.thermo_params.component('OH_-'), Anion)
    assert 'HCO3_-' in model.fs.thermo_params.component_list
    assert isinstance(model.fs.thermo_params.component('HCO3_-'), Anion)
    assert 'CO3_2-' in model.fs.thermo_params.component_list
    assert isinstance(model.fs.thermo_params.component('CO3_2-'), Anion)
    assert 'H2CO3' in model.fs.thermo_params.component_list
    assert isinstance(model.fs.thermo_params.H2CO3, Solute)

    assert hasattr(model.fs.thermo_params, 'phase_list')
    assert len(model.fs.thermo_params.phase_list) == 2
    assert isinstance(model.fs.thermo_params.Liq, AqueousPhase)
    assert isinstance(model.fs.thermo_params.Vap, VaporPhase)
