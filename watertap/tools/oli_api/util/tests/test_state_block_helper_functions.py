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
#
#################################################################################

import pytest
from pyomo.environ import units as pyunits
from watertap.tools.oli_api.util.state_block_helper_functions import (
    create_state_block,
    extract_state_vars,
)


@pytest.fixture
def chemistry_source():
    return {
        "temperature": 298.15,  # temperature in K
        "pressure": 101325,  # pressure in Pa
        "components": {  # concentrations in mg/L
            "Cl_-": 870,
            "Na_+": 739,
            "SO4_2-": 1011,
            "Mg_2+": 90,
            "Ca_2+": 258,
            "K_+": 9,
            "HCO3_-": 385,
        },
        "units": {
            "temperature": pyunits.K,
            "pressure": pyunits.Pa,
            "components": pyunits.mg / pyunits.L,
        },
    }


@pytest.fixture
def state_block(chemistry_source):
    return create_state_block(chemistry_source)


@pytest.mark.unit
def test_extract_state_vars(chemistry_source, state_block):
    conc_var = state_block.fs.stream[0].conc_mass_phase_comp
    units = chemistry_source["units"]
    extract_state_vars(state_block.fs.stream[0], conc_var, chemistry_source["units"])
