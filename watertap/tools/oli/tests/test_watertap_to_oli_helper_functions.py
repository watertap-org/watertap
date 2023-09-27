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
from watertap.tools.oli.watertap_to_oli_helper_functions import watertap_to_oli, oli_reverse_lookup, names_db

@pytest.mark.unit
def test_watertap_to_oli():
    test_conditions = []
    # neutral, no brackets
    test_conditions.append('NaOH')
    # neutral, bracket
    test_conditions.append('B[OH]3')
    # charged, bracket
    test_conditions.append('B[OH]4_-')
    # cation, z = 1
    test_conditions.append('K_+')
    # anion, z = -1
    test_conditions.append('Cl_-')
    # cation, z > 1
    test_conditions.append('Mg_2+')
    # anion Z < -1
    test_conditions.append('HCO3_2-')

    for test_name in test_conditions:
        watertap_to_oli(test_name)
        
@pytest.mark.unit
def test_reverse_lookup():
    oli_reverse_lookup("NAION", names_db)