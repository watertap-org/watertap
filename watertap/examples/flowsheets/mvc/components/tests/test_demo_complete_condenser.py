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
from watertap.examples.flowsheets.mvc.components.demo_complete_condenser import main as main_complete_condenser
from idaes.core.util import get_solver

solver = get_solver()

# -----------------------------------------------------------------------------
@pytest.mark.component
def test_heat_exchanger(capsys):
    main_complete_condenser()
    captured = capsys.readouterr()

    assert captured.out == \
"""
====================================================================================
Unit : fs.unit                                                             Time: 0.0
------------------------------------------------------------------------------------
    Unit Performance

    Variables: 

    Key       : Value       : Fixed : Bounds
    Heat duty : -4.8715e+06 : False : (None, None)

------------------------------------------------------------------------------------
    Stream Table
                                           Inlet     Outlet  
    flow_mass_phase_comp ('Liq', 'H2O') 1.0000e-08     1.0000
    flow_mass_phase_comp ('Vap', 'H2O')     1.0000 1.0000e-10
    temperature                             400.00     340.00
    pressure                                50000.     50000.
====================================================================================
"""