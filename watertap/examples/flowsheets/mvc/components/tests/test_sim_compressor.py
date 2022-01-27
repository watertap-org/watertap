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
from watertap.examples.flowsheets.mvc.components.sim_compressor import main as main_compressor
from idaes.core.util import get_solver

solver = get_solver()

# -----------------------------------------------------------------------------
@pytest.mark.component
def test_compressor(capsys):
    main_compressor()
    captured = capsys.readouterr()

    assert captured.out == \
"""
====================================================================================
Unit : fs.compressor                                                       Time: 0.0
------------------------------------------------------------------------------------
    Unit Performance

    Variables: 

    Key            : Value      : Fixed : Bounds
        Efficiency :    0.80000 :  True : (1e-08, 1)
    Pressure ratio :     2.0000 :  True : (1, 10)
              Work : 4.0616e+06 : False : (1, 5000000.0)

------------------------------------------------------------------------------------
    Stream Table
                                           Inlet     Outlet  
    flow_mass_phase_comp ('Liq', 'H2O') 1.0000e-08 1.0000e-08
    flow_mass_phase_comp ('Vap', 'H2O')     35.214     35.214
    temperature                             350.00     429.57
    pressure                                50000. 1.0000e+05
====================================================================================
"""