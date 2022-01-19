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
from watertap.examples.flowsheets.mvc.components.demo_heat_exchanger import main as main_heat_exchanger
from idaes.core.util import get_solver

solver = get_solver()

# -----------------------------------------------------------------------------
@pytest.mark.component
def test_heat_exchanger(capsys):
    main_heat_exchanger()
    captured = capsys.readouterr()

    assert captured.out == \
"""
====================================================================================
Unit : fs.unit                                                             Time: 0.0
------------------------------------------------------------------------------------
    Unit Performance

    Variables: 

    Key            : Value  : Fixed : Bounds
           HX Area : 5.0000 :  True : (0, None)
    HX Coefficient : 1000.0 :  True : (0, None)
         Heat Duty : 89050. : False : (None, None)

    Expressions: 

    Key             : Value
    Delta T Driving : 17.810
         Delta T In : 9.2239
        Delta T Out : 30.689

------------------------------------------------------------------------------------
    Stream Table
                                         Hot Inlet  Hot Outlet  Cold Inlet  Cold Outlet
    flow_mass_phase_comp ('Liq', 'H2O')     1.0000      1.0000     0.50000     0.50000 
    flow_mass_phase_comp ('Liq', 'TDS')   0.010000    0.010000    0.010000    0.010000 
    temperature                             350.00      328.69      298.00      340.78 
    pressure                            2.0000e+05  2.0000e+05  2.0000e+05  2.0000e+05 
====================================================================================
"""