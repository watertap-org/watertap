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
from io import StringIO
from watertap.examples.case_studies.santa_barbara.santa_barbara import main
from idaes.core.util import get_solver

# -----------------------------------------------------------------------------
@pytest.mark.component
def test_santa_barbara():
    m = main()

    report_io = StringIO()
    m.fs.feed.report(ostream=report_io)
    output = \
        """
====================================================================================
Unit : fs.feed                                                             Time: 0.0
------------------------------------------------------------------------------------
    Stream Table
                            Outlet 
    Mass Concentration H2O   988.47
    Mass Concentration tds   35.000
    Mass Concentration tss 0.030000
    Volumetric Flowrate     0.30920
====================================================================================
"""
    assert output == report_io.getvalue()

    report_io = StringIO()
    m.fs.desalination.RO.report(ostream=report_io)
    output = \
        """
====================================================================================
Unit : fs.desalination.RO                                                  Time: 0.0
------------------------------------------------------------------------------------
    Unit Performance

    Variables: 

    Key                        : Value   : Fixed : Bounds
                 Membrane Area :  13914. :  True : (0.1, 1000.0)
    Solvent Mass Recovery Rate : 0.43681 : False : (0.01, 0.999999)
      Volumetric Recovery Rate : 0.43293 : False : (0.01, 0.999999)

------------------------------------------------------------------------------------
    Stream Table
                                         Feed Inlet  Feed Outlet  Permeate Outlet
    flow_mass_phase_comp ('Liq', 'H2O')      305.57      172.09         133.48   
    flow_mass_phase_comp ('Liq', 'TDS')      10.822      10.792       0.029557   
    temperature                              298.00      298.02         298.02   
    pressure                             7.0000e+06  6.7759e+06     1.0132e+05   
====================================================================================
"""
    assert output == report_io.getvalue()