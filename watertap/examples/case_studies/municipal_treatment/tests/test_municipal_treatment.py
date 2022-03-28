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
from pyomo.environ import value
from watertap.examples.case_studies.municipal_treatment.municipal_treatment import main
from idaes.core.util import get_solver

# -----------------------------------------------------------------------------
@pytest.mark.component
def test_municipal_treatment():
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
    Mass Concentration H2O    999.36
    Mass Concentration tds   0.63000
    Mass Concentration toc 0.0040000
    Mass Concentration tss 0.0065250
    Volumetric Flowrate      0.92240
====================================================================================
"""
    assert output == report_io.getvalue()

    report_io = StringIO()
    m.fs.intake_pump.report(ostream=report_io)
    output = \
        """
====================================================================================
Unit : fs.intake_pump                                                      Time: 0.0
------------------------------------------------------------------------------------
    Unit Performance

    Variables: 

    Key                : Value  : Fixed : Bounds
    Electricity Demand : 93.200 :  True : (0, None)

------------------------------------------------------------------------------------
    Stream Table
                              Inlet    Outlet 
    Volumetric Flowrate      0.92240   0.92240
    Mass Concentration H2O    999.36    999.36
    Mass Concentration tds   0.63000   0.63000
    Mass Concentration tss 0.0065250 0.0065250
    Mass Concentration toc 0.0040000 0.0040000
====================================================================================
"""
    assert output == report_io.getvalue()

    report_io = StringIO()
    m.fs.gac.report(ostream=report_io)
    output = \
        """
====================================================================================
Unit : fs.gac                                                              Time: 0.0
------------------------------------------------------------------------------------
    Unit Performance

    Variables: 

    Key                     : Value     : Fixed : Bounds
    Activated Carbon Demand :    182.91 : False : (0, None)
         Electricity Demand :    28.873 : False : (0, None)
      Electricity Intensity : 0.0086958 : False : (None, None)
     Empty Bed Contact Time :   0.16667 :  True : (0, None)
       Solute Removal [tds] :    0.0000 :  True : (0, None)
       Solute Removal [toc] :    0.0000 :  True : (0, None)
       Solute Removal [tss] :   0.97000 :  True : (0, None)
             Water Recovery :   0.96000 :  True : (1e-08, 1.0000001)

------------------------------------------------------------------------------------
    Stream Table
                              Inlet     Treated   Byproduct
    Volumetric Flowrate       0.92230    0.88543  0.036869 
    Mass Concentration H2O     999.37     999.34    1000.0 
    Mass Concentration tds    0.63007    0.65630    0.0000 
    Mass Concentration tss 5.9867e-05 1.8708e-06 0.0014527 
    Mass Concentration toc  0.0010869  0.0011322    0.0000 
====================================================================================
"""
    assert output == report_io.getvalue()

    report_io = StringIO()
    m.fs.recharge_pump.report(ostream=report_io)
    output = \
        """
====================================================================================
Unit : fs.recharge_pump                                                    Time: 0.0
------------------------------------------------------------------------------------
    Unit Performance

    Variables: 

    Key                : Value  : Fixed : Bounds
    Electricity Demand : 186.40 :  True : (0, None)

------------------------------------------------------------------------------------
    Stream Table
                              Inlet     Outlet  
    Volumetric Flowrate       0.88491    0.88491
    Mass Concentration H2O     999.93     999.93
    Mass Concentration tds   0.065669   0.065669
    Mass Concentration tss 1.8719e-06 1.8719e-06
    Mass Concentration toc  0.0010710  0.0010710
====================================================================================
"""
    assert output == report_io.getvalue()
#
#     assert value(m.LCOW) == pytest.approx(VALUE, rel=1e-3) # TODO: add test for LCOW
