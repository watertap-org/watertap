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

__author__ = "Alexander V. Dudchenko"


from watertap.core.util import report

from pyomo.environ import ConcreteModel, Var, units as pyunits
from io import StringIO


def test_report_table():
    """Test that the report table can be built without error."""
    # Create a dummy model with a variable

    m = ConcreteModel()
    m.x = Var(["t1", "t2", "t3"], initialize=1.0, units=pyunits.m)
    os = StringIO()
    # Build the report table
    # should look like this:
    # ------------------------------------------------------------------------------------
    #     test_var state

    #     test_value:
    #     Key : Value  : Units : Fixed : Bounds
    #      t1 : 1.0000 :     m : False : (None, None)
    #      t2 : 1.0000 :     m : False : (None, None)
    #      t3 : 1.0000 :     m : False : (None, None)

    # ------------------------------------------------------------------------------------
    report.build_report_table("test_var", {"test_value": m.x}, os)
    result = "\n------------------------------------------------------------------------------------\n    test_var state\n\n    test_value: \n    Key : Value  : Units : Fixed : Bounds\n     t1 : 1.0000 :     m : False : (None, None)\n     t2 : 1.0000 :     m : False : (None, None)\n     t3 : 1.0000 :     m : False : (None, None)\n\n------------------------------------------------------------------------------------\n"
    print(os.getvalue())
    assert os.getvalue() == result, "Report table does not match expected output"
    os = StringIO()
    # Build the report table
    report.build_report_table(
        "test_var", {"test_value": m.x}, os, use_default_units=True
    )
    # should look like this:
    # ------------------------------------------------------------------------------------
    #     test_var state

    #     test_value:
    #     Key : Value  : Units : Fixed : Bounds
    #      t1 : 1.0000 :     meter : False : (None, None)
    #      t2 : 1.0000 :     meter : False : (None, None)
    #      t3 : 1.0000 :     meter : False : (None, None)

    # ------------------------------------------------------------------------------------
    result = "\n------------------------------------------------------------------------------------\n    test_var state\n\n    test_value: \n    Key : Value  : Units : Fixed : Bounds\n     t1 : 1.0000 : meter : False : (None, None)\n     t2 : 1.0000 : meter : False : (None, None)\n     t3 : 1.0000 : meter : False : (None, None)\n\n------------------------------------------------------------------------------------\n"
    print(os.getvalue())
    assert os.getvalue() == result, "Report table does not match expected output"
