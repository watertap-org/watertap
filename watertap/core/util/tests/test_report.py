from watertap.core.util import report

from pyomo.environ import ConcreteModel, Var, units as pyunits
import pytest

from io import StringIO


import pyomo.common.unittest as unittest


def test_report_table():
    """Test that the report table can be built without error."""
    # Create a dummy model with a variable
    from pyomo.environ import ConcreteModel, Var, units as pyunits

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
