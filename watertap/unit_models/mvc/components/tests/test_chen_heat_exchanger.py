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
import sys
import pytest
from io import StringIO

from pyomo.environ import ConcreteModel, assert_optimal_termination
from pyomo.util.check_units import assert_units_consistent
from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.generic_models.unit_models.heat_exchanger import (
    HeatExchanger,
    HeatExchangerFlowPattern,
)
import idaes.core.util.scaling as iscale

from watertap.unit_models.mvc.components.lmtd_chen_callback import (
    delta_temperature_chen_callback,
)
import watertap.property_models.seawater_prop_pack as props

solver = get_solver()

# -----------------------------------------------------------------------------
@pytest.mark.requires_idaes_solver
@pytest.mark.component
def test_heat_exchanger():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = props.SeawaterParameterBlock()
    m.fs.unit = HeatExchanger(
        default={
            "hot_side_name": "hot",
            "cold_side_name": "cold",
            "hot": {"property_package": m.fs.properties},
            "cold": {"property_package": m.fs.properties},
            "delta_temperature_callback": delta_temperature_chen_callback,
            "flow_pattern": HeatExchangerFlowPattern.countercurrent,
        }
    )

    # scaling
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "TDS")
    )
    iscale.set_scaling_factor(m.fs.unit.hot.heat, 1e-3)
    iscale.set_scaling_factor(m.fs.unit.cold.heat, 1e-3)
    iscale.set_scaling_factor(m.fs.unit.overall_heat_transfer_coefficient, 1e-3)
    iscale.set_scaling_factor(m.fs.unit.area, 1)
    iscale.calculate_scaling_factors(m)

    # ---specifications---
    # state variables
    m.fs.unit.hot_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(1)
    m.fs.unit.hot_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(0.01)
    m.fs.unit.hot_inlet.temperature[0].fix(350)
    m.fs.unit.hot_inlet.pressure[0].fix(2e5)

    m.fs.unit.cold_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(0.5)
    m.fs.unit.cold_inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(0.01)
    m.fs.unit.cold_inlet.temperature[0].fix(298)
    m.fs.unit.cold_inlet.pressure[0].fix(2e5)

    m.fs.unit.area.fix(5)
    m.fs.unit.overall_heat_transfer_coefficient.fix(1000)

    # solving
    assert_units_consistent(m)
    degrees_of_freedom(m)

    m.fs.unit.initialize()

    solver = get_solver()
    results = solver.solve(m, tee=False)
    assert_optimal_termination(results)

    report_io = StringIO()
    m.fs.unit.report(ostream=report_io)
    output = """
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
    assert output == report_io.getvalue()
