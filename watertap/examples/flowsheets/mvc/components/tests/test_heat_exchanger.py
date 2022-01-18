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
from pyomo.environ import (ConcreteModel,
                           Block,
                           Var,
                           Constraint,
                           TerminationCondition,
                           SolverStatus,
                           value,
                           SolverFactory,
                           Expression,
                           TransformationFactory,
                           units as pyunits)
from watertap.examples.flowsheets.mvc.components.heat_exchanger import main as main_heat_exchanger
from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import (solve_indexed_blocks,
                                            propagate_state)
from idaes.generic_models.unit_models import Mixer, Separator, Product, Feed
from idaes.generic_models.unit_models.mixer import MomentumMixingType
from pyomo.util.check_units import assert_units_consistent
from idaes.core.util.scaling import (unscaled_variables_generator,
                                     unscaled_constraints_generator)

import watertap.property_models.NaCl_prop_pack as props
from watertap.unit_models.reverse_osmosis_0D import ReverseOsmosis0D
from watertap.unit_models.pressure_exchanger import PressureExchanger
from watertap.unit_models.pump_isothermal import Pump
from watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery import (
build, set_operating_conditions, initialize_system, solve, optimize_set_up, optimize, display_system, display_state, display_design)


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