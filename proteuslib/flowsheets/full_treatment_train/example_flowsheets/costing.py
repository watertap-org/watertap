###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################

from pyomo.environ import (
    Block, Constraint, Expression, Var, Param, Reals, NonNegativeReals, units as pyunits)
from idaes.core.util.exceptions import ConfigurationError
import proteuslib.flowsheets.RO_with_energy_recovery.financials as financials
import proteuslib.flowsheets.full_treatment_train.example_flowsheets.financials as financials
from proteuslib.unit_models.pump_isothermal import Pump, PumpData
from proteuslib.flowsheets.full_treatment_train.example_flowsheets.flowsheet_limited import *
def build_costing_flowsheet(m, module=None):
    '''
    Add costing to a given flowsheet. This function will
        1) add a block which contains cost parameters,
        2) call the get_costing method for each unit model (note: unit model must have a get_costing method
        to be detected), and
        3) call get_system_costing which will tally up all capex and opex for each process
    m : model
    module : financials module
    '''

    # Add the costing parameter block
    module.add_costing_param_block(m.fs)

    # call get_costing for each unit model
    get_costing_sweep(m.fs, module=financials)

    # call get_system_costing for whole flowsheet
    # module.get_system_costing(m.fs)

def get_costing_sweep(self, module= None):
    for b_unit in self.component_objects(Block, descend_into=True):
        # print(b_unit)
        if hasattr(b_unit, 'get_costing') and callable(b_unit.get_costing):
            name = getattr(b_unit, 'local_name')
            # if getattr(b_unit, '__class__') == 'idaes.core.process_block._ScalarPump':
            if isinstance(b_unit, PumpData):
                print(f"We got ourselves a pump called {name}!")
            else:
                print(f"We got ourselves a {name}!")
                # b_unit.get_costing(module=module)

if __name__ == "__main__":
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    build_flowsheet_limited_NF(m, has_bypass=True, has_desal_feed=False, is_twostage=False,
                               NF_type='ZO', NF_base='ion',
                               RO_type='Sep', RO_base='TDS', RO_level='simple')
    build_costing_flowsheet(m, module=financials)