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
                           value,
                           Param,
                           Var,
                           Expression,
                           Constraint,
                           assert_optimal_termination)
from pyomo.network import Port
from idaes.core import (FlowsheetBlock,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        ControlVolume0DBlock)
from watertap.unit_models.reverse_osmosis_0D import (ReverseOsmosis0D,
                                                     ConcentrationPolarizationType,
                                                     MassTransferCoefficient,
                                                     PressureChangeType)
import watertap.property_models.NaCl_prop_pack as props
from watertap.unit_models.tests.unit_test_harness import UnitTestHarness

from idaes.core.util import get_solver
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              number_unused_variables)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (calculate_scaling_factors,
                                     unscaled_variables_generator,
                                     unscaled_constraints_generator,
                                     badly_scaled_var_generator)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver(options={'bound_push': 1e-8})

# -----------------------------------------------------------------------------
class TestReverseOsmosis0D_default(UnitTestHarness):
    def configure(self):
        # build unit
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.properties = props.NaClParameterBlock()
        m.fs.unit = ReverseOsmosis0D(default={
            "property_package": m.fs.properties,
            "has_pressure_change": True,
            "pressure_change_type": PressureChangeType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.calculated,
            "concentration_polarization_type": ConcentrationPolarizationType.calculated,
        })

        # specify unit
        m.fs.unit.feed_side.properties_in[0].flow_mass_phase_comp['Liq', 'NaCl'].fix(0.035)
        m.fs.unit.feed_side.properties_in[0].flow_mass_phase_comp['Liq', 'H2O'].fix(0.965)
        m.fs.unit.feed_side.properties_in[0].pressure.fix(50e5)
        m.fs.unit.feed_side.properties_in[0].temperature.fix(298.15)
        m.fs.unit.area.fix(50)
        m.fs.unit.A_comp.fix(4.2e-12)
        m.fs.unit.B_comp.fix(3.5e-8)
        m.fs.unit.permeate.pressure[0].fix(101325)
        m.fs.unit.channel_height.fix(0.002)
        m.fs.unit.spacer_porosity.fix(0.85)
        m.fs.unit.length.fix(20)

        # scale unit
        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'NaCl'))
        calculate_scaling_factors(m)

        self.unit_model = m.fs.unit  # TODO: change to unit_model_block
        self.unit_statistics = {'number_variables': 128,
                                'number_total_constraints': 117,
                                'number_unused_variables': 0}
        self.unit_solution = {}  # TODO: change to report output
