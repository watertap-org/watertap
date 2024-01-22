#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pytest
from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock
from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
import watertap.property_models.NaCl_prop_pack as props
from watertap.unit_models.tests.unit_test_harness import UnitTestHarness

from idaes.core.solvers import get_solver
from idaes.core.util.scaling import calculate_scaling_factors

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver(options={"bound_push": 1e-8})


# -----------------------------------------------------------------------------
class TestReverseOsmosis0D_default(UnitTestHarness):
    def configure(self):
        # build unit
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = props.NaClParameterBlock()
        m.fs.unit = ReverseOsmosis0D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            pressure_change_type=PressureChangeType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
        )

        # specify unit
        m.fs.unit.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
            0.035
        )
        m.fs.unit.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "H2O"].fix(
            0.965
        )
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
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
        )
        calculate_scaling_factors(m)

        self.unit_model_block = m.fs.unit
        self.unit_statistics = {
            "number_variables": 128,
            "number_total_constraints": 117,
            "number_unused_variables": 0,
        }


#         self.unit_report = \
# """
# ====================================================================================
# Unit : fs.unit                                                             Time: 0.0
# ------------------------------------------------------------------------------------
#     Unit Performance
#
#     Variables:
#
#     Key                                            : Value       : Fixed : Bounds
#                                      Membrane Area :      50.000 :  True : (0.1, 1000.0)
#                                    Membrane Length :      20.000 :  True : (0.1, 500.0)
#                                     Membrane Width :      2.5000 : False : (0.1, 500.0)
#                    NaCl Concentration @Inlet,Bulk  :      35.751 : False : (0.001, 2000.0)
#      NaCl Concentration @Inlet,Membrane-Interface  :      43.381 : False : (0.001, 2000.0)
#                   NaCl Concentration @Outlet,Bulk  :      46.102 : False : (0.001, 2000.0)
#     NaCl Concentration @Outlet,Membrane-Interface  :      50.832 : False : (0.001, 2000.0)
#                       NaCl Permeate Concentration  :     0.36827 : False : (0.001, 2000.0)
#         Osmotic Pressure @Inlet,Membrane-Interface :  3.4820e+06 : False : (500.0, 50000000.0)
#                      Osmotic Pressure @Outlet,Bulk :  3.7086e+06 : False : (500.0, 50000000.0)
#       Osmotic Pressure @Outlet,Membrane-Interface  :  4.1055e+06 : False : (500.0, 50000000.0)
#                                    Pressure Change : -1.7701e+05 : False : (None, None)
#                             Reynolds Number @Inlet :      473.82 : False : (10, 5000.0)
#                            Reynolds Number @Outlet :      362.00 : False : (10, 5000.0)
#                         Solvent Mass Recovery Rate :     0.22864 : False : (0.01, 0.999999)
#                                    Velocity @Inlet :     0.23035 : False : (0.01, 5)
#                                   Velocity @Outlet :     0.17821 : False : (0.01, 5)
#                         Volumetric Flowrate @Inlet :  0.00097899 : False : (1e-08, None)
#                        Volumetric Flowrate @Outlet :  0.00075741 : False : (1e-08, None)
#                           Volumetric Recovery Rate :     0.22653 : False : (0.01, 0.999999)
#
# ------------------------------------------------------------------------------------
#     Stream Table
#                                           Feed Inlet  Feed Outlet  Permeate Outlet
#     flow_mass_phase_comp ('Liq', 'H2O')      0.96500     0.74436        0.22064
#     flow_mass_phase_comp ('Liq', 'NaCl')    0.035000    0.034918     8.1671e-05
#     temperature                               298.15      298.15         298.15
#     pressure                              5.0000e+06  4.8230e+06     1.0132e+05
# ====================================================================================
# """
