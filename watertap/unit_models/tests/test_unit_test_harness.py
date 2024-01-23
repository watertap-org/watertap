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
        self.prop_pack = props.NaClParameterBlock
        self.param_args = {}
        # build unit
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = props.NaClParameterBlock()
        m.fs.unit = ReverseOsmosis0D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
        )

        # fully specify system
        feed_flow_mass = 1
        feed_mass_frac_NaCl = 0.035
        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
        feed_pressure = 50e5
        feed_temperature = 273.15 + 25
        membrane_pressure_drop = 3e5
        length = 20
        membrane_area = 50
        A = 4.2e-12
        B = 3.5e-8
        pressure_atmospheric = 101325

        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )
        m.fs.unit.inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.deltaP.fix(-membrane_pressure_drop)
        m.fs.unit.area.fix(membrane_area)
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)

        m.fs.unit.feed_side.channel_height.fix(0.002)
        m.fs.unit.feed_side.spacer_porosity.fix(0.75)
        m.fs.unit.length.fix(length)

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
            "number_variables": 121,
            "number_total_constraints": 109,
            "number_unused_variables": 0,
        }
        self.unit_solution = {
            # "feed_side.properties_out[0].conc_mass_phase_comp[Liq, NaCl]": 1.5,
            ("flow_mass_phase_comp", ("Liq", "NaCl")): 0.035,
            # "flux_mass_phase_comp_avg[0, Liq, H2O]": 4.722e-3,
            # "flux_mass_phase_comp_avg[0, Liq, NaCl]": 1.576e-6,
            # "mixed_permeate[0].flow_mass_phase_comp[Liq, H2O]": 0.2361,
            # "mixed_permeate[0].flow_mass_phase_comp[Liq, NaCl]": 7.879e-5,
            # "feed_side.mixed_permeate[0].flow_mass_phase_comp[Liq, H2O]": 0.95,
            # "feed_side.mixed_permeate[0].flow_mass_phase_comp[Liq, NaCl]": 0.05,
            # ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.05,
            # ("dens_mass_phase", "Liq"): 1032.8,
            # ("flow_vol_phase", "Liq"): 9.682e-4,
            # ("conc_mass_phase_comp", ("Liq", "H2O")): 981.1,
            # ("conc_mass_phase_comp", ("Liq", "NaCl")): 51.64,
            # ("flow_mol_phase_comp", ("Liq", "H2O")): 52.73,
            # ("flow_mol_phase_comp", ("Liq", "NaCl")): 0.8556,
            # ("mole_frac_phase_comp", ("Liq", "H2O")): 0.9840,
            # ("mole_frac_phase_comp", ("Liq", "NaCl")): 1.597e-2,
            # ("molality_phase_comp", ("Liq", "NaCl")): 0.9006,
            # ("diffus_phase_comp", ("Liq", "NaCl")): 1.471e-9,
            # ("visc_d_phase", "Liq"): 1.0875e-3,
            # ("osm_coeff", None): 0.9347,
            # ("pressure_osm_phase", "Liq"): 4.174e6,
            # ("enth_mass_phase", "Liq"): 1.093e5,
        }


# @pytest.mark.component
# class TestUnitSolution(UnitRegressionTest):
#     def configure(self):
#         self.prop_pack = props.NaClParameterBlock
#         self.param_args = {}
#
#         self.solver = "ipopt"
#         self.optarg = {"nlp_scaling_method": "user-scaling"}
#
#         self.scaling_args = {
#             ("flow_mass_phase_comp", ("Liq", "H2O")): 1,
#             ("flow_mass_phase_comp", ("Liq", "NaCl")): 1e2,
#         }
#         self.state_args = {
#             ("flow_mass_phase_comp", ("Liq", "H2O")): 0.95,
#             ("flow_mass_phase_comp", ("Liq", "NaCl")): 0.05,
#             ("temperature", None): 273.15 + 25,
#             ("pressure", None): 50e5,
#         }
#         self.unit_solution = {
#             "flux_mass_phase_comp_avg[0, Liq, H2O]": 4.722e-3,
#             "flux_mass_phase_comp_avg[0, Liq, NaCl]": 1.576e-6,
#             "mixed_permeate[0].flow_mass_phase_comp[Liq, H2O]": 0.2361,
#             "mixed_permeate[0].flow_mass_phase_comp[Liq, NaCl]": 7.879e-5,
#             "feed_side.mixed_permeate[0].flow_mass_phase_comp[Liq, H2O]": 0.95,
#             "feed_side.mixed_permeate[0].flow_mass_phase_comp[Liq, NaCl]": 0.05,
# ("mass_frac_phase_comp", ("Liq", "NaCl")): 0.05,
# ("dens_mass_phase", "Liq"): 1032.8,
# ("flow_vol_phase", "Liq"): 9.682e-4,
# ("conc_mass_phase_comp", ("Liq", "H2O")): 981.1,
# ("conc_mass_phase_comp", ("Liq", "NaCl")): 51.64,
# ("flow_mol_phase_comp", ("Liq", "H2O")): 52.73,
# ("flow_mol_phase_comp", ("Liq", "NaCl")): 0.8556,
# ("mole_frac_phase_comp", ("Liq", "H2O")): 0.9840,
# ("mole_frac_phase_comp", ("Liq", "NaCl")): 1.597e-2,
# ("molality_phase_comp", ("Liq", "NaCl")): 0.9006,
# ("diffus_phase_comp", ("Liq", "NaCl")): 1.471e-9,
# ("visc_d_phase", "Liq"): 1.0875e-3,
# ("osm_coeff", None): 0.9347,
# ("pressure_osm_phase", "Liq"): 4.174e6,
# ("enth_mass_phase", "Liq"): 1.093e5,
# }

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
