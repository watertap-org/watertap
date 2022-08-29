###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# 'https://github.com/watertap-org/watertap/'
#
###############################################################################
import pytest
import pyomo.environ as pyo
from pyomo.environ import (
    ConcreteModel,
    check_optimal_termination,
    value,
)
from pyomo.network import Port
from pyomo.util.check_units import assert_units_consistent

from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    MomentumBalanceType,
)
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    badly_scaled_var_generator,
)
from idaes.core import UnitModelCostingBlock

from watertap.core.util.scaling import variable_sens_generator
from watertap.property_models.ion_DSPMDE_prop_pack import (
    DSPMDEParameterBlock,
)
from watertap.unit_models.gac import (
    GAC,
    FilmTransferCoefficientType,
    SurfaceDiffusionCoefficientType,
)
from watertap.unit_models.electrodialysis_0D import Electrodialysis0D
from watertap.costing import WaterTAPCosting

__author__ = "Hunter Barber"

solver = get_solver()


# -----------------------------------------------------------------------------
# pull gac_frame_simplified from test_gac
class TestGACSimplified:
    @pytest.fixture(scope="class")
    def gac_frame_simplified(self):
        m_gac = ConcreteModel()
        m_gac.fs = FlowsheetBlock(default={"dynamic": False})

        m_gac.fs.properties = DSPMDEParameterBlock(
            default={"solute_list": ["DCE"], "mw_data": {"H2O": 18e-3, "DCE": 98.96e-3}}
        )

        m_gac.fs.unit = GAC(
            default={
                "property_package": m_gac.fs.properties,
                "film_transfer_coefficient_type": "fixed",
                "surface_diffusion_coefficient_type": "fixed",
            }
        )

        # feed specifications
        m_gac.fs.unit.process_flow.properties_in[0].pressure.fix(
            101325
        )  # feed pressure [Pa]
        m_gac.fs.unit.process_flow.properties_in[0].temperature.fix(
            273.15 + 25
        )  # feed temperature [K]
        m_gac.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
            "Liq", "H2O"
        ].fix(55555.55426666667)
        m_gac.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
            "Liq", "DCE"
        ].fix(0.0002344381568310428)

        # trial problem from Hand, 1984 for removal of trace DCE
        m_gac.fs.unit.conc_ratio_replace.fix(0.50)
        m_gac.fs.unit.freund_k.fix(37.9e-6 * (1e6**0.8316))
        m_gac.fs.unit.freund_ninv.fix(0.8316)
        m_gac.fs.unit.ebct.fix(300)  # seconds
        m_gac.fs.unit.bed_voidage.fix(0.449)
        m_gac.fs.unit.bed_length.fix(6)  # assumed
        m_gac.fs.unit.particle_porosity.fix(0.5)
        m_gac.fs.unit.particle_dens_app.fix(722)
        m_gac.fs.unit.particle_dia.fix(0.00106)
        m_gac.fs.unit.kf.fix(3.29e-5)
        m_gac.fs.unit.ds.fix(1.77e-13)
        m_gac.fs.unit.a0.fix(3.68421)
        m_gac.fs.unit.a1.fix(13.1579)
        m_gac.fs.unit.b0.fix(0.784576)
        m_gac.fs.unit.b1.fix(0.239663)
        m_gac.fs.unit.b2.fix(0.484422)
        m_gac.fs.unit.b3.fix(0.003206)
        m_gac.fs.unit.b4.fix(0.134987)

        return m_gac

    @pytest.mark.unit
    def test_standard_gac(self, gac_frame_simplified):
        # check specification for flowsheet in frame are solvable
        m_gac = gac_frame_simplified

        m_gac.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e-4, index=("Liq", "H2O")
        )
        m_gac.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e4, index=("Liq", "DCE")
        )
        calculate_scaling_factors(m_gac)
        initialization_tester(gac_frame_simplified)

        # Check var scaling before solve
        badly_scaled_var_lst = list(
            badly_scaled_var_generator(m_gac, large=1e2, small=1e-2, zero=1e-8)
        )
        assert badly_scaled_var_lst == []

        results = solver.solve(m_gac)

        # Check for optimal solution
        assert check_optimal_termination(results)

        # Check var scaling after solve
        badly_scaled_var_lst = list(
            badly_scaled_var_generator(m_gac, large=1e2, small=1e-2, zero=1e-8)
        )
        assert badly_scaled_var_lst == []

    @pytest.mark.unit
    def test_variable_sens_generator_gac(self, gac_frame_simplified):
        m_gac = gac_frame_simplified

        # test variable sens generator with gac model
        sens_var_lst = list(variable_sens_generator(m_gac, zero=1e-8))
        for i in sens_var_lst:
            print(i)

        assert sens_var_lst == []


# -----------------------------------------------------------------------------
# pull electrodialysis_cell_frame from test_electrodialysis_0D.py
class TestElectrodialysisVoltageConst:
    @pytest.fixture(scope="class")
    def electrodialysis_cell_frame(self):
        m_ed = ConcreteModel()

        m_ed.fs = FlowsheetBlock(default={"dynamic": False})
        ion_dict = {
            "solute_list": ["Na_+", "Cl_-"],
            "mw_data": {"H2O": 18e-3, "Na_+": 23e-3, "Cl_-": 35.5e-3},
            "electrical_mobility_data": {"Na_+": 5.19e-8, "Cl_-": 7.92e-8},
            "charge": {"Na_+": 1, "Cl_-": -1},
        }

        m_ed.fs.properties = DSPMDEParameterBlock(default=ion_dict)
        m_ed.fs.unit = Electrodialysis0D(
            default={
                "property_package": m_ed.fs.properties,
                "operation_mode": "Constant_Voltage",
            }
        )

        # Specify a system
        # Note: Testing scenarios in this file are primarily in accord with an experimental
        # setup reported by Campione et al. in Desalination 465 (2019): 79-93.
        # set the operational parameters
        m_ed.fs.unit.water_trans_number_membrane["cem"].fix(5.8)
        m_ed.fs.unit.water_trans_number_membrane["aem"].fix(4.3)
        m_ed.fs.unit.water_permeability_membrane["cem"].fix(2.16e-14)
        m_ed.fs.unit.water_permeability_membrane["aem"].fix(1.75e-14)
        m_ed.fs.unit.voltage.fix(0.5)
        m_ed.fs.unit.electrodes_resistance.fix(0)
        m_ed.fs.unit.cell_pair_num.fix(10)
        m_ed.fs.unit.current_utilization.fix(1)
        m_ed.fs.unit.spacer_thickness.fix(2.7e-4)
        m_ed.fs.unit.membrane_areal_resistance["cem"].fix(1.89e-4)
        m_ed.fs.unit.membrane_areal_resistance["aem"].fix(1.77e-4)
        m_ed.fs.unit.cell_width.fix(0.1)
        m_ed.fs.unit.cell_length.fix(0.79)
        m_ed.fs.unit.membrane_thickness["aem"].fix(1.3e-4)
        m_ed.fs.unit.membrane_thickness["cem"].fix(1.3e-4)
        m_ed.fs.unit.solute_diffusivity_membrane["cem", "Na_+"].fix(1.8e-10)
        m_ed.fs.unit.solute_diffusivity_membrane["aem", "Na_+"].fix(1.25e-10)
        m_ed.fs.unit.solute_diffusivity_membrane["cem", "Cl_-"].fix(1.8e-10)
        m_ed.fs.unit.solute_diffusivity_membrane["aem", "Cl_-"].fix(1.25e-10)
        m_ed.fs.unit.ion_trans_number_membrane["cem", "Na_+"].fix(1)
        m_ed.fs.unit.ion_trans_number_membrane["aem", "Na_+"].fix(0)
        m_ed.fs.unit.ion_trans_number_membrane["cem", "Cl_-"].fix(0)
        m_ed.fs.unit.ion_trans_number_membrane["aem", "Cl_-"].fix(1)

        # set the inlet stream
        m_ed.fs.unit.inlet_diluate.pressure.fix(101325)
        m_ed.fs.unit.inlet_diluate.temperature.fix(298.15)
        m_ed.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
        m_ed.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-4)
        m_ed.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-4)
        m_ed.fs.unit.inlet_concentrate.pressure.fix(101325)
        m_ed.fs.unit.inlet_concentrate.temperature.fix(298.15)
        m_ed.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
        m_ed.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(
            7.38e-4
        )
        m_ed.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(
            7.38e-4
        )

        return m_ed

    @pytest.mark.component
    def test_initialization_scaling(self, electrodialysis_cell_frame):
        m_ed = electrodialysis_cell_frame

        assert degrees_of_freedom(m_ed) == 0

        # set default scaling for state vars
        m_ed.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "H2O")
        )
        m_ed.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e4, index=("Liq", "Na_+")
        )
        m_ed.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e4, index=("Liq", "Cl_-")
        )

        calculate_scaling_factors(m_ed.fs)
        initialization_tester(m_ed)

        badly_scaled_var_values = {
            var.name: val for (var, val) in badly_scaled_var_generator(m_ed)
        }
        assert not badly_scaled_var_values

        results = solver.solve(m_ed)
        check_optimal_termination(results)

        badly_scaled_var_values = {
            var.name: val for (var, val) in badly_scaled_var_generator(m_ed)
        }
        assert not badly_scaled_var_values

    @pytest.mark.unit
    def test_variable_sens_generator_electrodialysis_0D(
        self, electrodialysis_cell_frame
    ):
        m_ed = electrodialysis_cell_frame

        # test variable sens generator with electrodialysis 0D model
        sens_var_lst = list(variable_sens_generator(m_ed))
        for i in sens_var_lst:
            print(i)

        assert sens_var_lst == []
