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
"""
Tests for anaerobic digester example.

Verified against results from:

Rosen, C. and Jeppsson, U., 2006.
Aspects on ADM1 Implementation within the BSM2 Framework.
Department of Industrial Electrical Engineering and Automation, Lund University, Lund, Sweden, pp.1-35.

"""
import pytest
from pyomo.environ import (
    ConcreteModel,
    Suffix,
    TransformationFactory,
    Var,
)

from idaes.core import (
    FlowsheetBlock,
)
from idaes.core.util.scaling import (
    get_jacobian,
    jacobian_cond,
)
from watertap.core.solvers import get_solver

from watertap.unit_models.anaerobic_digester import AD, ADScaler
from watertap.property_models.unit_specific.anaerobic_digestion.adm1_properties import (
    ADM1ParameterBlock,
    ADM1PropertiesScaler,
)
from watertap.property_models.unit_specific.anaerobic_digestion.adm1_properties_vapor import (
    ADM1_vaporParameterBlock,
)
from watertap.property_models.unit_specific.anaerobic_digestion.adm1_reactions import (
    ADM1ReactionParameterBlock,
    ADM1ReactionScaler,
)

from watertap.unit_models.tests.unit_test_harness import UnitTestHarness
import idaes.core.util.scaling as iscale
from idaes.core.scaling.scaling_base import ScalerBase

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------


# TODO: Refine testing once iscale functionality has been deprecated/removed
def build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props = ADM1ParameterBlock()
    m.fs.props_vap = ADM1_vaporParameterBlock()
    m.fs.rxn_props = ADM1ReactionParameterBlock(property_package=m.fs.props)

    m.fs.unit = AD(
        liquid_property_package=m.fs.props,
        vapor_property_package=m.fs.props_vap,
        reaction_package=m.fs.rxn_props,
        has_heat_transfer=True,
        has_pressure_change=False,
    )

    # Set the operating conditions
    m.fs.unit.inlet.flow_vol.fix(170 / 24 / 3600)
    m.fs.unit.inlet.temperature.fix(308.15)
    m.fs.unit.inlet.pressure.fix(101325)

    m.fs.unit.inlet.conc_mass_comp[0, "S_su"].fix(0.01)
    m.fs.unit.inlet.conc_mass_comp[0, "S_aa"].fix(0.001)
    m.fs.unit.inlet.conc_mass_comp[0, "S_fa"].fix(0.001)
    m.fs.unit.inlet.conc_mass_comp[0, "S_va"].fix(0.001)
    m.fs.unit.inlet.conc_mass_comp[0, "S_bu"].fix(0.001)
    m.fs.unit.inlet.conc_mass_comp[0, "S_pro"].fix(0.001)
    m.fs.unit.inlet.conc_mass_comp[0, "S_ac"].fix(0.001)
    m.fs.unit.inlet.conc_mass_comp[0, "S_h2"].fix(1e-8)
    m.fs.unit.inlet.conc_mass_comp[0, "S_ch4"].fix(1e-5)
    m.fs.unit.inlet.conc_mass_comp[0, "S_IC"].fix(0.48)
    m.fs.unit.inlet.conc_mass_comp[0, "S_IN"].fix(0.14)
    m.fs.unit.inlet.conc_mass_comp[0, "S_I"].fix(0.02)

    m.fs.unit.inlet.conc_mass_comp[0, "X_c"].fix(2)
    m.fs.unit.inlet.conc_mass_comp[0, "X_ch"].fix(5)
    m.fs.unit.inlet.conc_mass_comp[0, "X_pr"].fix(20)
    m.fs.unit.inlet.conc_mass_comp[0, "X_li"].fix(5)
    m.fs.unit.inlet.conc_mass_comp[0, "X_su"].fix(0.0)
    m.fs.unit.inlet.conc_mass_comp[0, "X_aa"].fix(0.010)
    m.fs.unit.inlet.conc_mass_comp[0, "X_fa"].fix(0.010)
    m.fs.unit.inlet.conc_mass_comp[0, "X_c4"].fix(0.010)
    m.fs.unit.inlet.conc_mass_comp[0, "X_pro"].fix(0.010)
    m.fs.unit.inlet.conc_mass_comp[0, "X_ac"].fix(0.010)
    m.fs.unit.inlet.conc_mass_comp[0, "X_h2"].fix(0.010)
    m.fs.unit.inlet.conc_mass_comp[0, "X_I"].fix(25)

    m.fs.unit.inlet.cations[0].fix(0.04)
    m.fs.unit.inlet.anions[0].fix(0.02)

    m.fs.unit.volume_liquid.fix(3400)
    m.fs.unit.volume_vapor.fix(300)
    m.fs.unit.liquid_outlet.temperature.fix(308.15)

    iscale.calculate_scaling_factors(m.fs.unit)

    # Set scaling factors for badly scaled variables
    iscale.set_scaling_factor(
        m.fs.unit.liquid_phase.mass_transfer_term[0, "Liq", "S_h2"], 1e7
    )
    iscale.set_scaling_factor(m.fs.unit.liquid_phase.heat[0], 1e3)

    return m


class TestAnaerobicDigester(UnitTestHarness):
    def configure(self):
        m = build()

        self.unit_solutions[m.fs.unit.liquid_outlet.pressure[0]] = 101325
        self.unit_solutions[m.fs.unit.liquid_outlet.temperature[0]] = 308.15
        self.unit_solutions[m.fs.unit.liquid_outlet.conc_mass_comp[0, "S_I"]] = (
            0.3287724
        )
        self.unit_solutions[m.fs.unit.liquid_outlet.conc_mass_comp[0, "S_aa"]] = (
            0.00531408
        )
        self.unit_solutions[m.fs.unit.liquid_outlet.conc_mass_comp[0, "S_ac"]] = (
            0.1977833
        )
        self.unit_solutions[m.fs.unit.liquid_outlet.conc_mass_comp[0, "S_bu"]] = (
            0.0132484
        )
        self.unit_solutions[m.fs.unit.liquid_outlet.conc_mass_comp[0, "S_ch4"]] = (
            0.0549707
        )
        self.unit_solutions[m.fs.unit.liquid_outlet.conc_mass_comp[0, "S_fa"]] = (
            0.0986058
        )
        self.unit_solutions[m.fs.unit.liquid_outlet.conc_mass_comp[0, "S_h2"]] = (
            2.35916e-07
        )
        self.unit_solutions[m.fs.unit.liquid_outlet.conc_mass_comp[0, "S_pro"]] = (
            0.0157813
        )
        self.unit_solutions[m.fs.unit.liquid_outlet.conc_mass_comp[0, "S_su"]] = (
            0.01195333
        )
        self.unit_solutions[m.fs.unit.liquid_outlet.conc_mass_comp[0, "S_va"]] = (
            0.011622969
        )
        self.unit_solutions[m.fs.unit.liquid_outlet.conc_mass_comp[0, "X_I"]] = 25.6217
        self.unit_solutions[m.fs.unit.liquid_outlet.conc_mass_comp[0, "X_aa"]] = (
            1.1793147
        )
        self.unit_solutions[m.fs.unit.liquid_outlet.conc_mass_comp[0, "X_ac"]] = (
            0.760653
        )
        self.unit_solutions[m.fs.unit.liquid_outlet.conc_mass_comp[0, "X_c"]] = 0.308718
        self.unit_solutions[m.fs.unit.liquid_outlet.conc_mass_comp[0, "X_c4"]] = (
            0.431974
        )
        self.unit_solutions[m.fs.unit.liquid_outlet.conc_mass_comp[0, "X_ch"]] = (
            0.027947465
        )
        self.unit_solutions[m.fs.unit.liquid_outlet.conc_mass_comp[0, "X_fa"]] = (
            0.2430681
        )
        self.unit_solutions[m.fs.unit.liquid_outlet.conc_mass_comp[0, "X_h2"]] = (
            0.3170629
        )
        self.unit_solutions[m.fs.unit.liquid_outlet.conc_mass_comp[0, "X_li"]] = (
            0.0294834
        )
        self.unit_solutions[m.fs.unit.liquid_outlet.conc_mass_comp[0, "X_pr"]] = (
            0.102574392
        )
        self.unit_solutions[m.fs.unit.liquid_outlet.conc_mass_comp[0, "X_pro"]] = (
            0.137323
        )
        self.unit_solutions[m.fs.unit.liquid_outlet.conc_mass_comp[0, "X_su"]] = (
            0.420219
        )
        self.unit_solutions[m.fs.unit.liquid_outlet.conc_mass_comp[0, "S_IC"]] = (
            1.8320212
        )
        self.unit_solutions[m.fs.unit.liquid_outlet.conc_mass_comp[0, "S_IN"]] = (
            1.8235307
        )
        self.unit_solutions[m.fs.unit.liquid_outlet.anions[0]] = 0.0200033
        self.unit_solutions[m.fs.unit.liquid_outlet.cations[0]] = 0.0400066
        self.unit_solutions[m.fs.unit.vapor_outlet.pressure[0]] = 106659.5225
        self.unit_solutions[m.fs.unit.vapor_outlet.temperature[0]] = 308.15
        self.unit_solutions[m.fs.unit.vapor_outlet.flow_vol[0]] = 0.03249637
        self.unit_solutions[m.fs.unit.vapor_outlet.conc_mass_comp[0, "S_ch4"]] = (
            1.6216465
        )
        self.unit_solutions[m.fs.unit.vapor_outlet.conc_mass_comp[0, "S_co2"]] = (
            0.169417
        )
        self.unit_solutions[m.fs.unit.KH_co2[0]] = 0.02714666
        self.unit_solutions[m.fs.unit.KH_ch4[0]] = 0.001161902
        self.unit_solutions[m.fs.unit.KH_h2[0]] = 0.0007384652
        self.unit_solutions[m.fs.unit.electricity_consumption[0]] = 23.7291667
        self.unit_solutions[m.fs.unit.hydraulic_retention_time[0]] = 1880470.588

        # Conservation check

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.inlet.flow_vol[0],
                "out": (
                    m.fs.unit.liquid_outlet.flow_vol[0] * m.fs.props.dens_mass
                    + m.fs.unit.vapor_outlet.flow_vol[0] * m.fs.props_vap.dens_mass
                )
                / m.fs.props.dens_mass,
            },
            "Check 2": {
                "in": (
                    m.fs.unit.inlet.flow_vol[0]
                    * m.fs.props.dens_mass
                    * m.fs.props.cp_mass
                    * (m.fs.unit.inlet.temperature[0] - m.fs.props.temperature_ref)
                )
                - (
                    m.fs.unit.liquid_outlet.flow_vol[0]
                    * m.fs.props.dens_mass
                    * m.fs.props.cp_mass
                    * (
                        m.fs.unit.liquid_outlet.temperature[0]
                        - m.fs.props.temperature_ref
                    )
                ),
                "out": -1 * m.fs.unit.liquid_phase.enthalpy_transfer[0],
            },
        }

        return m


class TestADScaler:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = ADM1ParameterBlock()
        m.fs.props_vap = ADM1_vaporParameterBlock()
        m.fs.rxn_props = ADM1ReactionParameterBlock(property_package=m.fs.props)

        m.fs.unit = AD(
            liquid_property_package=m.fs.props,
            vapor_property_package=m.fs.props_vap,
            reaction_package=m.fs.rxn_props,
            has_heat_transfer=True,
            has_pressure_change=False,
        )

        # Set the operating conditions
        m.fs.unit.inlet.flow_vol.fix(170 / 24 / 3600)
        m.fs.unit.inlet.temperature.fix(308.15)
        m.fs.unit.inlet.pressure.fix(101325)

        m.fs.unit.inlet.conc_mass_comp[0, "S_su"].fix(0.01)
        m.fs.unit.inlet.conc_mass_comp[0, "S_aa"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_fa"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_va"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_bu"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_pro"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_ac"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_h2"].fix(1e-8)
        m.fs.unit.inlet.conc_mass_comp[0, "S_ch4"].fix(1e-5)
        m.fs.unit.inlet.conc_mass_comp[0, "S_IC"].fix(0.48)
        m.fs.unit.inlet.conc_mass_comp[0, "S_IN"].fix(0.14)
        m.fs.unit.inlet.conc_mass_comp[0, "S_I"].fix(0.02)

        m.fs.unit.inlet.conc_mass_comp[0, "X_c"].fix(2)
        m.fs.unit.inlet.conc_mass_comp[0, "X_ch"].fix(5)
        m.fs.unit.inlet.conc_mass_comp[0, "X_pr"].fix(20)
        m.fs.unit.inlet.conc_mass_comp[0, "X_li"].fix(5)
        m.fs.unit.inlet.conc_mass_comp[0, "X_su"].fix(0.0)
        m.fs.unit.inlet.conc_mass_comp[0, "X_aa"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_fa"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_c4"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_pro"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_ac"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_h2"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_I"].fix(25)

        m.fs.unit.inlet.cations[0].fix(0.04)
        m.fs.unit.inlet.anions[0].fix(0.02)

        m.fs.unit.volume_liquid.fix(3400)
        m.fs.unit.volume_vapor.fix(300)
        m.fs.unit.liquid_outlet.temperature.fix(308.15)

        return m

    @pytest.mark.component
    def test_variable_scaling_routine(self, model):
        scaler = model.fs.unit.default_scaler()

        assert isinstance(scaler, ADScaler)

        scaler.variable_scaling_routine(model.fs.unit)

        # Inlet state
        sfx_in = model.fs.unit.liquid_phase.properties_in[0].scaling_factor
        assert isinstance(sfx_in, Suffix)
        # Scaling factors for FTP and rate reactions
        assert len(sfx_in) == 27

        # Outlet state - should be the same as the inlet
        sfx_out = model.fs.unit.liquid_phase.properties_out[0].scaling_factor
        assert isinstance(sfx_out, Suffix)
        # Scaling factors for FTP and rate reactions
        assert len(sfx_out) == 27

        # Reaction block
        sfx_rxn = model.fs.unit.liquid_phase.reactions[0].scaling_factor
        assert isinstance(sfx_rxn, Suffix)
        # Scaling factors for rxn rate and process inhibition terms
        assert len(sfx_rxn) == 38

        # Check that unit model has scaling factors
        sfx_cv = model.fs.unit.liquid_phase.scaling_factor
        assert isinstance(sfx_cv, Suffix)
        # Scaling factors for volume and rate reactions
        assert len(sfx_cv) == 47

    #
    @pytest.mark.component
    def test_constraint_scaling_routine(self, model):
        scaler = model.fs.unit.default_scaler()

        assert isinstance(scaler, ADScaler)

        scaler.constraint_scaling_routine(model.fs.unit)

        sfx_out = model.fs.unit.liquid_phase.properties_out[0].scaling_factor
        assert isinstance(sfx_out, Suffix)
        assert len(sfx_out) == 0

        sfx_rxn = model.fs.unit.liquid_phase.reactions[0].scaling_factor
        assert isinstance(sfx_rxn, Suffix)
        # Scaling factors for rate expression, process inhibition functions, and other rxn constraints
        assert len(sfx_rxn) == 51

        sfx_unit = model.fs.unit.scaling_factor
        assert isinstance(sfx_unit, Suffix)
        # Scaling factors for material balance, performance equations and other unit model constraints
        assert len(sfx_unit) == 58

    @pytest.mark.component
    def test_scale_model(self, model):
        scaler = model.fs.unit.default_scaler()

        assert isinstance(scaler, ADScaler)

        scaler.scale_model(model.fs.unit)

        # Inlet state
        sfx_in = model.fs.unit.liquid_phase.properties_in[0].scaling_factor
        assert isinstance(sfx_in, Suffix)
        # Scaling factors for FTP and rate reactions
        assert len(sfx_in) == 27

        # Outlet state - should be the same as the inlet
        sfx_out = model.fs.unit.liquid_phase.properties_out[0].scaling_factor
        assert isinstance(sfx_out, Suffix)
        # Scaling factors for FTP and rate reactions
        assert len(sfx_out) == 27

        # Reaction block
        sfx_rxn = model.fs.unit.liquid_phase.reactions[0].scaling_factor
        assert isinstance(sfx_rxn, Suffix)
        # Scaling factors for the combination of constraint and variables sfx_rxn
        assert len(sfx_rxn) == 89

    # TODO: Remove test once iscale is deprecated
    @pytest.mark.integration
    def test_example_case_iscale(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = ADM1ParameterBlock()
        m.fs.props_vap = ADM1_vaporParameterBlock()
        m.fs.rxn_props = ADM1ReactionParameterBlock(property_package=m.fs.props)

        m.fs.unit = AD(
            liquid_property_package=m.fs.props,
            vapor_property_package=m.fs.props_vap,
            reaction_package=m.fs.rxn_props,
            has_heat_transfer=True,
            has_pressure_change=False,
        )

        # Set the operating conditions
        m.fs.unit.inlet.flow_vol.fix(170 / 24 / 3600)
        m.fs.unit.inlet.temperature.fix(308.15)
        m.fs.unit.inlet.pressure.fix(101325)

        m.fs.unit.inlet.conc_mass_comp[0, "S_su"].fix(0.01)
        m.fs.unit.inlet.conc_mass_comp[0, "S_aa"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_fa"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_va"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_bu"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_pro"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_ac"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_h2"].fix(1e-8)
        m.fs.unit.inlet.conc_mass_comp[0, "S_ch4"].fix(1e-5)
        m.fs.unit.inlet.conc_mass_comp[0, "S_IC"].fix(0.48)
        m.fs.unit.inlet.conc_mass_comp[0, "S_IN"].fix(0.14)
        m.fs.unit.inlet.conc_mass_comp[0, "S_I"].fix(0.02)

        m.fs.unit.inlet.conc_mass_comp[0, "X_c"].fix(2)
        m.fs.unit.inlet.conc_mass_comp[0, "X_ch"].fix(5)
        m.fs.unit.inlet.conc_mass_comp[0, "X_pr"].fix(20)
        m.fs.unit.inlet.conc_mass_comp[0, "X_li"].fix(5)
        m.fs.unit.inlet.conc_mass_comp[0, "X_su"].fix(0.0)
        m.fs.unit.inlet.conc_mass_comp[0, "X_aa"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_fa"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_c4"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_pro"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_ac"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_h2"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_I"].fix(25)

        m.fs.unit.inlet.cations[0].fix(0.04)
        m.fs.unit.inlet.anions[0].fix(0.02)

        m.fs.unit.volume_liquid.fix(3400)
        m.fs.unit.volume_vapor.fix(300)
        m.fs.unit.liquid_outlet.temperature.fix(308.15)

        iscale.calculate_scaling_factors(m.fs.unit)

        # Check condition number to confirm scaling
        sm = TransformationFactory("core.scale_model").create_using(m, rename=False)
        jac, _ = get_jacobian(sm, scaled=False)
        assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(
            2.36919186521693e14, rel=1e-3
        )

    @pytest.mark.integration
    def test_example_case_scaler_scaling_default(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = ADM1ParameterBlock()
        m.fs.props_vap = ADM1_vaporParameterBlock()
        m.fs.rxn_props = ADM1ReactionParameterBlock(property_package=m.fs.props)

        m.fs.unit = AD(
            liquid_property_package=m.fs.props,
            vapor_property_package=m.fs.props_vap,
            reaction_package=m.fs.rxn_props,
            has_heat_transfer=True,
            has_pressure_change=False,
        )

        # Set the operating conditions
        m.fs.unit.inlet.flow_vol.fix(170 / 24 / 3600)
        m.fs.unit.inlet.temperature.fix(308.15)
        m.fs.unit.inlet.pressure.fix(101325)

        m.fs.unit.inlet.conc_mass_comp[0, "S_su"].fix(0.01)
        m.fs.unit.inlet.conc_mass_comp[0, "S_aa"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_fa"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_va"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_bu"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_pro"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_ac"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_h2"].fix(1e-8)
        m.fs.unit.inlet.conc_mass_comp[0, "S_ch4"].fix(1e-5)
        m.fs.unit.inlet.conc_mass_comp[0, "S_IC"].fix(0.48)
        m.fs.unit.inlet.conc_mass_comp[0, "S_IN"].fix(0.14)
        m.fs.unit.inlet.conc_mass_comp[0, "S_I"].fix(0.02)

        m.fs.unit.inlet.conc_mass_comp[0, "X_c"].fix(2)
        m.fs.unit.inlet.conc_mass_comp[0, "X_ch"].fix(5)
        m.fs.unit.inlet.conc_mass_comp[0, "X_pr"].fix(20)
        m.fs.unit.inlet.conc_mass_comp[0, "X_li"].fix(5)
        m.fs.unit.inlet.conc_mass_comp[0, "X_su"].fix(0.0)
        m.fs.unit.inlet.conc_mass_comp[0, "X_aa"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_fa"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_c4"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_pro"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_ac"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_h2"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_I"].fix(25)

        m.fs.unit.inlet.cations[0].fix(0.04)
        m.fs.unit.inlet.anions[0].fix(0.02)

        m.fs.unit.volume_liquid.fix(3400)
        m.fs.unit.volume_vapor.fix(300)
        m.fs.unit.liquid_outlet.temperature.fix(308.15)

        scaler = ADScaler()
        scaler.scale_model(
            m.fs.unit,
            submodel_scalers={
                m.fs.unit.liquid_phase.properties_in: ADM1PropertiesScaler,
                m.fs.unit.liquid_phase.properties_out: ADM1PropertiesScaler,
                m.fs.unit.liquid_phase.reactions: ADM1ReactionScaler,
            },
        )

        # Check condition number to confirm scaling
        sm = TransformationFactory("core.scale_model").create_using(m, rename=False)
        jac, _ = get_jacobian(sm, scaled=False)
        assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(
            2.608771605737e12, rel=1e-3
        )

    @pytest.mark.integration
    def test_example_case_scaler_scaling(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = ADM1ParameterBlock()
        m.fs.props_vap = ADM1_vaporParameterBlock()
        m.fs.rxn_props = ADM1ReactionParameterBlock(property_package=m.fs.props)

        m.fs.unit = AD(
            liquid_property_package=m.fs.props,
            vapor_property_package=m.fs.props_vap,
            reaction_package=m.fs.rxn_props,
            has_heat_transfer=True,
            has_pressure_change=False,
        )

        # Set the operating conditions
        m.fs.unit.inlet.flow_vol.fix(170 / 24 / 3600)
        m.fs.unit.inlet.temperature.fix(308.15)
        m.fs.unit.inlet.pressure.fix(101325)

        m.fs.unit.inlet.conc_mass_comp[0, "S_su"].fix(0.01)
        m.fs.unit.inlet.conc_mass_comp[0, "S_aa"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_fa"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_va"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_bu"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_pro"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_ac"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[0, "S_h2"].fix(1e-8)
        m.fs.unit.inlet.conc_mass_comp[0, "S_ch4"].fix(1e-5)
        m.fs.unit.inlet.conc_mass_comp[0, "S_IC"].fix(0.48)
        m.fs.unit.inlet.conc_mass_comp[0, "S_IN"].fix(0.14)
        m.fs.unit.inlet.conc_mass_comp[0, "S_I"].fix(0.02)

        m.fs.unit.inlet.conc_mass_comp[0, "X_c"].fix(2)
        m.fs.unit.inlet.conc_mass_comp[0, "X_ch"].fix(5)
        m.fs.unit.inlet.conc_mass_comp[0, "X_pr"].fix(20)
        m.fs.unit.inlet.conc_mass_comp[0, "X_li"].fix(5)
        m.fs.unit.inlet.conc_mass_comp[0, "X_su"].fix(0.0)
        m.fs.unit.inlet.conc_mass_comp[0, "X_aa"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_fa"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_c4"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_pro"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_ac"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_h2"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[0, "X_I"].fix(25)

        m.fs.unit.inlet.cations[0].fix(0.04)
        m.fs.unit.inlet.anions[0].fix(0.02)

        m.fs.unit.volume_liquid.fix(3400)
        m.fs.unit.volume_vapor.fix(300)
        m.fs.unit.liquid_outlet.temperature.fix(308.15)

        sb = ScalerBase()

        # Apply scaling to unscaled variables
        for var in m.fs.component_data_objects(Var, descend_into=True):
            if "conc_mass_comp" in var.name:
                sb.set_variable_scaling_factor(var, 1e1)
            if "conc_mol" in var.name:
                sb.set_variable_scaling_factor(var, 1e2)
            if "reaction_rate" in var.name:
                sb.set_variable_scaling_factor(var, 1e6)

            sb.set_variable_scaling_factor(m.fs.unit.hydraulic_retention_time[0], 1e-6)

        scaler = ADScaler()
        scaler.scale_model(
            m.fs.unit,
            submodel_scalers={
                m.fs.unit.liquid_phase.properties_in: ADM1PropertiesScaler,
                m.fs.unit.liquid_phase.properties_out: ADM1PropertiesScaler,
                m.fs.unit.liquid_phase.reactions: ADM1ReactionScaler,
            },
        )

        # Check condition number to confirm scaling
        sm = TransformationFactory("core.scale_model").create_using(m, rename=False)
        jac, _ = get_jacobian(sm, scaled=False)
        assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(
            7.0054790813389e11, rel=1e-3
        )
