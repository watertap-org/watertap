.. _how_to_use_unit_test_harness:

How to use the unit test harness
================================

Overview
--------

This guide shows you how to use the unit test harness to generate tests for WaterTAP unit models. The purpose of this
tool is to standardize testing so developers don't need to write tests from the ground-up for each unit model.

How To
------

Begin by importing the following essential functions. This example
assumes a test file is being created for an anaerobic digester.

.. testsetup::

    # quiet idaes logs
    import idaes.logger as idaeslogger
    idaeslogger.getLogger('ideas.core').setLevel('CRITICAL')
    idaeslogger.getLogger('ideas.core.util.scaling').setLevel('CRITICAL')
    idaeslogger.getLogger('idaes.init').setLevel('CRITICAL')

.. testcode::

    import pytest

    from pyomo.environ import ConcreteModel

    from idaes.core import FlowsheetBlock, UnitModelCostingBlock
    from idaes.core.solvers import get_solver
    import idaes.core.util.scaling as iscale

    from watertap.costing import WaterTAPCosting

    # The following imports are unit-model specific
    from watertap.unit_models.anaerobic_digester import AD
    from watertap.property_models.anaerobic_digestion.adm1_properties import ADM1ParameterBlock
    from watertap.property_models.anaerobic_digestion.adm1_properties_vapor import ADM1_vaporParameterBlock
    from watertap.property_models.anaerobic_digestion.adm1_reactions import ADM1ReactionParameterBlock
    from watertap.unit_models.tests.unit_test_harness import UnitTestHarness

    # Get the default solver for testing
    solver = get_solver()

Next, setup the configure function which will create the flowsheet, specify the property and reaction packages,
specify the unit model configuration, set the operating conditions, add the unit model costing, and
set the scaling factors for any variables that are badly scaled. Then, iterate through any variables on the unit model that you'd like to confirm the value of.
Failures may arise at this stage, at which point an error message will be displayed that prompts you
to adjust something in the configure function and/or address the discrepancy between the
expected value for a variable (user-input) and its actual value.

.. testcode::

    class TestUnitDefault(UnitTestHarness):
        @pytest.mark.unit
        def configure(self):
            # Create the ConcreteModel and FlowsheetBlock
            m = ConcreteModel()
            m.fs = FlowsheetBlock(dynamic=False)

            # Set up the property package and a reaction package, if relevant
            m.fs.props = ADM1ParameterBlock()
            m.fs.props_vap = ADM1_vaporParameterBlock()
            m.fs.rxn_props = ADM1ReactionParameterBlock(property_package=m.fs.props)

            # Create the unit model and specify configration options
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

            # Add unit model costing
            m.fs.costing = WaterTAPCosting()

            m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
            m.fs.costing.cost_process()

            # Set scaling factors for badly scaled variables
            iscale.set_scaling_factor(
            m.fs.unit.liquid_phase.mass_transfer_term[0, "Liq", "S_h2"], 1e7
            )
            iscale.set_scaling_factor(m.fs.unit.costing.capital_cost, 1e-6)

            # Specify the unit model being tested
            self.unit_model_block = m.fs.unit

            # Check the expected unit model outputs
            self.unit_solutions[m.fs.unit.liquid_outlet.pressure[0]] = 101325
            self.unit_solutions[m.fs.unit.liquid_outlet.temperature[0]] = 308.15
            self.unit_solutions[
                m.fs.unit.liquid_outlet.conc_mass_comp[0, "S_I"]
            ] = 0.3287724
            self.unit_solutions[
                m.fs.unit.liquid_outlet.conc_mass_comp[0, "S_aa"]
            ] = 0.00531408
            self.unit_solutions[
                m.fs.unit.liquid_outlet.conc_mass_comp[0, "S_ac"]
            ] = 0.1977833
            self.unit_solutions[
                m.fs.unit.liquid_outlet.conc_mass_comp[0, "S_bu"]
            ] = 0.0132484
            self.unit_solutions[
                m.fs.unit.liquid_outlet.conc_mass_comp[0, "S_ch4"]
            ] = 0.0549707
            self.unit_solutions[
                m.fs.unit.liquid_outlet.conc_mass_comp[0, "S_fa"]
            ] = 0.0986058
            self.unit_solutions[
                m.fs.unit.liquid_outlet.conc_mass_comp[0, "S_h2"]
            ] = 2.35916e-07
            self.unit_solutions[
                m.fs.unit.liquid_outlet.conc_mass_comp[0, "S_pro"]
            ] = 0.0157813
            self.unit_solutions[
                m.fs.unit.liquid_outlet.conc_mass_comp[0, "S_su"]
            ] = 0.01195333
            self.unit_solutions[
                m.fs.unit.liquid_outlet.conc_mass_comp[0, "S_va"]
            ] = 0.011622969
            self.unit_solutions[m.fs.unit.liquid_outlet.conc_mass_comp[0, "X_I"]] = 25.6217
            self.unit_solutions[
                m.fs.unit.liquid_outlet.conc_mass_comp[0, "X_aa"]
            ] = 1.1793147
            self.unit_solutions[
                m.fs.unit.liquid_outlet.conc_mass_comp[0, "X_ac"]
            ] = 0.760653
            self.unit_solutions[m.fs.unit.liquid_outlet.conc_mass_comp[0, "X_c"]] = 0.308718
            self.unit_solutions[
                m.fs.unit.liquid_outlet.conc_mass_comp[0, "X_c4"]
            ] = 0.431974
            self.unit_solutions[
                m.fs.unit.liquid_outlet.conc_mass_comp[0, "X_ch"]
            ] = 0.027947465
            self.unit_solutions[
                m.fs.unit.liquid_outlet.conc_mass_comp[0, "X_fa"]
            ] = 0.2430681
            self.unit_solutions[
                m.fs.unit.liquid_outlet.conc_mass_comp[0, "X_h2"]
            ] = 0.3170629
            self.unit_solutions[
                m.fs.unit.liquid_outlet.conc_mass_comp[0, "X_li"]
            ] = 0.0294834
            self.unit_solutions[
                m.fs.unit.liquid_outlet.conc_mass_comp[0, "X_pr"]
            ] = 0.102574392
            self.unit_solutions[
                m.fs.unit.liquid_outlet.conc_mass_comp[0, "X_pro"]
            ] = 0.137323
            self.unit_solutions[
                m.fs.unit.liquid_outlet.conc_mass_comp[0, "X_su"]
            ] = 0.420219
            self.unit_solutions[
                m.fs.unit.liquid_outlet.conc_mass_comp[0, "S_IC"]
            ] = 1.8320212
            self.unit_solutions[
                m.fs.unit.liquid_outlet.conc_mass_comp[0, "S_IN"]
            ] = 1.8235307
            self.unit_solutions[m.fs.unit.liquid_outlet.anions[0]] = 0.0200033
            self.unit_solutions[m.fs.unit.liquid_outlet.cations[0]] = 0.0400066
            self.unit_solutions[m.fs.unit.vapor_outlet.pressure[0]] = 106659.5225
            self.unit_solutions[m.fs.unit.vapor_outlet.temperature[0]] = 308.15
            self.unit_solutions[m.fs.unit.vapor_outlet.flow_vol[0]] = 0.03249637
            self.unit_solutions[
                m.fs.unit.vapor_outlet.conc_mass_comp[0, "S_ch4"]
            ] = 1.6216465
            self.unit_solutions[
                m.fs.unit.vapor_outlet.conc_mass_comp[0, "S_co2"]
            ] = 0.169417
            self.unit_solutions[m.fs.unit.KH_co2[0]] = 0.02714666
            self.unit_solutions[m.fs.unit.KH_ch4[0]] = 0.001161902
            self.unit_solutions[m.fs.unit.KH_h2[0]] = 0.0007384652
            self.unit_solutions[m.fs.unit.electricity_consumption[0]] = 23.7291667
            self.unit_solutions[m.fs.unit.hydraulic_retention_time[0]] = 1880470.588
            self.unit_solutions[m.fs.unit.costing.capital_cost] = 2166581.415

