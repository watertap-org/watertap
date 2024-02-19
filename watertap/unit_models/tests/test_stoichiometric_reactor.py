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
from pyomo.environ import (
    ConcreteModel,
    value,
    check_optimal_termination,
    units as pyunits,
)

from pyomo.util.check_units import assert_units_consistent
from idaes.core import (
    FlowsheetBlock,
    UnitModelCostingBlock,
)
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
)
from idaes.core.util.scaling import (
    calculate_scaling_factors,
)
from idaes.core.util.testing import initialization_tester

from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
    ActivityCoefficientModel,
    DensityCalculation,
    MaterialFlowBasis,
)

from watertap.unit_models.stoichiometric_reactor import StoichiometricReactor
from watertap.costing import WaterTAPCosting

__author__ = "Timothy Bartholomew, Alexander Dudchenko"

solver = get_solver()


class TestStoichiometricReactor:
    @pytest.fixture(scope="class")
    def basic_unit_mass(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        default_dict = {
            "solute_list": [
                "Ca_2+",
                "HCO3_-",
                "Na_+",
                "Cl_-",
                "Mg_2+",
            ],
            "diffusivity_data": {
                ("Liq", "Ca_2+"): 9.2e-10,
                ("Liq", "HCO3_-"): 1.19e-9,
                ("Liq", "Na_+"): 1.33e-9,
                ("Liq", "Cl_-"): 2.03e-9,
                ("Liq", "Mg_2+"): 0.706e-9,
            },
            "mw_data": {
                "H2O": 18e-3,
                "Ca_2+": 40e-3,
                "HCO3_-": 61.0168e-3,
                "Na_+": 23e-3,
                "Cl_-": 35e-3,
                "Mg_2+": 24.305e-3,
            },
            "stokes_radius_data": {
                "Ca_2+": 0.309e-9,
                "HCO3_-": 2.06e-10,
                "Cl_-": 0.121e-9,
                "Na_+": 0.184e-9,
                "Mg_2+": 0.347e-9,
            },
            "charge": {
                "Ca_2+": 2,
                "HCO3_-": -1,
                "Na_+": 1,
                "Cl_-": -1,
                "Mg_2+": 2,
            },
            "activity_coefficient_model": ActivityCoefficientModel.ideal,
            "density_calculation": DensityCalculation.constant,
            "material_flow_basis": MaterialFlowBasis.mass,
        }
        m.fs.properties = MCASParameterBlock(**default_dict)
        # Define precipitatnts
        precipitants = {
            "Calcite": {
                "mw": 100.09 * pyunits.g / pyunits.mol,
                "precipitation_stoichiometric": {"Ca_2+": 1, "HCO3_-": 1},
            },
            "Brucite": {
                "mw": 58.3197 * pyunits.g / pyunits.mol,
                "precipitation_stoichiometric": {"Mg_2+": 1, "H2O": 1},
            },
        }
        # Define reagents
        reagents = {
            "Na2CO3": {
                "mw": 105.99 * pyunits.g / pyunits.mol,
                "dissolution_stoichiometric": {"Na_+": 2, "HCO3_-": 1},
            },
            "CaO": {
                "mw": 56.0774 * pyunits.g / pyunits.mol,
                "dissolution_stoichiometric": {"Ca_2+": 1, "H2O": 1},
            },
        }
        m.fs.unit = StoichiometricReactor(
            property_package=m.fs.properties,
            reagent=reagents,
            precipitate=precipitants,
        )

        m.fs.unit.waste_mass_frac_precipitate.fix(0.2)
        # the reactor us assumed performance model so we are simply testing that
        # mass balances are correct
        m.fs.unit.reagent_dose["Na2CO3"].fix(1e-3)
        m.fs.unit.reagent_dose["CaO"].fix(1e-3)

        m.fs.unit.flow_mass_precipitate["Calcite"].fix(1e-3)
        m.fs.unit.flow_mass_precipitate["Brucite"].fix(1e-4)

        # fix_fix_comp
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(1)
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "Ca_2+"].fix(0.001)
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "HCO3_-"].fix(0.001)
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "Na_+"].fix(0.001)
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "Cl_-"].fix(0.001)
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "Mg_2+"].fix(0.001)

        for key, obj in m.fs.unit.inlet.flow_mass_phase_comp.items():
            val = value(obj)
            m.fs.properties.set_default_scaling(
                "flow_mass_phase_comp", 1 / val, index=("Liq", key[-1])
            )

        m.fs.unit.inlet.pressure[0].fix(101325)
        m.fs.unit.inlet.temperature[0].fix(273.15 + 20)

        return m

    @pytest.fixture(scope="class")
    def basic_unit_molar(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        default_dict = {
            "solute_list": [
                "Ca_2+",
                "HCO3_-",
                "Na_+",
                "Cl_-",
                "Mg_2+",
            ],
            "diffusivity_data": {
                ("Liq", "Ca_2+"): 9.2e-10,
                ("Liq", "HCO3_-"): 1.19e-9,
                ("Liq", "Na_+"): 1.33e-9,
                ("Liq", "Cl_-"): 2.03e-9,
                ("Liq", "Mg_2+"): 0.706e-9,
            },
            "mw_data": {
                "H2O": 18e-3,
                "Ca_2+": 40e-3,
                "HCO3_-": 61.0168e-3,
                "Na_+": 23e-3,
                "Cl_-": 35e-3,
                "Mg_2+": 24.305e-3,
            },
            "stokes_radius_data": {
                "Ca_2+": 0.309e-9,
                "HCO3_-": 2.06e-10,
                "Cl_-": 0.121e-9,
                "Na_+": 0.184e-9,
                "Mg_2+": 0.347e-9,
            },
            "charge": {
                "Ca_2+": 2,
                "HCO3_-": -1,
                "Na_+": 1,
                "Cl_-": -1,
                "Mg_2+": 2,
            },
            "activity_coefficient_model": ActivityCoefficientModel.ideal,
            "density_calculation": DensityCalculation.constant,
            "material_flow_basis": MaterialFlowBasis.molar,
        }
        m.fs.properties = MCASParameterBlock(**default_dict)
        # Define precipitatnts
        precipitants = {
            "Calcite": {
                "mw": 100.09 * pyunits.g / pyunits.mol,
                "precipitation_stoichiometric": {"Ca_2+": 1, "HCO3_-": 1},
            },
            "Brucite": {
                "mw": 58.3197 * pyunits.g / pyunits.mol,
                "precipitation_stoichiometric": {"Mg_2+": 1, "H2O": 1},
            },
        }
        # Define reagents
        reagents = {
            "Na2CO3": {
                "mw": 105.99 * pyunits.g / pyunits.mol,
                "dissolution_stoichiometric": {"Na_+": 2, "HCO3_-": 1},
            },
            "CaO": {
                "mw": 56.0774 * pyunits.g / pyunits.mol,
                "dissolution_stoichiometric": {"Ca_2+": 1, "H2O": 1},
            },
        }
        m.fs.unit = StoichiometricReactor(
            property_package=m.fs.properties,
            reagent=reagents,
            precipitate=precipitants,
        )

        m.fs.unit.waste_mass_frac_precipitate.fix(0.2)
        # the reactor us assumed performance model so we are simply testing that
        # mass balances are correct
        m.fs.unit.reagent_dose["Na2CO3"].fix(1e-3)
        m.fs.unit.reagent_dose["CaO"].fix(1e-3)

        m.fs.unit.flow_mass_precipitate["Calcite"].fix(1e-3)
        m.fs.unit.flow_mass_precipitate["Brucite"].fix(1e-4)

        # fix_fix_comp
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(55)
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Ca_2+"].fix(0.1)
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "HCO3_-"].fix(0.1)
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(0.1)
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(0.1)
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Mg_2+"].fix(0.1)

        for key, obj in m.fs.unit.inlet.flow_mol_phase_comp.items():
            val = value(obj)
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1 / val, index=("Liq", key[-1])
            )

        m.fs.unit.inlet.pressure[0].fix(101325)
        m.fs.unit.inlet.temperature[0].fix(273.15 + 20)

        return m

    @pytest.fixture(scope="class")
    def dissolution_reactor(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        default_dict = {
            "solute_list": [
                "Ca_2+",
                "HCO3_-",
                "Na_+",
                "Cl_-",
                "Mg_2+",
            ],
            "diffusivity_data": {
                ("Liq", "Ca_2+"): 9.2e-10,
                ("Liq", "HCO3_-"): 1.19e-9,
                ("Liq", "Na_+"): 1.33e-9,
                ("Liq", "Cl_-"): 2.03e-9,
                ("Liq", "Mg_2+"): 0.706e-9,
            },
            "mw_data": {
                "H2O": 18e-3,
                "Ca_2+": 40e-3,
                "HCO3_-": 61.0168e-3,
                "Na_+": 23e-3,
                "Cl_-": 35e-3,
                "Mg_2+": 24.305e-3,
            },
            "stokes_radius_data": {
                "Ca_2+": 0.309e-9,
                "HCO3_-": 2.06e-10,
                "Cl_-": 0.121e-9,
                "Na_+": 0.184e-9,
                "Mg_2+": 0.347e-9,
            },
            "charge": {
                "Ca_2+": 2,
                "HCO3_-": -1,
                "Na_+": 1,
                "Cl_-": -1,
                "Mg_2+": 2,
            },
            "activity_coefficient_model": ActivityCoefficientModel.ideal,
            "density_calculation": DensityCalculation.constant,
            "material_flow_basis": MaterialFlowBasis.molar,
        }
        m.fs.properties = MCASParameterBlock(**default_dict)
        # Define reagents
        reagents = {
            "Na2CO3": {
                "mw": 105.99 * pyunits.g / pyunits.mol,
                "dissolution_stoichiometric": {"Na_+": 2, "HCO3_-": 1},
                "reagent_density": 1.2 * pyunits.kg / pyunits.L,
            },
            "CaO": {
                "mw": 56.0774 * pyunits.g / pyunits.mol,
                "dissolution_stoichiometric": {"Ca_2+": 1, "H2O": 1},
                "reagent_density": 1.2 * pyunits.kg / pyunits.L,
            },
        }
        m.fs.unit = StoichiometricReactor(
            property_package=m.fs.properties,
            reagent=reagents,
        )

        # the reactor us assumed performance model so we are simply testing that
        # mass balances are correct
        m.fs.unit.reagent_dose["Na2CO3"].fix(1e-3)
        m.fs.unit.reagent_dose["CaO"].fix(1e-3)
        # fix_fix_comp
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(55)
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Ca_2+"].fix(0.1)
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "HCO3_-"].fix(0.1)
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(0.1)
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(0.1)
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Mg_2+"].fix(0.1)

        for key, obj in m.fs.unit.inlet.flow_mol_phase_comp.items():
            val = value(obj)
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1 / val, index=("Liq", key[-1])
            )

        m.fs.unit.inlet.pressure[0].fix(101325)
        m.fs.unit.inlet.temperature[0].fix(273.15 + 20)

        return m

    @pytest.fixture(scope="class")
    def precipitation_reactor(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        default_dict = {
            "solute_list": [
                "Ca_2+",
                "HCO3_-",
                "Na_+",
                "Cl_-",
                "Mg_2+",
            ],
            "diffusivity_data": {
                ("Liq", "Ca_2+"): 9.2e-10,
                ("Liq", "HCO3_-"): 1.19e-9,
                ("Liq", "Na_+"): 1.33e-9,
                ("Liq", "Cl_-"): 2.03e-9,
                ("Liq", "Mg_2+"): 0.706e-9,
            },
            "mw_data": {
                "H2O": 18e-3,
                "Ca_2+": 40e-3,
                "HCO3_-": 61.0168e-3,
                "Na_+": 23e-3,
                "Cl_-": 35e-3,
                "Mg_2+": 24.305e-3,
            },
            "stokes_radius_data": {
                "Ca_2+": 0.309e-9,
                "HCO3_-": 2.06e-10,
                "Cl_-": 0.121e-9,
                "Na_+": 0.184e-9,
                "Mg_2+": 0.347e-9,
            },
            "charge": {
                "Ca_2+": 2,
                "HCO3_-": -1,
                "Na_+": 1,
                "Cl_-": -1,
                "Mg_2+": 2,
            },
            "activity_coefficient_model": ActivityCoefficientModel.ideal,
            "density_calculation": DensityCalculation.constant,
            "material_flow_basis": MaterialFlowBasis.molar,
        }
        m.fs.properties = MCASParameterBlock(**default_dict)
        # Define reagents
        # Define reagents
        precipitants = {
            "Calcite": {
                "mw": 100.09 * pyunits.g / pyunits.mol,
                "precipitation_stoichiometric": {"Ca_2+": 1, "HCO3_-": 1},
            },
            "Brucite": {
                "mw": 58.3197 * pyunits.g / pyunits.mol,
                "precipitation_stoichiometric": {"Mg_2+": 1, "H2O": 1},
            },
        }
        m.fs.unit = StoichiometricReactor(
            property_package=m.fs.properties,
            precipitate=precipitants,
        )

        m.fs.unit.waste_mass_frac_precipitate.fix(0.2)
        # the reactor us assumed performance model so we are simply testing that
        # mass balances are correct

        m.fs.unit.flow_mass_precipitate["Calcite"].fix(1e-3)
        m.fs.unit.flow_mass_precipitate["Brucite"].fix(1e-4)
        m.fs.unit.waste_mass_frac_precipitate.fix(0.2)

        # fix_fix_comp
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(55)
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Ca_2+"].fix(0.1)
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "HCO3_-"].fix(0.1)
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(0.1)
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(0.1)
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Mg_2+"].fix(0.1)

        for key, obj in m.fs.unit.inlet.flow_mol_phase_comp.items():
            val = value(obj)
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1 / val, index=("Liq", key[-1])
            )

        m.fs.unit.inlet.pressure[0].fix(101325)
        m.fs.unit.inlet.temperature[0].fix(273.15 + 20)

        return m

    @pytest.mark.unit
    def test_dof_units(self, basic_unit_molar, basic_unit_mass):
        m = basic_unit_molar
        assert degrees_of_freedom(m) == 0
        assert assert_units_consistent(m) is None
        m = basic_unit_mass
        assert degrees_of_freedom(m) == 0
        assert assert_units_consistent(m) is None

    @pytest.mark.unit
    def test_molar(self, basic_unit_molar):
        m = basic_unit_molar
        calculate_scaling_factors(m)

        initialization_tester(m)

        flow_ca_in_reagent = (
            1e-3
            / (56.0774 * 1e-3)
            * value(
                m.fs.unit.dissolution_reactor.properties_in[0].flow_vol_phase["Liq"]
            )
        )
        expected_ca_mol_flow = 0.1 + flow_ca_in_reagent

        assert pytest.approx(flow_ca_in_reagent, rel=1e-5) == value(
            m.fs.unit.flow_mass_reagent["CaO"] / (56.0774 * 1e-3)
        )

        assert pytest.approx(expected_ca_mol_flow, rel=1e-5) == value(
            m.fs.unit.dissolution_reactor.properties_out[0].flow_mol_phase_comp[
                "Liq", "Ca_2+"
            ]
        )

        expected_ca_mol_flow = 0.1 + flow_ca_in_reagent - 1e-3 / (100.09 * 1e-3)
        assert pytest.approx(expected_ca_mol_flow, rel=1e-5) == value(
            m.fs.unit.precipitation_reactor.properties_out[0].flow_mol_phase_comp[
                "Liq", "Ca_2+"
            ]
        )

        flow_na_in_reagent = (
            1e-3
            / (105.99 * 1e-3)
            * value(
                m.fs.unit.dissolution_reactor.properties_in[0].flow_vol_phase["Liq"]
            )
        )

        expected_na_mol_flow = 0.1 + flow_na_in_reagent * 2

        assert pytest.approx(flow_na_in_reagent, rel=1e-5) == value(
            m.fs.unit.flow_mass_reagent["Na2CO3"] / (105.99 * 1e-3)
        )

        assert pytest.approx(expected_na_mol_flow, rel=1e-5) == value(
            m.fs.unit.dissolution_reactor.properties_out[0].flow_mol_phase_comp[
                "Liq", "Na_+"
            ]
        )

        expected_mg_mol_flow = 0.1 - 1e-4 / (58.3197 * 1e-3)

        assert pytest.approx(expected_mg_mol_flow, rel=1e-5) == value(
            m.fs.unit.precipitation_reactor.properties_out[0].flow_mol_phase_comp[
                "Liq", "Mg_2+"
            ]
        )

    @pytest.mark.unit
    def test_mass(self, basic_unit_mass):
        m = basic_unit_mass
        calculate_scaling_factors(m)

        initialization_tester(m)

        flow_ca_in_reagent = 1e-3 * value(
            m.fs.unit.dissolution_reactor.properties_in[0].flow_vol_phase["Liq"]
        )
        expected_ca_flow = (
            0.001 + flow_ca_in_reagent / 56.0774 * 40
        )  # - 1e-3 / 100.9 * 40

        assert pytest.approx(flow_ca_in_reagent, rel=1e-2) == value(
            m.fs.unit.flow_mass_reagent["CaO"]
        )

        assert pytest.approx(expected_ca_flow, rel=1e-2) == value(
            m.fs.unit.dissolution_reactor.properties_out[0].flow_mass_phase_comp[
                "Liq", "Ca_2+"
            ]
        )
        expected_ca_flow = (
            0.001 + flow_ca_in_reagent / 56.0774 * 40
        ) - 1e-3 / 100.9 * 40
        assert pytest.approx(expected_ca_flow, rel=1e-2) == value(
            m.fs.unit.precipitation_reactor.properties_out[0].flow_mass_phase_comp[
                "Liq", "Ca_2+"
            ]
        )
        flow_na_in_reagent = 1e-3 * value(
            m.fs.unit.dissolution_reactor.properties_in[0].flow_vol_phase["Liq"]
        )

        expected_na_flow = 0.001 + flow_na_in_reagent / 105.99 * 2 * 23

        assert pytest.approx(flow_na_in_reagent, rel=1e-5) == value(
            m.fs.unit.flow_mass_reagent["Na2CO3"]
        )

        assert pytest.approx(expected_na_flow, rel=1e-5) == value(
            m.fs.unit.dissolution_reactor.properties_out[0].flow_mass_phase_comp[
                "Liq", "Na_+"
            ]
        )

        expected_mg_flow = 0.001 - 1e-4 / 58.3197 * 24.305

        assert pytest.approx(expected_mg_flow, rel=1e-5) == value(
            m.fs.unit.precipitation_reactor.properties_out[0].flow_mass_phase_comp[
                "Liq", "Mg_2+"
            ]
        )

    @pytest.mark.unit
    def test_dissolution(self, dissolution_reactor):
        m = dissolution_reactor
        calculate_scaling_factors(m)

        initialization_tester(m)

        flow_ca_in_reagent = (
            1e-3
            / (56.0774 * 1e-3)
            * value(
                m.fs.unit.dissolution_reactor.properties_in[0].flow_vol_phase["Liq"]
            )
        )
        expected_ca_mol_flow = 0.1 + flow_ca_in_reagent

        assert pytest.approx(flow_ca_in_reagent, rel=1e-5) == value(
            m.fs.unit.flow_mass_reagent["CaO"] / (56.0774 * 1e-3)
        )

        assert pytest.approx(expected_ca_mol_flow, rel=1e-5) == value(
            m.fs.unit.dissolution_reactor.properties_out[0].flow_mol_phase_comp[
                "Liq", "Ca_2+"
            ]
        )

        flow_na_in_reagent = (
            1e-3
            / (105.99 * 1e-3)
            * value(
                m.fs.unit.dissolution_reactor.properties_in[0].flow_vol_phase["Liq"]
            )
        )

        expected_na_mol_flow = 0.1 + flow_na_in_reagent * 2

        assert pytest.approx(flow_na_in_reagent, rel=1e-5) == value(
            m.fs.unit.flow_mass_reagent["Na2CO3"] / (105.99 * 1e-3)
        )

        assert pytest.approx(expected_na_mol_flow, rel=1e-5) == value(
            m.fs.unit.dissolution_reactor.properties_out[0].flow_mol_phase_comp[
                "Liq", "Na_+"
            ]
        )

    @pytest.mark.unit
    def test_precipitation(self, precipitation_reactor):
        m = precipitation_reactor
        calculate_scaling_factors(m)

        initialization_tester(m)

        expected_ca_mol_flow = 0.1 - 1e-3 / (100.09 * 1e-3)
        assert pytest.approx(expected_ca_mol_flow, rel=1e-5) == value(
            m.fs.unit.precipitation_reactor.properties_out[0].flow_mol_phase_comp[
                "Liq", "Ca_2+"
            ]
        )

        expected_mg_mol_flow = 0.1 - 1e-4 / (58.3197 * 1e-3)

        assert pytest.approx(expected_mg_mol_flow, rel=1e-5) == value(
            m.fs.unit.precipitation_reactor.properties_out[0].flow_mol_phase_comp[
                "Liq", "Mg_2+"
            ]
        )

    @pytest.mark.unit
    def test_failed_build(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        default_dict = {
            "solute_list": [
                "Ca_2+",
                "HCO3_-",
                "Na_+",
                "Cl_-",
                "Mg_2+",
            ],
            "diffusivity_data": {
                ("Liq", "Ca_2+"): 9.2e-10,
                ("Liq", "HCO3_-"): 1.19e-9,
                ("Liq", "Na_+"): 1.33e-9,
                ("Liq", "Cl_-"): 2.03e-9,
                ("Liq", "Mg_2+"): 0.706e-9,
            },
            "mw_data": {
                "H2O": 18e-3,
                "Ca_2+": 40e-3,
                "HCO3_-": 61.0168e-3,
                "Na_+": 23e-3,
                "Cl_-": 35e-3,
                "Mg_2+": 24.305e-3,
            },
            "stokes_radius_data": {
                "Ca_2+": 0.309e-9,
                "HCO3_-": 2.06e-10,
                "Cl_-": 0.121e-9,
                "Na_+": 0.184e-9,
                "Mg_2+": 0.347e-9,
            },
            "charge": {
                "Ca_2+": 2,
                "HCO3_-": -1,
                "Na_+": 1,
                "Cl_-": -1,
                "Mg_2+": 2,
            },
            "activity_coefficient_model": ActivityCoefficientModel.ideal,
            "density_calculation": DensityCalculation.constant,
            "material_flow_basis": MaterialFlowBasis.molar,
        }
        m.fs.properties = MCASParameterBlock(**default_dict)
        with pytest.raises(TypeError):
            m.fs.unit = StoichiometricReactor(
                property_package=m.fs.properties,
            )

    @pytest.mark.unit
    def test_costing_softening(self, basic_unit_molar):
        m = basic_unit_molar
        m.fs.costing = WaterTAPCosting()
        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.initialize()
        calculate_scaling_factors(m)
        result = solver.solve(m)
        assert check_optimal_termination(result)
        assert assert_units_consistent(m) is None
        assert degrees_of_freedom(m) == 0

        assert (
            value(m.fs.costing.stoichiometric_reactor.capital_cost_softening) == 374.9
        )
        assert pytest.approx(value(m.fs.costing.total_capital_cost), rel=1e-3) == 245.34
        assert pytest.approx(value(m.fs.costing.total_operating_cost), rel=1e-3) == 7.36

    @pytest.mark.unit
    def test_costing_acid_addition(self, dissolution_reactor):
        # NOTE: testing costing for a dissultion reactor only - we are useing
        # soda ash adn lime addition, but costing is for HCl addition
        # just a test
        m = dissolution_reactor
        m.fs.costing = WaterTAPCosting()
        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.initialize()
        calculate_scaling_factors(m)
        result = solver.solve(m)
        assert check_optimal_termination(result)
        assert assert_units_consistent(m) is None
        assert degrees_of_freedom(m) == 0

        assert (
            value(m.fs.costing.stoichiometric_reactor.capital_cost_acid_addition)
            == 127.8
        )
        assert pytest.approx(value(m.fs.costing.total_capital_cost), rel=1e-3) == 10.02
        assert (
            pytest.approx(value(m.fs.costing.total_operating_cost), rel=1e-3) == 0.30065
        )
