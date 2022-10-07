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

# Imports from idaes core
from idaes.core import AqueousPhase
from idaes.core.base.components import Solvent, Solute, Cation, Anion
from idaes.core.base.phases import PhaseType as PT

# Imports from idaes generic models
from idaes.models.properties.modular_properties.pure.ConstantProperties import Constant
from idaes.models.properties.modular_properties.state_definitions import FpcTP
from idaes.models.properties.modular_properties.eos.ideal import Ideal

# Importing the enum for concentration unit basis used in the 'get_concentration_term' function
from idaes.models.properties.modular_properties.base.generic_reaction import (
    ConcentrationForm,
)

# Import the object/function for heat of reaction
from idaes.models.properties.modular_properties.reactions.dh_rxn import constant_dh_rxn

# Import safe log power law equation
from idaes.models.properties.modular_properties.reactions.equilibrium_forms import (
    log_power_law_equil,
)

# Import k-value functions
from idaes.models.properties.modular_properties.reactions.equilibrium_constant import (
    van_t_hoff,
)

# Import the idaes objects for Generic Properties and Reactions
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)

from watertap.property_models.ion_DSPMDE_prop_pack import DSPMDEParameterBlock

from watertap.unit_models.boron_removal import BoronRemoval
from pyomo.environ import (
    ConcreteModel,
    assert_optimal_termination,
    value,
    Set,
    Param,
    Var,
    units as pyunits,
    Suffix,
    Constraint,
    SolverFactory,
    SolverStatus,
    TerminationCondition,
    log10,
)
from idaes.core import (
    FlowsheetBlock,
)
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
import idaes.core.util.scaling as iscale
from idaes.core.util.scaling import badly_scaled_var_generator
from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog
import re

__author__ = "Austin Ladshaw"

solver = get_solver()

# Helper function for multiple test setup
def model_setup(
    m,
    chem_add=1.8e-5,
    state={
        "H2O": 100,
        "H_+": 1e-7,
        "OH_-": 1e-7,
        "B[OH]3": 2e-4,
        "B[OH]4_-": 1e-6,
        "Na_+": 1e-4,
        "HCO3_-": 1e-4,
    },
):

    m.fs.unit.inlet.pressure.fix(101325)
    m.fs.unit.inlet.temperature.fix(298.15)
    for j in state:
        idx = (0, "Liq", j)
        if idx in m.fs.unit.inlet.flow_mol_phase_comp:
            m.fs.unit.inlet.flow_mol_phase_comp[idx].fix(state[j])
    m.fs.unit.caustic_dose_rate.fix(chem_add)
    m.fs.unit.reactor_volume.fix(1)


# Helper function to automate scaling
def scaling_setup(
    m,
    state={
        "H2O": 100,
        "H_+": 5e-5,
        "OH_-": 5e-5,
        "B[OH]3": 2e-4,
        "B[OH]4_-": 1e-6,
        "Na_+": 1e-4,
        "HCO3_-": 1e-4,
    },
):
    # Set some scaling factors and look for 'bad' scaling
    for j in state:
        idx = (0, "Liq", j)
        if idx in m.fs.unit.inlet.flow_mol_phase_comp:
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1 / state[j], index=("Liq", j)
            )

    iscale.calculate_scaling_factors(m.fs)


# -----------------------------------------------------------------------------
# Start test class
class TestBoronRemoval_IonPropPack_Min:
    @pytest.fixture(scope="class")
    def min_boron_removal_model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        # create dict to define ions (the prop pack requires this)
        ion_dict = {
            "solute_list": ["B[OH]3", "B[OH]4_-", "Na_+"],
            "mw_data": {
                "H2O": 18e-3,
                "B[OH]3": 61.83e-3,
                "B[OH]4_-": 78.83e-3,
                "Na_+": 23e-3,
            },
            "charge": {
                "B[OH]4_-": -1,
                "Na_+": 1,
            },
        }

        # attach prop pack to flowsheet
        m.fs.properties = DSPMDEParameterBlock(**ion_dict)

        map = {
            "boron_name": "B[OH]3",  # [is required]
            "borate_name": "B[OH]4_-",  # [is required]
            "caustic_additive": {
                "cation_name": "Na_+",  # [is optional]
                "mw_additive": (40, pyunits.g / pyunits.mol),  # [is required]
                "moles_cation_per_additive": 1,  # [is required]
            },
        }
        m.fs.unit = BoronRemoval(
            property_package=m.fs.properties, chemical_mapping_data=map
        )

        return m

    @pytest.mark.unit
    def test_build_model(self, min_boron_removal_model):
        m = min_boron_removal_model

        assert isinstance(m.fs.unit.ion_charge, Param)
        assert len(m.fs.unit.ion_charge) == 2

        assert hasattr(m.fs.unit, "boron_name_id")
        assert hasattr(m.fs.unit, "borate_name_id")
        assert hasattr(m.fs.unit, "proton_name_id")
        assert hasattr(m.fs.unit, "hydroxide_name_id")
        assert hasattr(m.fs.unit, "cation_name_id")

        assert hasattr(m.fs.unit, "caustic_chem_name")

        assert hasattr(m.fs.unit, "control_volume")

        assert isinstance(m.fs.unit.caustic_mw, Param)
        assert isinstance(m.fs.unit.additive_molar_ratio, Param)
        assert isinstance(m.fs.unit.caustic_dose_rate, Var)

        assert isinstance(m.fs.unit.Kw_0, Param)
        assert isinstance(m.fs.unit.dH_w, Param)
        assert isinstance(m.fs.unit.Ka_0, Param)
        assert isinstance(m.fs.unit.dH_a, Param)

        assert isinstance(m.fs.unit.conc_mol_H, Var)
        assert isinstance(m.fs.unit.conc_mol_OH, Var)
        assert isinstance(m.fs.unit.conc_mol_Boron, Var)
        assert isinstance(m.fs.unit.conc_mol_Borate, Var)

        assert isinstance(m.fs.unit.eq_mass_transfer_term, Constraint)
        assert isinstance(m.fs.unit.eq_electroneutrality, Constraint)
        assert isinstance(m.fs.unit.eq_total_boron, Constraint)
        assert isinstance(m.fs.unit.eq_water_dissociation, Constraint)
        assert isinstance(m.fs.unit.eq_boron_dissociation, Constraint)

    @pytest.mark.unit
    def test_stats(self, min_boron_removal_model):
        m = min_boron_removal_model

        # Check to make sure we have the correct DOF and
        #   check to make sure the units are correct
        assert_units_consistent(m)
        assert degrees_of_freedom(m) == 8

        # set the variables
        state = {
            "H2O": 100,
            "H_+": 1e-7,
            "OH_-": 1e-7,
            "B[OH]3": 2e-4,
            "B[OH]4_-": 1e-6,
            "Na_+": 1e-15,
            "HCO3_-": 1e-4,
        }
        # For this test, we don't actually know how much
        #   of the base to add. Instead, we have a target
        #   exit flow of boron (see below). We will initially
        #   guess that we want 5 mg/L of additive, but the
        #   actual solution is 10 mg/L.
        model_setup(m, chem_add=0.9e-5, state=state)

        # Modified this test to fix the desired outlet flow
        #   of Boron and unfix the dosage needed to get that outlet
        m.fs.unit.outlet.flow_mol_phase_comp[(0, "Liq", "B[OH]3")].fix(1.98677e-5)
        m.fs.unit.caustic_dose_rate.unfix()

        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_scaling(self, min_boron_removal_model):
        m = min_boron_removal_model

        # Set some scaling factors and look for 'bad' scaling
        scaling_setup(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(iscale.unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

        # check if any variables are badly scaled
        badly_scaled_var_values = {
            var.name: val
            for (var, val) in iscale.badly_scaled_var_generator(
                m, large=1e3, small=1e-3
            )
        }
        assert not badly_scaled_var_values

    @pytest.mark.component
    def test_initialization(self, min_boron_removal_model):
        m = min_boron_removal_model
        initialization_tester(m)

        # check to make sure DOF does not change
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_solve(self, min_boron_removal_model):
        m = min_boron_removal_model

        # first, check to make sure that after initialized, the scaling is still good
        badly_scaled_var_values = {
            var.name: val
            for (var, val) in iscale.badly_scaled_var_generator(
                m, large=1e3, small=1e-3
            )
        }
        assert not badly_scaled_var_values

        # run solver and check for optimal solution
        results = solver.solve(m, tee=True)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, min_boron_removal_model):
        m = min_boron_removal_model

        assert value(
            m.fs.unit.outlet.flow_mol_phase_comp[0, "Liq", m.fs.unit.boron_name_id]
        ) == pytest.approx(1.98677e-5, rel=1e-4)
        assert value(
            m.fs.unit.outlet.flow_mol_phase_comp[0, "Liq", m.fs.unit.borate_name_id]
        ) == pytest.approx(1.81133e-4, rel=1e-4)
        assert value(m.fs.unit.outlet_pH()) == pytest.approx(10.171, rel=1e-4)
        assert value(m.fs.unit.outlet_pOH()) == pytest.approx(3.8257, rel=1e-4)
        assert value(m.fs.unit.caustic_dose_rate[0].value) == pytest.approx(
            1.8e-5, rel=1e-4
        )


# -----------------------------------------------------------------------------
# Start test class
class TestBoronRemoval_IonPropPack_with_ResAlk:
    @pytest.fixture(scope="class")
    def alk_boron_removal_model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        # create dict to define ions (the prop pack requires this)
        ion_dict = {
            "solute_list": ["B[OH]3", "B[OH]4_-", "Na_+", "HCO3_-"],
            "mw_data": {
                "H2O": 18e-3,
                "B[OH]3": 61.83e-3,
                "B[OH]4_-": 78.83e-3,
                "Na_+": 23e-3,
                "HCO3_-": 61e-3,
            },
            "charge": {"B[OH]4_-": -1, "Na_+": 1, "HCO3_-": -1},
        }

        # attach prop pack to flowsheet
        m.fs.properties = DSPMDEParameterBlock(**ion_dict)

        map = {
            "boron_name": "B[OH]3",  # [is required]
            "borate_name": "B[OH]4_-",  # [is required]
            "caustic_additive": {
                "cation_name": "Na_+",  # [is optional]
                "mw_additive": (40, pyunits.g / pyunits.mol),  # [is required]
                "moles_cation_per_additive": 1,  # [is required]
            },
        }
        m.fs.unit = BoronRemoval(
            property_package=m.fs.properties, chemical_mapping_data=map
        )

        return m

    @pytest.mark.unit
    def test_build_model(self, alk_boron_removal_model):
        m = alk_boron_removal_model

        assert isinstance(m.fs.unit.ion_charge, Param)
        assert len(m.fs.unit.ion_charge) == 3

        assert hasattr(m.fs.unit, "boron_name_id")
        assert hasattr(m.fs.unit, "borate_name_id")
        assert hasattr(m.fs.unit, "proton_name_id")
        assert hasattr(m.fs.unit, "hydroxide_name_id")
        assert hasattr(m.fs.unit, "cation_name_id")

        assert hasattr(m.fs.unit, "caustic_chem_name")

        assert hasattr(m.fs.unit, "control_volume")

        assert isinstance(m.fs.unit.caustic_mw, Param)
        assert isinstance(m.fs.unit.additive_molar_ratio, Param)
        assert isinstance(m.fs.unit.caustic_dose_rate, Var)

        assert isinstance(m.fs.unit.Kw_0, Param)
        assert isinstance(m.fs.unit.dH_w, Param)
        assert isinstance(m.fs.unit.Ka_0, Param)
        assert isinstance(m.fs.unit.dH_a, Param)

        assert isinstance(m.fs.unit.conc_mol_H, Var)
        assert isinstance(m.fs.unit.conc_mol_OH, Var)
        assert isinstance(m.fs.unit.conc_mol_Boron, Var)
        assert isinstance(m.fs.unit.conc_mol_Borate, Var)

        assert isinstance(m.fs.unit.eq_mass_transfer_term, Constraint)
        assert isinstance(m.fs.unit.eq_electroneutrality, Constraint)
        assert isinstance(m.fs.unit.eq_total_boron, Constraint)
        assert isinstance(m.fs.unit.eq_water_dissociation, Constraint)
        assert isinstance(m.fs.unit.eq_boron_dissociation, Constraint)

    @pytest.mark.unit
    def test_stats(self, alk_boron_removal_model):
        m = alk_boron_removal_model

        # Check to make sure we have the correct DOF and
        #   check to make sure the units are correct
        assert_units_consistent(m)
        assert degrees_of_freedom(m) == 9

        # set the variables
        model_setup(m, chem_add=1.8e-4)

        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_scaling(self, alk_boron_removal_model):
        m = alk_boron_removal_model

        # Set some scaling factors and look for 'bad' scaling
        scaling_setup(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(iscale.unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

        # check if any variables are badly scaled
        badly_scaled_var_values = {
            var.name: val
            for (var, val) in iscale.badly_scaled_var_generator(
                m, large=1e3, small=1e-3
            )
        }
        assert not badly_scaled_var_values

    @pytest.mark.component
    def test_initialization(self, alk_boron_removal_model):
        m = alk_boron_removal_model
        initialization_tester(m)

        # check to make sure DOF does not change
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_solve(self, alk_boron_removal_model):
        m = alk_boron_removal_model

        # first, check to make sure that after initialized, the scaling is still good
        badly_scaled_var_values = {
            var.name: val
            for (var, val) in iscale.badly_scaled_var_generator(
                m, large=1e3, small=1e-3
            )
        }
        assert not badly_scaled_var_values

        # run solver and check for optimal solution
        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, alk_boron_removal_model):
        m = alk_boron_removal_model

        assert value(
            m.fs.unit.outlet.flow_mol_phase_comp[0, "Liq", m.fs.unit.boron_name_id]
        ) == pytest.approx(1.36911e-6, rel=1e-4)
        assert value(
            m.fs.unit.outlet.flow_mol_phase_comp[0, "Liq", m.fs.unit.borate_name_id]
        ) == pytest.approx(1.99642e-4, rel=1e-4)
        assert value(m.fs.unit.outlet_pH()) == pytest.approx(11.375, rel=1e-4)
        assert value(m.fs.unit.outlet_pOH()) == pytest.approx(2.6218, rel=1e-4)


# -----------------------------------------------------------------------------
# Start test class
class TestBoronRemoval_IonPropPack_with_ResBase:
    @pytest.fixture(scope="class")
    def base_boron_removal_model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        # create dict to define ions (the prop pack requires this)
        ion_dict = {
            "solute_list": ["B[OH]3", "B[OH]4_-", "Na_+"],
            "mw_data": {
                "H2O": 18e-3,
                "B[OH]3": 61.83e-3,
                "B[OH]4_-": 78.83e-3,
                "Na_+": 23e-3,
            },
            "charge": {
                "B[OH]4_-": -1,
                "Na_+": 1,
            },
        }

        # attach prop pack to flowsheet
        m.fs.properties = DSPMDEParameterBlock(**ion_dict)

        map = {
            "boron_name": "B[OH]3",  # [is required]
            "borate_name": "B[OH]4_-",  # [is required]
            "caustic_additive": {
                "cation_name": "Na_+",  # [is optional]
                "mw_additive": (40, pyunits.g / pyunits.mol),  # [is required]
                "moles_cation_per_additive": 1,  # [is required]
            },
        }
        m.fs.unit = BoronRemoval(
            property_package=m.fs.properties, chemical_mapping_data=map
        )

        return m

    @pytest.mark.unit
    def test_build_model(self, base_boron_removal_model):
        m = base_boron_removal_model

        assert isinstance(m.fs.unit.ion_charge, Param)
        assert len(m.fs.unit.ion_charge) == 2

        assert hasattr(m.fs.unit, "boron_name_id")
        assert hasattr(m.fs.unit, "borate_name_id")
        assert hasattr(m.fs.unit, "proton_name_id")
        assert hasattr(m.fs.unit, "hydroxide_name_id")
        assert hasattr(m.fs.unit, "cation_name_id")

        assert hasattr(m.fs.unit, "caustic_chem_name")

        assert hasattr(m.fs.unit, "control_volume")

        assert isinstance(m.fs.unit.caustic_mw, Param)
        assert isinstance(m.fs.unit.additive_molar_ratio, Param)
        assert isinstance(m.fs.unit.caustic_dose_rate, Var)

        assert isinstance(m.fs.unit.Kw_0, Param)
        assert isinstance(m.fs.unit.dH_w, Param)
        assert isinstance(m.fs.unit.Ka_0, Param)
        assert isinstance(m.fs.unit.dH_a, Param)

        assert isinstance(m.fs.unit.conc_mol_H, Var)
        assert isinstance(m.fs.unit.conc_mol_OH, Var)
        assert isinstance(m.fs.unit.conc_mol_Boron, Var)
        assert isinstance(m.fs.unit.conc_mol_Borate, Var)

        assert isinstance(m.fs.unit.eq_mass_transfer_term, Constraint)
        assert isinstance(m.fs.unit.eq_electroneutrality, Constraint)
        assert isinstance(m.fs.unit.eq_total_boron, Constraint)
        assert isinstance(m.fs.unit.eq_water_dissociation, Constraint)
        assert isinstance(m.fs.unit.eq_boron_dissociation, Constraint)

    @pytest.mark.unit
    def test_stats(self, base_boron_removal_model):
        m = base_boron_removal_model

        # Check to make sure we have the correct DOF and
        #   check to make sure the units are correct
        assert_units_consistent(m)
        assert degrees_of_freedom(m) == 8

        # set the variables
        model_setup(m, chem_add=0)

        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_scaling(self, base_boron_removal_model):
        m = base_boron_removal_model

        # Set some scaling factors and look for 'bad' scaling
        scaling_setup(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(iscale.unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

        # check if any variables are badly scaled
        badly_scaled_var_values = {
            var.name: val
            for (var, val) in iscale.badly_scaled_var_generator(
                m, large=1e3, small=1e-3
            )
        }
        assert not badly_scaled_var_values

    @pytest.mark.component
    def test_initialization(self, base_boron_removal_model):
        m = base_boron_removal_model
        initialization_tester(m)

        # check to make sure DOF does not change
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_solve(self, base_boron_removal_model):
        m = base_boron_removal_model

        # first, check to make sure that after initialized, the scaling is still good
        badly_scaled_var_values = {
            var.name: val
            for (var, val) in iscale.badly_scaled_var_generator(
                m, large=1e3, small=1e-3
            )
        }
        assert not badly_scaled_var_values

        # run solver and check for optimal solution
        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, base_boron_removal_model):
        m = base_boron_removal_model

        assert value(
            m.fs.unit.outlet.flow_mol_phase_comp[0, "Liq", m.fs.unit.boron_name_id]
        ) == pytest.approx(1.20642e-4, rel=1e-4)
        assert value(
            m.fs.unit.outlet.flow_mol_phase_comp[0, "Liq", m.fs.unit.borate_name_id]
        ) == pytest.approx(8.03579e-5, rel=1e-4)
        assert value(m.fs.unit.outlet_pH()) == pytest.approx(9.0343, rel=1e-4)
        assert value(m.fs.unit.outlet_pOH()) == pytest.approx(4.9621, rel=1e-4)


# -----------------------------------------------------------------------------
# Start test class
class TestBoronRemoval_GenPropPack:
    @pytest.fixture(scope="class")
    def gen_boron_removal_model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        # Configuration dictionary for generic
        thermo_config = {
            "components": {
                "H2O": {
                    "type": Solvent,
                    # Define the methods used to calculate the following properties
                    "dens_mol_liq_comp": Constant,
                    "enth_mol_liq_comp": Constant,
                    "cp_mol_liq_comp": Constant,
                    "entr_mol_liq_comp": Constant,
                    # Parameter data is always associated with the methods defined above
                    "parameter_data": {
                        "mw": (18.0153, pyunits.g / pyunits.mol),
                        "dens_mol_liq_comp_coeff": (
                            55.2,
                            pyunits.kmol * pyunits.m**-3,
                        ),
                        "cp_mol_liq_comp_coeff": (
                            75.312,
                            pyunits.J / pyunits.mol / pyunits.K,
                        ),
                        "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                        "entr_mol_form_liq_comp_ref": (
                            0,
                            pyunits.J / pyunits.K / pyunits.mol,
                        ),
                    },
                    # End parameter_data
                },
                "H_+": {
                    "type": Cation,
                    "charge": 1,
                    # Define the methods used to calculate the following properties
                    "dens_mol_liq_comp": Constant,
                    "enth_mol_liq_comp": Constant,
                    "cp_mol_liq_comp": Constant,
                    "entr_mol_liq_comp": Constant,
                    # Parameter data is always associated with the methods defined above
                    "parameter_data": {
                        "mw": (1, pyunits.g / pyunits.mol),
                        "dens_mol_liq_comp_coeff": (
                            55.2,
                            pyunits.kmol * pyunits.m**-3,
                        ),
                        "cp_mol_liq_comp_coeff": (
                            75.312,
                            pyunits.J / pyunits.mol / pyunits.K,
                        ),
                        "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                        "entr_mol_form_liq_comp_ref": (
                            0,
                            pyunits.J / pyunits.K / pyunits.mol,
                        ),
                    },
                    # End parameter_data
                },
                "Na_+": {
                    "type": Cation,
                    "charge": 1,
                    # Define the methods used to calculate the following properties
                    "dens_mol_liq_comp": Constant,
                    "enth_mol_liq_comp": Constant,
                    "cp_mol_liq_comp": Constant,
                    "entr_mol_liq_comp": Constant,
                    # Parameter data is always associated with the methods defined above
                    "parameter_data": {
                        "mw": (23, pyunits.g / pyunits.mol),
                        "dens_mol_liq_comp_coeff": (
                            55.2,
                            pyunits.kmol * pyunits.m**-3,
                        ),
                        "cp_mol_liq_comp_coeff": (
                            75.312,
                            pyunits.J / pyunits.mol / pyunits.K,
                        ),
                        "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                        "entr_mol_form_liq_comp_ref": (
                            0,
                            pyunits.J / pyunits.K / pyunits.mol,
                        ),
                    },
                    # End parameter_data
                },
                "OH_-": {
                    "type": Anion,
                    "charge": -1,
                    # Define the methods used to calculate the following properties
                    "dens_mol_liq_comp": Constant,
                    "enth_mol_liq_comp": Constant,
                    "cp_mol_liq_comp": Constant,
                    "entr_mol_liq_comp": Constant,
                    # Parameter data is always associated with the methods defined above
                    "parameter_data": {
                        "mw": (17, pyunits.g / pyunits.mol),
                        "dens_mol_liq_comp_coeff": (
                            55.2,
                            pyunits.kmol * pyunits.m**-3,
                        ),
                        "cp_mol_liq_comp_coeff": (
                            75.312,
                            pyunits.J / pyunits.mol / pyunits.K,
                        ),
                        "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                        "entr_mol_form_liq_comp_ref": (
                            0,
                            pyunits.J / pyunits.K / pyunits.mol,
                        ),
                    },
                    # End parameter_data
                },
                "B[OH]4_-": {
                    "type": Anion,
                    "charge": -1,
                    # Define the methods used to calculate the following properties
                    "dens_mol_liq_comp": Constant,
                    "enth_mol_liq_comp": Constant,
                    "cp_mol_liq_comp": Constant,
                    "entr_mol_liq_comp": Constant,
                    # Parameter data is always associated with the methods defined above
                    "parameter_data": {
                        "mw": (78.83, pyunits.g / pyunits.mol),
                        "dens_mol_liq_comp_coeff": (
                            55.2,
                            pyunits.kmol * pyunits.m**-3,
                        ),
                        "cp_mol_liq_comp_coeff": (
                            75.312,
                            pyunits.J / pyunits.mol / pyunits.K,
                        ),
                        "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                        "entr_mol_form_liq_comp_ref": (
                            0,
                            pyunits.J / pyunits.K / pyunits.mol,
                        ),
                    },
                    # End parameter_data
                },
                "B[OH]3": {
                    "type": Solute,
                    "valid_phase_types": PT.aqueousPhase,
                    # Define the methods used to calculate the following properties
                    "dens_mol_liq_comp": Constant,
                    "enth_mol_liq_comp": Constant,
                    "cp_mol_liq_comp": Constant,
                    "entr_mol_liq_comp": Constant,
                    # Parameter data is always associated with the methods defined above
                    "parameter_data": {
                        "mw": (61.83, pyunits.g / pyunits.mol),
                        "dens_mol_liq_comp_coeff": (
                            55.2,
                            pyunits.kmol * pyunits.m**-3,
                        ),
                        "cp_mol_liq_comp_coeff": (
                            75.312,
                            pyunits.J / pyunits.mol / pyunits.K,
                        ),
                        "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                        "entr_mol_form_liq_comp_ref": (
                            0,
                            pyunits.J / pyunits.K / pyunits.mol,
                        ),
                    },
                    # End parameter_data
                },
            },
            # End Component list
            "phases": {
                "Liq": {"type": AqueousPhase, "equation_of_state": Ideal},
            },
            "state_definition": FpcTP,
            "state_bounds": {
                "temperature": (273.15, 300, 650),
                "pressure": (5e4, 1e5, 1e6),
            },
            "pressure_ref": 1e5,
            "temperature_ref": 300,
            "base_units": {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            },
        }
        # End thermo_config definition

        # attach prop pack to flowsheet
        m.fs.properties = GenericParameterBlock(**thermo_config)

        map = {
            "boron_name": "B[OH]3",  # [is required]
            "borate_name": "B[OH]4_-",  # [is required]
            "proton_name": "H_+",  # [is optional]
            "hydroxide_name": "OH_-",  # [is optional]
            "caustic_additive": {
                "cation_name": "Na_+",  # [is required]
                "mw_additive": (40, pyunits.g / pyunits.mol),  # [is required]
                "moles_cation_per_additive": 1,  # [is required]
            },
        }
        m.fs.unit = BoronRemoval(
            property_package=m.fs.properties, chemical_mapping_data=map
        )

        return m

    @pytest.mark.unit
    def test_build_model(self, gen_boron_removal_model):
        m = gen_boron_removal_model

        assert isinstance(m.fs.unit.ion_charge, Param)
        assert len(m.fs.unit.ion_charge) == 4

        assert hasattr(m.fs.unit, "boron_name_id")
        assert hasattr(m.fs.unit, "borate_name_id")
        assert hasattr(m.fs.unit, "proton_name_id")
        assert hasattr(m.fs.unit, "hydroxide_name_id")
        assert hasattr(m.fs.unit, "cation_name_id")

        assert hasattr(m.fs.unit, "caustic_chem_name")

        assert hasattr(m.fs.unit, "control_volume")

        assert isinstance(m.fs.unit.caustic_mw, Param)
        assert isinstance(m.fs.unit.additive_molar_ratio, Param)
        assert isinstance(m.fs.unit.caustic_dose_rate, Var)

        assert isinstance(m.fs.unit.Kw_0, Param)
        assert isinstance(m.fs.unit.dH_w, Param)
        assert isinstance(m.fs.unit.Ka_0, Param)
        assert isinstance(m.fs.unit.dH_a, Param)

        assert isinstance(m.fs.unit.conc_mol_H, Var)
        assert isinstance(m.fs.unit.conc_mol_OH, Var)
        assert isinstance(m.fs.unit.conc_mol_Boron, Var)
        assert isinstance(m.fs.unit.conc_mol_Borate, Var)

        assert isinstance(m.fs.unit.eq_mass_transfer_term, Constraint)
        assert isinstance(m.fs.unit.eq_electroneutrality, Constraint)
        assert isinstance(m.fs.unit.eq_total_boron, Constraint)
        assert isinstance(m.fs.unit.eq_water_dissociation, Constraint)
        assert isinstance(m.fs.unit.eq_boron_dissociation, Constraint)

    @pytest.mark.unit
    def test_stats(self, gen_boron_removal_model):
        m = gen_boron_removal_model

        # Check to make sure we have the correct DOF and
        #   check to make sure the units are correct
        assert_units_consistent(m)
        assert degrees_of_freedom(m) == 10

        # set the variables
        model_setup(m, chem_add=1.8e-6)

        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_scaling(self, gen_boron_removal_model):
        m = gen_boron_removal_model

        # Set some scaling factors and look for 'bad' scaling
        scaling_setup(m)

        ## TODO: Revisit scaling for the generic package
        #   The unscaled vars come from the 'true-to-apparent'
        #   species vars in the generic package (which we don't
        #   ever use). May need to revisit how to scale
        # check that all variables have scaling factors
        unscaled_var_list = list(iscale.unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 8

    @pytest.mark.component
    def test_initialization(self, gen_boron_removal_model):
        m = gen_boron_removal_model
        initialization_tester(m)

        # check to make sure DOF does not change
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_solve(self, gen_boron_removal_model):
        m = gen_boron_removal_model

        # run solver and check for optimal solution
        results = solver.solve(m)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, gen_boron_removal_model):
        m = gen_boron_removal_model

        assert value(
            m.fs.unit.outlet.flow_mol_phase_comp[0, "Liq", m.fs.unit.boron_name_id]
        ) == pytest.approx(1.20783e-4, rel=1e-4)
        assert value(
            m.fs.unit.outlet.flow_mol_phase_comp[0, "Liq", m.fs.unit.borate_name_id]
        ) == pytest.approx(8.02173e-5, rel=1e-4)
        assert value(m.fs.unit.outlet_pH()) == pytest.approx(9.0330, rel=1e-4)
        assert value(m.fs.unit.outlet_pOH()) == pytest.approx(4.9633, rel=1e-4)


# -----------------------------------------------------------------------------
# Start test class with bad config
class TestBoronRemoval_BadConfigs:
    @pytest.fixture(scope="class")
    def boron_removal_bad_configs(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        # create dict to define ions (the prop pack requires this)
        ion_dict = {
            "solute_list": ["B[OH]3", "B[OH]4_-", "H_+", "OH_-", "Na_+"],
            "mw_data": {"H2O": 18e-3, "B[OH]3": 61.83e-3, "B[OH]4_-": 78.83e-3},
            "charge": {
                "B[OH]4_-": -1,
            },
        }

        # attach prop pack to flowsheet
        m.fs.properties = DSPMDEParameterBlock(**ion_dict)

        return m

    @pytest.mark.unit
    def test_build_failures(self, boron_removal_bad_configs):
        m = boron_removal_bad_configs

        # Invalid name of boron
        map = {
            "boron_name": "Boron",  # [is required]
            "borate_name": "B[OH]4_-",  # [is required]
            "caustic_additive": {
                "cation_name": "Na_+",  # [is optional]
                "mw_additive": (40, pyunits.g / pyunits.mol),  # [is required]
                "moles_cation_per_additive": 1,  # [is required]
            },
        }
        error_msg = (
            "Given 'boron_name' {Boron} does not match any species "
            "name from the property package "
        )
        with pytest.raises(ConfigurationError, match=re.escape(error_msg)):
            m.fs.unit = BoronRemoval(
                property_package=m.fs.properties, chemical_mapping_data=map
            )

        # Invalid name of borate
        map = {
            "boron_name": "B[OH]3",  # [is required]
            "borate_name": "Borate",  # [is required]
            "caustic_additive": {
                "cation_name": "Na_+",  # [is optional]
                "mw_additive": (40, pyunits.g / pyunits.mol),  # [is required]
                "moles_cation_per_additive": 1,  # [is required]
            },
        }
        error_msg = (
            "Given 'borate_name' {Borate} does not match any species "
            "name from the property package "
        )
        with pytest.raises(ConfigurationError, match=re.escape(error_msg)):
            m.fs.unit = BoronRemoval(
                property_package=m.fs.properties, chemical_mapping_data=map
            )

        # Empty dict
        map = {}
        error_msg = "Did not provide a 'dict' for 'chemical_mapping_data' "
        with pytest.raises(ConfigurationError, match=re.escape(error_msg)):
            m.fs.unit = BoronRemoval(
                property_package=m.fs.properties, chemical_mapping_data=map
            )

        # Invalid name of hydroxide
        map = {
            "boron_name": "B[OH]3",  # [is required]
            "borate_name": "B[OH]4_-",  # [is required]
            "hydroxide_name": "Raccoon_City",
            "caustic_additive": {
                "cation_name": "Na_+",  # [is optional]
                "mw_additive": (40, pyunits.g / pyunits.mol),  # [is required]
                "moles_cation_per_additive": 1,  # [is required]
            },
        }
        error_msg = (
            "Given 'hydroxide_name' {Raccoon_City} does not match any species "
            "name from the property package "
        )
        with pytest.raises(ConfigurationError, match=re.escape(error_msg)):
            m.fs.unit = BoronRemoval(
                property_package=m.fs.properties, chemical_mapping_data=map
            )

        # Invalid name of protons
        map = {
            "boron_name": "B[OH]3",  # [is required]
            "borate_name": "B[OH]4_-",  # [is required]
            "proton_name": "Pulled_Pork",
            "caustic_additive": {
                "cation_name": "Na_+",  # [is optional]
                "mw_additive": (40, pyunits.g / pyunits.mol),  # [is required]
                "moles_cation_per_additive": 1,  # [is required]
            },
        }
        error_msg = (
            "Given 'proton_name' {Pulled_Pork} does not match any species "
            "name from the property package "
        )
        with pytest.raises(ConfigurationError, match=re.escape(error_msg)):
            m.fs.unit = BoronRemoval(
                property_package=m.fs.properties, chemical_mapping_data=map
            )

        # Invalid name of cation
        map = {
            "boron_name": "B[OH]3",  # [is required]
            "borate_name": "B[OH]4_-",  # [is required]
            "caustic_additive": {
                "cation_name": "Im_Batman",  # [is optional]
                "mw_additive": (40, pyunits.g / pyunits.mol),  # [is required]
                "moles_cation_per_additive": 1,  # [is required]
            },
        }
        error_msg = (
            "Given 'cation_name' {Im_Batman} does not match any species "
            "name from the property package "
        )
        with pytest.raises(ConfigurationError, match=re.escape(error_msg)):
            m.fs.unit = BoronRemoval(
                property_package=m.fs.properties, chemical_mapping_data=map
            )

        # Missing information (borate_name)
        map = {
            "boron_name": "B[OH]3",  # [is required]
            "caustic_additive": {
                "mw_additive": (40, pyunits.g / pyunits.mol),  # [is required]
                "moles_cation_per_additive": 1,  # [is required]
            },
        }
        error_msg = "Missing some required information in 'chemical_mapping_data' "
        with pytest.raises(ConfigurationError, match=re.escape(error_msg)):
            m.fs.unit = BoronRemoval(
                property_package=m.fs.properties, chemical_mapping_data=map
            )

        # Missing information (mw_additive)
        map = {
            "boron_name": "B[OH]3",  # [is required]
            "borate_name": "B[OH]4_-",  # [is required]
            "caustic_additive": {
                "moles_cation_per_additive": 1,  # [is required]
            },
        }
        error_msg = "Missing some required information in 'chemical_mapping_data' "
        with pytest.raises(ConfigurationError, match=re.escape(error_msg)):
            m.fs.unit = BoronRemoval(
                property_package=m.fs.properties, chemical_mapping_data=map
            )

        # Improper data type (mw_additive must be a tuple with value and units)
        map = {
            "boron_name": "B[OH]3",  # [is required]
            "borate_name": "B[OH]4_-",  # [is required]
            "caustic_additive": {
                "cation_name": "Na_+",  # [is required]
                "mw_additive": 40,  # [is required]
                "moles_cation_per_additive": 1,  # [is required]
            },
        }
        error_msg = "Did not provide a tuple for 'mw_additive' "
        with pytest.raises(ConfigurationError, match=re.escape(error_msg)):
            m.fs.unit = BoronRemoval(
                property_package=m.fs.properties, chemical_mapping_data=map
            )


# -----------------------------------------------------------------------------
# Start test class with bad config
class TestBoronRemoval_BadConfigs_Generic:
    @pytest.fixture(scope="class")
    def boron_removal_bad_configs_gen(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        # Configuration dictionary for generic
        thermo_config = {
            "components": {
                "H2O": {
                    "type": Solvent,
                    # Define the methods used to calculate the following properties
                    "dens_mol_liq_comp": Constant,
                    "enth_mol_liq_comp": Constant,
                    "cp_mol_liq_comp": Constant,
                    "entr_mol_liq_comp": Constant,
                    # Parameter data is always associated with the methods defined above
                    "parameter_data": {
                        "mw": (18.0153, pyunits.g / pyunits.mol),
                        "dens_mol_liq_comp_coeff": (
                            55.2,
                            pyunits.kmol * pyunits.m**-3,
                        ),
                        "cp_mol_liq_comp_coeff": (
                            75.312,
                            pyunits.J / pyunits.mol / pyunits.K,
                        ),
                        "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                        "entr_mol_form_liq_comp_ref": (
                            0,
                            pyunits.J / pyunits.K / pyunits.mol,
                        ),
                    },
                    # End parameter_data
                },
                "H_+": {
                    "type": Cation,
                    "charge": 1,
                    # Define the methods used to calculate the following properties
                    "dens_mol_liq_comp": Constant,
                    "enth_mol_liq_comp": Constant,
                    "cp_mol_liq_comp": Constant,
                    "entr_mol_liq_comp": Constant,
                    # Parameter data is always associated with the methods defined above
                    "parameter_data": {
                        "mw": (1, pyunits.g / pyunits.mol),
                        "dens_mol_liq_comp_coeff": (
                            55.2,
                            pyunits.kmol * pyunits.m**-3,
                        ),
                        "cp_mol_liq_comp_coeff": (
                            75.312,
                            pyunits.J / pyunits.mol / pyunits.K,
                        ),
                        "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                        "entr_mol_form_liq_comp_ref": (
                            0,
                            pyunits.J / pyunits.K / pyunits.mol,
                        ),
                    },
                    # End parameter_data
                },
                "Na_+": {
                    "type": Cation,
                    "charge": 1,
                    # Define the methods used to calculate the following properties
                    "dens_mol_liq_comp": Constant,
                    "enth_mol_liq_comp": Constant,
                    "cp_mol_liq_comp": Constant,
                    "entr_mol_liq_comp": Constant,
                    # Parameter data is always associated with the methods defined above
                    "parameter_data": {
                        "mw": (23, pyunits.g / pyunits.mol),
                        "dens_mol_liq_comp_coeff": (
                            55.2,
                            pyunits.kmol * pyunits.m**-3,
                        ),
                        "cp_mol_liq_comp_coeff": (
                            75.312,
                            pyunits.J / pyunits.mol / pyunits.K,
                        ),
                        "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                        "entr_mol_form_liq_comp_ref": (
                            0,
                            pyunits.J / pyunits.K / pyunits.mol,
                        ),
                    },
                    # End parameter_data
                },
                "OH_-": {
                    "type": Anion,
                    "charge": -1,
                    # Define the methods used to calculate the following properties
                    "dens_mol_liq_comp": Constant,
                    "enth_mol_liq_comp": Constant,
                    "cp_mol_liq_comp": Constant,
                    "entr_mol_liq_comp": Constant,
                    # Parameter data is always associated with the methods defined above
                    "parameter_data": {
                        "mw": (17, pyunits.g / pyunits.mol),
                        "dens_mol_liq_comp_coeff": (
                            55.2,
                            pyunits.kmol * pyunits.m**-3,
                        ),
                        "cp_mol_liq_comp_coeff": (
                            75.312,
                            pyunits.J / pyunits.mol / pyunits.K,
                        ),
                        "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                        "entr_mol_form_liq_comp_ref": (
                            0,
                            pyunits.J / pyunits.K / pyunits.mol,
                        ),
                    },
                    # End parameter_data
                },
                "B[OH]4_-": {
                    "type": Anion,
                    "charge": -1,
                    # Define the methods used to calculate the following properties
                    "dens_mol_liq_comp": Constant,
                    "enth_mol_liq_comp": Constant,
                    "cp_mol_liq_comp": Constant,
                    "entr_mol_liq_comp": Constant,
                    # Parameter data is always associated with the methods defined above
                    "parameter_data": {
                        "mw": (78.83, pyunits.g / pyunits.mol),
                        "dens_mol_liq_comp_coeff": (
                            55.2,
                            pyunits.kmol * pyunits.m**-3,
                        ),
                        "cp_mol_liq_comp_coeff": (
                            75.312,
                            pyunits.J / pyunits.mol / pyunits.K,
                        ),
                        "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                        "entr_mol_form_liq_comp_ref": (
                            0,
                            pyunits.J / pyunits.K / pyunits.mol,
                        ),
                    },
                    # End parameter_data
                },
                "B[OH]3": {
                    "type": Solute,
                    "valid_phase_types": PT.aqueousPhase,
                    # Define the methods used to calculate the following properties
                    "dens_mol_liq_comp": Constant,
                    "enth_mol_liq_comp": Constant,
                    "cp_mol_liq_comp": Constant,
                    "entr_mol_liq_comp": Constant,
                    # Parameter data is always associated with the methods defined above
                    "parameter_data": {
                        "mw": (61.83, pyunits.g / pyunits.mol),
                        "dens_mol_liq_comp_coeff": (
                            55.2,
                            pyunits.kmol * pyunits.m**-3,
                        ),
                        "cp_mol_liq_comp_coeff": (
                            75.312,
                            pyunits.J / pyunits.mol / pyunits.K,
                        ),
                        "enth_mol_form_liq_comp_ref": (0, pyunits.kJ / pyunits.mol),
                        "entr_mol_form_liq_comp_ref": (
                            0,
                            pyunits.J / pyunits.K / pyunits.mol,
                        ),
                    },
                    # End parameter_data
                },
            },
            # End Component list
            "phases": {
                "Liq": {"type": AqueousPhase, "equation_of_state": Ideal},
            },
            "state_definition": FpcTP,
            "state_bounds": {
                "temperature": (273.15, 300, 650),
                "pressure": (5e4, 1e5, 1e6),
            },
            "pressure_ref": 1e5,
            "temperature_ref": 300,
            "base_units": {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            },
            "inherent_reactions": {
                "H2O_Kw": {
                    "stoichiometry": {
                        ("Liq", "H2O"): -1,
                        ("Liq", "H_+"): 1,
                        ("Liq", "OH_-"): 1,
                    },
                    "heat_of_reaction": constant_dh_rxn,
                    "equilibrium_constant": van_t_hoff,
                    "equilibrium_form": log_power_law_equil,
                    "concentration_form": ConcentrationForm.moleFraction,
                    "parameter_data": {
                        "dh_rxn_ref": (55.830, pyunits.J / pyunits.mol),
                        "k_eq_ref": (10**-14 / 55.2 / 55.2, pyunits.dimensionless),
                        "T_eq_ref": (298, pyunits.K),
                        # By default, reaction orders follow stoichiometry
                        #    manually set reaction order here to override
                        "reaction_order": {
                            ("Liq", "H2O"): 0,
                            ("Liq", "H_+"): 1,
                            ("Liq", "OH_-"): 1,
                        },
                    }
                    # End parameter_data
                }
                # End R1
            },
        }
        # End thermo_config definition

        # attach prop pack to flowsheet
        m.fs.properties = GenericParameterBlock(**thermo_config)

        return m

    @pytest.mark.unit
    def test_build_failures(self, boron_removal_bad_configs_gen):
        m = boron_removal_bad_configs_gen

        # Invalid name of boron
        map = {
            "boron_name": "B[OH]3",  # [is required]
            "borate_name": "B[OH]4_-",  # [is required]
            "caustic_additive": {
                "cation_name": "Na_+",  # [is required]
                "mw_additive": (40, pyunits.g / pyunits.mol),  # [is required]
                "moles_cation_per_additive": 1,  # [is required]
            },
        }
        error_msg = "Property Package CANNOT contain 'inherent_reactions' "
        with pytest.raises(ConfigurationError, match=re.escape(error_msg)):
            m.fs.unit = BoronRemoval(
                property_package=m.fs.properties, chemical_mapping_data=map
            )
