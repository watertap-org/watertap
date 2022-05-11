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
from idaes.core.components import Solvent, Solute, Cation, Anion
from idaes.core.phases import PhaseType as PT

# Imports from idaes generic models
from idaes.generic_models.properties.core.pure.ConstantProperties import Constant
from idaes.generic_models.properties.core.state_definitions import FpcTP
from idaes.generic_models.properties.core.eos.ideal import Ideal

# Import the idaes objects for Generic Properties and Reactions
from idaes.generic_models.properties.core.generic.generic_property import (
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
from idaes.core.util import get_solver
import idaes.logger as idaeslog
import re

__author__ = "Austin Ladshaw"

solver = get_solver()

# Helper function for multiple test setup
def model_setup(m, chem_add = 10,
                    state={"H2O": 100, "H_+": 1e-7, "OH_-": 1e-7,
                          "B[OH]3": 2e-4, "B[OH]4_-": 1e-6,
                          "Na_+": 1e-4, "HCO3_-": 1e-4}):

    m.fs.unit.inlet.pressure.fix(101325)
    m.fs.unit.inlet.temperature.fix(298.15)
    for j in state:
        idx = (0, "Liq", j)
        if idx in m.fs.unit.inlet.flow_mol_phase_comp:
            m.fs.unit.inlet.flow_mol_phase_comp[idx].fix(state[j])
    m.fs.unit.caustic_dose.fix(chem_add)

# Helper function to automate scaling
def scaling_setup(m, state={"H2O": 100, "H_+": 5e-5, "OH_-": 5e-5,
                          "B[OH]3": 2e-4, "B[OH]4_-": 1e-6,
                          "Na_+": 1e-4, "HCO3_-": 1e-4}):
    # Set some scaling factors and look for 'bad' scaling
    for j in state:
        idx = (0, "Liq", j)
        if idx in m.fs.unit.inlet.flow_mol_phase_comp:
            m.fs.properties.set_default_scaling(
                "flow_mol_phase_comp", 1/state[j], index=("Liq", j)
            )

    iscale.calculate_scaling_factors(m.fs)

# -----------------------------------------------------------------------------
# Start test class
class TestBoronRemoval_IonPropPack_Min:
    @pytest.fixture(scope="class")
    def min_boron_removal_model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})
        # create dict to define ions (the prop pack requires this)
        ion_dict = {
            "solute_list": ["B[OH]3", "B[OH]4_-"],
            "mw_data": {"H2O": 18e-3, "B[OH]3": 61.83e-3, "B[OH]4_-": 78.83e-3},
            "charge": {"B[OH]3": 0, "B[OH]4_-": -1,},
        }

        # attach prop pack to flowsheet
        m.fs.properties = DSPMDEParameterBlock(default=ion_dict)

        map = {'boron_name': 'B[OH]3',
                'borate_name': 'B[OH]4_-',
                'caustic_additive':
                    {
                        'additive_name': 'NaOH',
                        'mw_additive': (40, pyunits.g/pyunits.mol),
                        'charge_additive': 1,
                    },
        }
        m.fs.unit = BoronRemoval(
            default={
                "property_package": m.fs.properties,
                "chemical_mapping_data": map,
            }
        )

        return m

    @pytest.mark.unit
    def test_build_model(self, min_boron_removal_model):
        m = min_boron_removal_model

        assert isinstance(m.fs.unit.ion_charge, Param)
        assert len(m.fs.unit.ion_charge) == 1

        assert hasattr(m.fs.unit, 'boron_name_id')
        assert hasattr(m.fs.unit, 'borate_name_id')
        assert hasattr(m.fs.unit, 'proton_name_id')
        assert hasattr(m.fs.unit, 'hydroxide_name_id')
        assert hasattr(m.fs.unit, 'cation_name_id')

        assert hasattr(m.fs.unit, 'caustic_chem_name')

        assert hasattr(m.fs.unit, 'control_volume')

        assert isinstance(m.fs.unit.caustic_mw, Param)
        assert isinstance(m.fs.unit.caustic_cation_charge, Param)
        assert isinstance(m.fs.unit.caustic_dose, Var)

        assert isinstance(m.fs.unit.Kw_0, Param)
        assert isinstance(m.fs.unit.dH_w, Param)
        assert isinstance(m.fs.unit.Ka_0, Param)
        assert isinstance(m.fs.unit.dH_a, Param)

        assert isinstance(m.fs.unit.mol_H, Var)
        assert isinstance(m.fs.unit.mol_OH, Var)
        assert isinstance(m.fs.unit.mol_Boron, Var)
        assert isinstance(m.fs.unit.mol_Borate, Var)

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
        assert degrees_of_freedom(m) == 6

        # set the variables
        model_setup(m)


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
        results = solver.solve(m)
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
        assert value(
            m.fs.unit.outlet_pH()
        ) == pytest.approx(10.171, rel=1e-4)
        assert value(
            m.fs.unit.outlet_pOH()
        ) == pytest.approx(3.8257, rel=1e-4)

# -----------------------------------------------------------------------------
# Start test class with bad config
class TestBoronRemoval_BadConfigs:
    @pytest.fixture(scope="class")
    def boron_removal_bad_configs(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        return m
