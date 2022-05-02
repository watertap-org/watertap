###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
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
import pytest
from pyomo.environ import (
    ConcreteModel,
    SolverFactory,
    TerminationCondition,
    value,
    Constraint,
    Var,
    Objective,
    Expression,
    Param,
    Set,
    assert_optimal_termination,
)
from pyomo.environ import units as pyunits
from pyomo.util.check_units import (
    assert_units_consistent,
    assert_units_equivalent,
    check_units_equivalent,
)
import pyomo.util.infeasible as infeas
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom, number_variables
import idaes.core.util.model_statistics as stats
from idaes.core.util.constants import Constants
import idaes.core.util.scaling as iscale
from idaes.core.util.testing import initialization_tester
import idaes.logger as idaeslog

from watertap.property_models.ion_DSPMDE_prop_pack import DSPMDEParameterBlock
from watertap.unit_models.electrodialysis_1D import Electrodialysis1D

from idaes.generic_models.properties.core.generic.generic_property import (
    GenericParameterBlock,
)

from idaes.core import AqueousPhase
from idaes.core.components import Solvent, Solute, Cation, Anion
from idaes.core.phases import PhaseType as PT

from idaes.generic_models.properties.core.pure.ConstantProperties import Constant
from idaes.generic_models.properties.core.state_definitions import FpcTP
from idaes.generic_models.properties.core.eos.ideal import Ideal

from idaes.core.util import get_solver

solver = get_solver()

__author__ = "Austin Ladshaw"

_log = idaeslog.getLogger(__name__)


def build_dspmde_model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    # create dict to define ions (the prop pack of Adam requires this)
    ion_dict = {
        "solute_list": ["Na_+", "Cl_-", "NaCl"],
        "diffusivity_data": {
            ("Liq", "Na_+"): 1.33e-9,
            ("Liq", "Cl_-"): 2.03e-9,
            ("Liq", "NaCl"): 1.70e-9,
        },
        "mw_data": {"H2O": 18e-3, "Na_+": 23e-3, "Cl_-": 35e-3, "NaCl": 58e-3},
        "stokes_radius_data": {"Na_+": 0.184e-9, "Cl_-": 0.121e-9, "NaCl": 0.305e-9},
        "charge": {"Na_+": 1, "Cl_-": -1, "NaCl": 0},
    }
    # attach prop pack to flowsheet
    m.fs.properties = DSPMDEParameterBlock(default=ion_dict)

    # build the unit model, pass prop pack to the model
    m.fs.unit = Electrodialysis1D(default={"property_package": m.fs.properties})

    return m


def build_generic_model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
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
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
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
                    "mw": (22.989769, pyunits.g / pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
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
            "Cl_-": {
                "type": Anion,
                "charge": -1,
                # Define the methods used to calculate the following properties
                "dens_mol_liq_comp": Constant,
                "enth_mol_liq_comp": Constant,
                "cp_mol_liq_comp": Constant,
                "entr_mol_liq_comp": Constant,
                # Parameter data is always associated with the methods defined above
                "parameter_data": {
                    "mw": (35.453, pyunits.g / pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
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
            "NaCl": {
                "type": Solute,
                "valid_phase_types": PT.aqueousPhase,
                # Define the methods used to calculate the following properties
                "dens_mol_liq_comp": Constant,
                "enth_mol_liq_comp": Constant,
                "cp_mol_liq_comp": Constant,
                "entr_mol_liq_comp": Constant,
                # Parameter data is always associated with the methods defined above
                "parameter_data": {
                    "mw": (58.442, pyunits.g / pyunits.mol),
                    "dens_mol_liq_comp_coeff": (55.2, pyunits.kmol * pyunits.m**-3),
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
    m.fs.properties = GenericParameterBlock(default=thermo_config)

    # build the unit model, pass prop pack to the model
    m.fs.unit = Electrodialysis1D(default={"property_package": m.fs.properties})

    return m


def fix_inlets_and_vars_no_current(m):
    # specify the feed for each inlet stream
    m.fs.unit.inlet_dilute.pressure.fix(101325)
    m.fs.unit.inlet_dilute.temperature.fix(298.15)
    m.fs.unit.inlet_dilute.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2)
    m.fs.unit.inlet_dilute.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(0.2)
    m.fs.unit.inlet_dilute.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(0.2)
    m.fs.unit.inlet_dilute.flow_mol_phase_comp[0, "Liq", "NaCl"].fix(0.002)

    m.fs.unit.inlet_concentrate.pressure.fix(101325)
    m.fs.unit.inlet_concentrate.temperature.fix(298.15)
    m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(1)
    m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(0.15)
    m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(0.15)
    m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "NaCl"].fix(0.0015)

    m.fs.unit.cell_length.fix(0.5)
    m.fs.unit.cell_width.fix(0.1)
    m.fs.unit.membrane_thickness["cem"].fix()
    m.fs.unit.membrane_thickness["aem"].fix()

    m.fs.unit.ion_diffusivity_membrane.fix()
    m.fs.unit.water_permeability_membrane.fix()


def scale_model(m):
    m.fs.properties.set_default_scaling("flow_mol_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e1, index=("Liq", "Na_+")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e1, index=("Liq", "Cl_-")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e3, index=("Liq", "NaCl")
    )

    iscale.calculate_scaling_factors(m)


def check_scaling(m):
    unscaled_constraint_list = list(iscale.unscaled_constraints_generator(m))
    if len(unscaled_constraint_list) > 0:
        print("List of unscaled constraints")
        print("----------------------------")
        for j in unscaled_constraint_list:
            print(j)
        print()

    # check that all variables have scaling factors
    unscaled_var_list = list(iscale.unscaled_variables_generator(m))
    if len(unscaled_var_list) > 0:
        print("List of unscaled variables")
        print("--------------------------")
        for j in unscaled_var_list:
            print(j)
        print()

    # check if any variables are badly scaled
    badly_scaled_var_values = {
        var.name: val
        for (var, val) in iscale.badly_scaled_var_generator(m, large=1e3, small=1e-3)
    }
    if len(badly_scaled_var_values) > 0:
        print("List of poorly scaled variables")
        print("-------------------------------")
        for j in badly_scaled_var_values:
            print(str(j) + "\t" + str(badly_scaled_var_values[j]))
        print()

    return (unscaled_constraint_list, unscaled_var_list, badly_scaled_var_values)


# -----------------------------------------------------------------------------
# Start test class
class TestElectrodialysis1D_withDSPMDEProps_noCurrent:
    @pytest.fixture(scope="class")
    def elec1d_dspmde(self):
        model = build_dspmde_model()
        return model

    @pytest.mark.unit
    def test_build_model(self, elec1d_dspmde):
        model = elec1d_dspmde

        assert len(model.fs.unit.config.property_package.component_list) == 4
        assert len(model.fs.unit.config.property_package.cation_set) == 1
        assert len(model.fs.unit.config.property_package.anion_set) == 1

        assert isinstance(model.fs.unit.ion_charge, Param)
        assert isinstance(model.fs.unit.membrane_set, Set)
        assert isinstance(model.fs.unit.cell_width, Var)
        assert isinstance(model.fs.unit.membrane_thickness, Var)
        assert isinstance(model.fs.unit.ion_diffusivity_membrane, Var)
        assert isinstance(model.fs.unit.water_permeability_membrane, Var)

        assert isinstance(model.fs.unit.dilute_side.length, Var)
        assert isinstance(model.fs.unit.concentrate_side.length, Var)

        assert isinstance(model.fs.unit.eq_equal_length, Constraint)

        assert isinstance(model.fs.unit.nonelec_flux, Var)
        assert isinstance(model.fs.unit.eq_nonelec_flux, Constraint)

        assert isinstance(model.fs.unit.elec_flux, Var)
        assert isinstance(model.fs.unit.eq_elec_flux, Constraint)

        assert isinstance(model.fs.unit.eq_mass_transfer_term_dilute, Constraint)
        assert isinstance(model.fs.unit.eq_mass_transfer_term_concentrate, Constraint)

    @pytest.mark.unit
    def test_stats(self, elec1d_dspmde):
        model = elec1d_dspmde

        # Check to make sure we have the correct DOF and
        #   check to make sure the units are correct
        assert_units_consistent(model)
        assert degrees_of_freedom(model) == 22

        fix_inlets_and_vars_no_current(model)
        assert degrees_of_freedom(model) == 0

    @pytest.mark.unit
    def test_scaling(self, elec1d_dspmde):
        model = elec1d_dspmde
        scale_model(model)
        (
            unscaled_constraint_list,
            unscaled_var_list,
            badly_scaled_var_values,
        ) = check_scaling(model)

        assert len(unscaled_constraint_list) == 0
        assert len(unscaled_var_list) == 0
        assert len(badly_scaled_var_values) == 0

    @pytest.mark.component
    def test_initialization(self, elec1d_dspmde):
        model = elec1d_dspmde
        model.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)

        # check to make sure DOF does not change
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_solve(self, elec1d_dspmde):
        model = elec1d_dspmde

        (
            unscaled_constraint_list,
            unscaled_var_list,
            badly_scaled_var_values,
        ) = check_scaling(model)

        assert len(unscaled_constraint_list) == 0
        assert len(unscaled_var_list) == 0
        assert len(badly_scaled_var_values) == 0

        # run solver and check for optimal solution
        results = solver.solve(model)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, elec1d_dspmde):
        model = elec1d_dspmde

        table = model.fs.unit._get_stream_table_contents()

        k1 = "Dilute Side Outlet"
        k2 = "flow_mol_phase_comp ('Liq', 'H2O')"
        assert value(table[k1][k2]) == pytest.approx(1.98393, rel=1e-4)

        k1 = "Dilute Side Outlet"
        k2 = "flow_mol_phase_comp ('Liq', 'NaCl')"
        assert value(table[k1][k2]) == pytest.approx(0.002, rel=1e-4)

        k1 = "Dilute Side Outlet"
        k2 = "flow_mol_phase_comp ('Liq', 'Na_+')"
        assert value(table[k1][k2]) == pytest.approx(0.20013, rel=1e-4)

        k1 = "Dilute Side Outlet"
        k2 = "flow_mol_phase_comp ('Liq', 'Cl_-')"
        assert value(table[k1][k2]) == pytest.approx(0.20013, rel=1e-4)

        k1 = "Concentrate Side Outlet"
        k2 = "flow_mol_phase_comp ('Liq', 'H2O')"
        assert value(table[k1][k2]) == pytest.approx(1.016069, rel=1e-4)

        k1 = "Concentrate Side Outlet"
        k2 = "flow_mol_phase_comp ('Liq', 'NaCl')"
        assert value(table[k1][k2]) == pytest.approx(0.0015, rel=1e-4)

        k1 = "Concentrate Side Outlet"
        k2 = "flow_mol_phase_comp ('Liq', 'Na_+')"
        assert value(table[k1][k2]) == pytest.approx(0.14986, rel=1e-4)

        k1 = "Concentrate Side Outlet"
        k2 = "flow_mol_phase_comp ('Liq', 'Cl_-')"
        assert value(table[k1][k2]) == pytest.approx(0.14986, rel=1e-4)


# -----------------------------------------------------------------------------
# Start test class
class TestElectrodialysis1D_withGenericProps_noCurrent:
    @pytest.fixture(scope="class")
    def elec1d_generic(self):
        model = build_generic_model()
        return model

    @pytest.mark.unit
    def test_build_model(self, elec1d_generic):
        model = elec1d_generic

        assert len(model.fs.unit.config.property_package.component_list) == 4
        assert len(model.fs.unit.config.property_package.cation_set) == 1
        assert len(model.fs.unit.config.property_package.anion_set) == 1

        assert isinstance(model.fs.unit.ion_charge, Param)
        assert isinstance(model.fs.unit.membrane_set, Set)
        assert isinstance(model.fs.unit.cell_width, Var)
        assert isinstance(model.fs.unit.membrane_thickness, Var)
        assert isinstance(model.fs.unit.ion_diffusivity_membrane, Var)
        assert isinstance(model.fs.unit.water_permeability_membrane, Var)

        assert isinstance(model.fs.unit.dilute_side.length, Var)
        assert isinstance(model.fs.unit.concentrate_side.length, Var)

        assert isinstance(model.fs.unit.eq_equal_length, Constraint)

        assert isinstance(model.fs.unit.nonelec_flux, Var)
        assert isinstance(model.fs.unit.eq_nonelec_flux, Constraint)

        assert isinstance(model.fs.unit.elec_flux, Var)
        assert isinstance(model.fs.unit.eq_elec_flux, Constraint)

        assert isinstance(model.fs.unit.eq_mass_transfer_term_dilute, Constraint)
        assert isinstance(model.fs.unit.eq_mass_transfer_term_concentrate, Constraint)

    @pytest.mark.unit
    def test_stats(self, elec1d_generic):
        model = elec1d_generic

        # Check to make sure we have the correct DOF and
        #   check to make sure the units are correct
        assert_units_consistent(model)
        assert degrees_of_freedom(model) == 22

        fix_inlets_and_vars_no_current(model)
        assert degrees_of_freedom(model) == 0

    @pytest.mark.unit
    def test_scaling(self, elec1d_generic):
        model = elec1d_generic
        scale_model(model)
        (
            unscaled_constraint_list,
            unscaled_var_list,
            badly_scaled_var_values,
        ) = check_scaling(model)

        assert len(unscaled_constraint_list) == 0
        assert len(unscaled_var_list) == 0
        assert len(badly_scaled_var_values) == 0

    @pytest.mark.component
    def test_initialization(self, elec1d_generic):
        model = elec1d_generic
        model.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)

        # check to make sure DOF does not change
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_solve(self, elec1d_generic):
        model = elec1d_generic

        (
            unscaled_constraint_list,
            unscaled_var_list,
            badly_scaled_var_values,
        ) = check_scaling(model)

        assert len(unscaled_constraint_list) == 0
        assert len(unscaled_var_list) == 0
        assert len(badly_scaled_var_values) == 0

        # run solver and check for optimal solution
        results = solver.solve(model)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, elec1d_generic):
        model = elec1d_generic

        table = model.fs.unit._get_stream_table_contents()

        k1 = "Dilute Side Outlet"
        k2 = "Molar Flowrate ('Liq', 'H2O')"
        assert value(table[k1][k2]) == pytest.approx(1.97252, rel=1e-4)

        k1 = "Dilute Side Outlet"
        k2 = "Molar Flowrate ('Liq', 'NaCl')"
        assert value(table[k1][k2]) == pytest.approx(0.002, rel=1e-4)

        k1 = "Dilute Side Outlet"
        k2 = "Molar Flowrate ('Liq', 'Na_+')"
        assert value(table[k1][k2]) == pytest.approx(0.20017, rel=1e-4)

        k1 = "Dilute Side Outlet"
        k2 = "Molar Flowrate ('Liq', 'Cl_-')"
        assert value(table[k1][k2]) == pytest.approx(0.20017, rel=1e-4)

        k1 = "Concentrate Side Outlet"
        k2 = "Molar Flowrate ('Liq', 'H2O')"
        assert value(table[k1][k2]) == pytest.approx(1.02747, rel=1e-4)

        k1 = "Concentrate Side Outlet"
        k2 = "Molar Flowrate ('Liq', 'NaCl')"
        assert value(table[k1][k2]) == pytest.approx(0.0015, rel=1e-4)

        k1 = "Concentrate Side Outlet"
        k2 = "Molar Flowrate ('Liq', 'Na_+')"
        assert value(table[k1][k2]) == pytest.approx(0.149834, rel=1e-4)

        k1 = "Concentrate Side Outlet"
        k2 = "Molar Flowrate ('Liq', 'Cl_-')"
        assert value(table[k1][k2]) == pytest.approx(0.149834, rel=1e-4)
