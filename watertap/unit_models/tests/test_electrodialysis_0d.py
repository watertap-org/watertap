###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
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
from watertap.property_models.ion_DSPMDE_prop_pack import DSPMDEParameterBlock 
from electrodialysis_0d import Electrodialysis0D
from pyomo.environ import (ConcreteModel,
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
                           TerminationCondition)
from idaes.core import (FlowsheetBlock,
                        MaterialFlowBasis,
                        MaterialBalanceType,
                        EnergyBalanceType)
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
import idaes.core.util.scaling as iscale
from idaes.core.util.scaling import badly_scaled_var_generator
from idaes.core.util.testing import initialization_tester
from idaes.core.util import get_solver
import re

__author__ = "Xiangyu Bi, Austin Ladshaw"

solver = get_solver()

# -----------------------------------------------------------------------------
# Start test class
class TestElectrodialysis():
    @pytest.fixture(scope="class")
    def electrodialysis_cell(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.properties =  DSPMDEParameterBlock()
        m.fs.unit = Electrodialysis0D(default = {"property_package": m.fs.properties})
        
        return m

    @pytest.mark.unit
    def test_build_model(self, electrodialysis_cell):
        m = electrodialysis_cell
        assert m.fs.unit.config.operation_mode == 'Constant Current' or m.fs.unit.config.operation_mode == 'Constant Voltage'
        assert isinstance(m.fs.unit.ion_set, Set)
        assert isinstance(m.fs.unit.membrane_set, Set)
        assert isinstance(m.fs.unit.water_density, Param)
        assert isinstance(m.fs.unit.water_MW, Param)
        assert isinstance(m.fs.unit.cell_width, Var)
        assert isinstance(m.fs.unit.cell_length, Var)
        assert isinstance(m.fs.unit.spacer_thickness, Var)
        assert isinstance(m.fs.unit.membrane_thickness, Var)
        assert isinstance(m.fs.unit.ion_diffusivity_membrane, Var)
        assert isinstance(m.fs.unit.ion_trans_number_membrane, Var)
        assert isinstance(m.fs.unit.water_trans_number_membrane, Var)
        assert isinstance(m.fs.unit.water_permeability_membrane, Var)
        assert isinstance(m.fs.unit.membrane_surface_resistence, Var)
        assert isinstance(m.fs.unit.current, Var)
        assert isinstance(m.fs.unit.voltage, Var)
        assert isinstance(m.fs.unit.current_utilization, Var)
        assert isinstance(m.fs.unit.elec_migration_flux_in, Var)
        assert isinstance(m.fs.unit.elec_migration_flux_out, Var)
        assert isinstance(m.fs.unit.nonelec_flux_in, Var)
        assert isinstance(m.fs.unit.nonelec_flux_out, Var)   
        assert isinstance(m.fs.unit.eq_current_voltage_relation, Constraint)
        assert isinstance(m.fs.unit.eq_elec_migration_flux_in, Constraint)
        assert isinstance(m.fs.unit.eq_elec_migration_flux_out, Constraint)
        assert isinstance(m.fs.unit.eq_nonelec_flux_in, Constraint)
        assert isinstance(m.fs.unit.eq_nonelec_flux_out, Constraint)
        assert isinstance(m.fs.unit.eq_mass_transfer_term_diluate, Constraint)
        assert isinstance(m.fs.unit.eq_mass_transfer_term_concentrate, Constraint)
        assert isinstance(m.fs.unit.eq_isothermal_diluate, Constraint)
        assert isinstance(m.fs.unit.eq_isothermal_concentrate, Constraint)
    
    @pytest.mark.unit
    def test_stats(self, electrodialysis_cell):
        m = electrodialysis_cell

        # Check to make sure we have the correct DOF and
        #   check to make sure the units are correct
        assert_units_consistent(m)
        assert degrees_of_freedom(m) == 31
    #Specify a system
        ion_dict = {
            "solute_list": ["Na_+", "Cl_-"],
            "mw_data": {"H2O": 18e-3,
                        "Na_+": 23e-3,
                        "Cl_-": 35.5e-3},
            "electrical_mobility_data": {"Na_+":5.19e-8, 
                                         "Cl_-":7.92e-8},            
            "charge": {"Na_+": 1,
                       "Cl_-": -1},
        } 
        m.fs.properties = DSPMDEParameterBlock(default = ion_dict)
        # set the operational parameters
        m.fs.unit.water_trans_number_membrane.fix(1.9)
        m.fs.unit.water_permeability_membrane['cem'].fix(2.16e-14)
        m.fs.unit.water_permeability_membrane['aem'].fix(1.75e-14)
        m.fs.unit.voltage.fix(0.5)
        #m.fs.unit.current.fix(15)
        m.fs.unit.current_utilization.fix(1)
        m.fs.unit.spacer_thickness.fix(1.5e-4)
        m.fs.unit.membrane_surface_resistence['cem'].fix(1.89e-4)
        m.fs.unit.membrane_surface_resistence['aem'].fix(1.77e-4)
        m.fs.unit.cell_width.fix(0.1)
        m.fs.unit.cell_length.fix(0.43)
        m.fs.unit.membrane_thickness['aem'].fix(1.3e-4)
        m.fs.unit.membrane_thickness['cem'].fix(1.3e-4)
        m.fs.unit.ion_diffusivity_membrane.fix(7e-9)
        m.fs.unit.ion_trans_number_membrane['cem','Na_+'].fix(1)
        m.fs.unit.ion_trans_number_membrane['aem','Na_+'].fix(0)
        m.fs.unit.ion_trans_number_membrane['cem','Cl_-'].fix(0)
        m.fs.unit.ion_trans_number_membrane['aem','Cl_-'].fix(1)

        assert sum(m.fs.unit.ion_trans_number_membrane['cem', j] for j in m.fs.properties.solute_set == 1)
        assert sum(m.fs.unit.ion_trans_number_membrane['aem', j] for j in m.fs.properties.solute_set == 1)
        


        # set the inlet stream
        m.fs.unit.inlet_diluate.pressure.fix(101325)
        m.fs.unit.inlet_diluate.temperature.fix(298.15)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, 'Liq', 'H2O'].fix(0.013)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, 'Liq', 'Na_+'].fix(2.46e-5)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, 'Liq', 'Cl_-'].fix(2.46e-5)
        m.fs.unit.inlet_concentrate.pressure.fix(101325)
        m.fs.unit.inlet_concentrate.temperature.fix(298.15)
        m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, 'Liq', 'H2O'].fix(0.013)
        m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, 'Liq', 'Na_+'].fix(2.46e-5)
        m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, 'Liq', 'Cl_-'].fix(2.46e-5)


        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_scaling(self, electrodialysis_cell):
         m = electrodialysis_cell

        # Set some scaling factors and look for 'bad' scaling
        model.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
        model.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'TSS'))
        model.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'TDS'))
        model.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'Sludge'))
        iscale.calculate_scaling_factors(model.fs)

        # check that all variables have scaling factors
        unscaled_var_list = list(iscale.unscaled_variables_generator(model))
        assert len(unscaled_var_list) == 0

        # check if any variables are badly scaled
        badly_scaled_var_values = {
        var.name: val
        for (var, val) in iscale.badly_scaled_var_generator(
            model, large=1e2, small=1e-2
            )
        }
        assert not badly_scaled_var_values

    @pytest.mark.component
    def test_initialization(self, coag_obj_w_chems):
        model = coag_obj_w_chems
        initialization_tester(model)

        # check to make sure DOF does not change
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_solve(self, coag_obj_w_chems):
        model = coag_obj_w_chems

        # first, check to make sure that after initialized, the scaling is still good
        badly_scaled_var_values = {
        var.name: val
        for (var, val) in iscale.badly_scaled_var_generator(
            model, large=1e2, small=1e-2
            )
        }
        assert not badly_scaled_var_values

        # run solver and check for optimal solution
        results = solver.solve(model)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, coag_obj_w_chems):
        model = coag_obj_w_chems

        assert value(model.fs.unit.outlet.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(1,  rel=1e-4)
        assert value(model.fs.unit.outlet.flow_mass_phase_comp[0, 'Liq', 'Sludge']) == pytest.approx(0.00999,  rel=1e-4)
        assert value(model.fs.unit.outlet.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(0.010015,  rel=1e-4)
        assert value(model.fs.unit.outlet.flow_mass_phase_comp[0, 'Liq', 'TSS']) == pytest.approx(9.36352e-06,  rel=1e-4)


# -----------------------------------------------------------------------------
# Start test class without chemicals added
class TestCoagulation_withNoChemicals():
    @pytest.fixture(scope="class")
    def coag_obj_wo_chems(self):
        model = ConcreteModel()
        model.fs = FlowsheetBlock(default={"dynamic": False})
        model.fs.properties = CoagulationParameterBlock()
        model.fs.unit = CoagulationFlocculation(default={
            "property_package": model.fs.properties})

        return model

    @pytest.mark.unit
    def test_build_model(self, coag_obj_wo_chems):
        model = coag_obj_wo_chems

        assert len(model.fs.unit.config.chemical_additives) == 0
        assert isinstance(model.fs.unit.slope, Var)
        assert isinstance(model.fs.unit.intercept, Var)
        assert isinstance(model.fs.unit.initial_turbidity_ntu, Var)
        assert isinstance(model.fs.unit.final_turbidity_ntu, Var)
        assert isinstance(model.fs.unit.chemical_doses, Var)
        assert len(model.fs.unit.chemical_doses) == 0
        assert isinstance(model.fs.unit.chemical_mw, Param)
        assert len(model.fs.unit.chemical_mw) == 0
        assert isinstance(model.fs.unit.salt_mw, Param)
        assert len(model.fs.unit.salt_mw) == 0
        assert isinstance(model.fs.unit.salt_from_additive_mole_ratio, Param)
        assert len(model.fs.unit.salt_from_additive_mole_ratio) == 0
        assert isinstance(model.fs.unit.tss_loss_rate, Var)
        assert isinstance(model.fs.unit.eq_tss_loss_rate, Constraint)

        assert not hasattr(model.fs.unit, "tds_gain_rate")
        assert not hasattr(model.fs.unit, "eq_tds_gain_rate")

        assert isinstance(model.fs.unit.eq_mass_transfer_term, Constraint)

    @pytest.mark.unit
    def test_stats(self, coag_obj_wo_chems):
        model = coag_obj_wo_chems

        # Check to make sure we have the correct DOF and
        #   check to make sure the units are correct
        assert_units_consistent(model)
        assert degrees_of_freedom(model) == 9

        # set the operational parameters
        model.fs.unit.fix_tss_turbidity_relation_defaults()
        model.fs.unit.initial_turbidity_ntu.fix(5000)
        model.fs.unit.final_turbidity_ntu.fix(100)

        # set the inlet streams
        assert degrees_of_freedom(model) == 6

        tss_in = value(model.fs.unit.compute_inlet_tss_mass_flow(0))
        assert tss_in == pytest.approx(0.0093,  rel=1e-4)
        model.fs.unit.inlet.pressure.fix(101325)
        model.fs.unit.inlet.temperature.fix(298.15)
        model.fs.unit.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(1)
        model.fs.unit.inlet.flow_mass_phase_comp[0, 'Liq', 'TDS'].fix(0.01)
        model.fs.unit.inlet.flow_mass_phase_comp[0, 'Liq', 'TSS'].fix(tss_in)
        model.fs.unit.inlet.flow_mass_phase_comp[0, 'Liq', 'Sludge'].fix(0.0)

        assert degrees_of_freedom(model) == 0

    @pytest.mark.unit
    def test_scaling(self, coag_obj_wo_chems):
        model = coag_obj_wo_chems

        # Set some scaling factors and look for 'bad' scaling
        model.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
        model.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'TSS'))
        model.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'TDS'))
        model.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'Sludge'))
        #iscale.set_scaling_factor(model.fs.unit.tss_loss_rate, 100)
        iscale.calculate_scaling_factors(model.fs)

        # check that all variables have scaling factors
        unscaled_var_list = list(iscale.unscaled_variables_generator(model))
        assert len(unscaled_var_list) == 0

        # check if any variables are badly scaled
        badly_scaled_var_values = {
        var.name: val
        for (var, val) in iscale.badly_scaled_var_generator(
            model, large=1e2, small=1e-2
            )
        }
        print(iscale.get_scaling_factor(model.fs.unit.tss_loss_rate))
        assert not badly_scaled_var_values

    @pytest.mark.component
    def test_initialization(self, coag_obj_wo_chems):
        model = coag_obj_wo_chems
        initialization_tester(model)

        # check to make sure DOF does not change
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_solve(self, coag_obj_wo_chems):
        model = coag_obj_wo_chems

        # first, check to make sure that after initialized, the scaling is still good
        badly_scaled_var_values = {
        var.name: val
        for (var, val) in iscale.badly_scaled_var_generator(
            model, large=1e2, small=1e-2
            )
        }
        assert not badly_scaled_var_values

        # run solver and check for optimal solution
        results = solver.solve(model)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, coag_obj_wo_chems):
        model = coag_obj_wo_chems

        assert value(model.fs.unit.outlet.flow_mass_phase_comp[0, 'Liq', 'H2O']) == pytest.approx(1,  rel=1e-4)
        assert value(model.fs.unit.outlet.flow_mass_phase_comp[0, 'Liq', 'Sludge']) == pytest.approx(0.0091127,  rel=1e-4)
        assert value(model.fs.unit.outlet.flow_mass_phase_comp[0, 'Liq', 'TDS']) == pytest.approx(0.01,  rel=1e-4)
        assert value(model.fs.unit.outlet.flow_mass_phase_comp[0, 'Liq', 'TSS']) == pytest.approx(0.00018725,  rel=1e-4)

# -----------------------------------------------------------------------------
# Start test class with bad config
class TestCoagulation_withBadConfig():
    @pytest.fixture(scope="class")
    def coag_obj_bad_config(self):
        model = ConcreteModel()
        model.fs = FlowsheetBlock(default={"dynamic": False})
        model.fs.properties = CoagulationParameterBlock()

        return model

    @pytest.mark.unit
    def test_build_model_catch_errors(self, coag_obj_bad_config):
        model = coag_obj_bad_config

        bad_dict1 = {'Alum':
                        {"foo_bar":
                            {"mw_additive": (200, pyunits.g/pyunits.mol),
                            "moles_salt_per_mole_additive": 3,
                            "mw_salt": (100, pyunits.g/pyunits.mol)}
                        }
                     }
        with pytest.raises(ConfigurationError, match="Did not provide a 'parameter_data' for chemical"):
            model.fs.unit = CoagulationFlocculation(default={
                "property_package": model.fs.properties,
                "chemical_additives": bad_dict1 })

        bad_dict2 = {'Alum':
                        {"parameter_data":
                            {"foo_bar": (200, pyunits.g/pyunits.mol),
                            "moles_salt_per_mole_additive": 3,
                            "mw_salt": (100, pyunits.g/pyunits.mol)}
                        }
                     }
        with pytest.raises(ConfigurationError, match="Did not provide a 'mw_additive' for chemical"):
            model.fs.unit = CoagulationFlocculation(default={
                "property_package": model.fs.properties,
                "chemical_additives": bad_dict2 })

        bad_dict3 = {'Alum':
                        {"parameter_data":
                            {"mw_additive": (200, pyunits.g/pyunits.mol),
                            "moles_salt_per_mole_additive": "foo-bar",
                            "mw_salt": (100, pyunits.g/pyunits.mol)}
                        }
                     }
        with pytest.raises(ConfigurationError, match="Did not provide a number for 'moles_salt_per_mole_additive'"):
            model.fs.unit = CoagulationFlocculation(default={
                "property_package": model.fs.properties,
                "chemical_additives": bad_dict3 })

        bad_dict4 = {'Alum':
                        {"parameter_data":
                            {"mw_additive": (200, pyunits.g/pyunits.mol),
                            "foo-bar": 3,
                            "mw_salt": (100, pyunits.g/pyunits.mol)}
                        }
                     }
        with pytest.raises(ConfigurationError, match="Did not provide a 'moles_salt_per_mole_additive' for chemical"):
            model.fs.unit = CoagulationFlocculation(default={
                "property_package": model.fs.properties,
                "chemical_additives": bad_dict4 })

        bad_dict5 = {'Alum':
                        {"parameter_data":
                            {"mw_additive": (200, pyunits.g/pyunits.mol),
                            "moles_salt_per_mole_additive": 3,
                            "foo-bar": (100, pyunits.g/pyunits.mol)}
                        }
                     }
        with pytest.raises(ConfigurationError, match="Did not provide a 'mw_salt' for chemical"):
            model.fs.unit = CoagulationFlocculation(default={
                "property_package": model.fs.properties,
                "chemical_additives": bad_dict5 })

        bad_dict6 = {'Alum':
                        {"parameter_data":
                            {"mw_additive": "not a tuple",
                            "moles_salt_per_mole_additive": 3,
                            "mw_salt": (100, pyunits.g/pyunits.mol)}
                        }
                     }
        with pytest.raises(ConfigurationError, match="Did not provide a tuple for 'mw_additive'"):
            model.fs.unit = CoagulationFlocculation(default={
                "property_package": model.fs.properties,
                "chemical_additives": bad_dict6 })

# -----------------------------------------------------------------------------
# Start test class with bad config
class TestCoagulation_withBadProperties():
    @pytest.fixture(scope="class")
    def coag_obj_bad_properties(self):
        model = ConcreteModel()
        model.fs = FlowsheetBlock(default={"dynamic": False})

        return model

    @pytest.mark.unit
    def test_build_model_catch_prop_errors(self, coag_obj_bad_properties):
        model = coag_obj_bad_properties

        error_msg = ("Coagulation-Flocculation model MUST contain ('Liq','TDS') "
                    "as a component, but the property package has only specified "
                    "the following components [('Liq', 'H2O'), ('Liq', 'NaCl')]")
        with pytest.raises(ConfigurationError, match=re.escape(error_msg)):
            model.fs.properties = NaClParameterBlock()
            model.fs.unit = CoagulationFlocculation(default={
                "property_package": model.fs.properties })

        error_msg = ("Coagulation-Flocculation model MUST contain ('Liq','Sludge') "
                    "as a component, but the property package has only specified "
                    "the following components [('Liq', 'H2O'), ('Liq', 'TDS')]")
        with pytest.raises(ConfigurationError, match=re.escape(error_msg)):
            model.fs.properties = SeawaterParameterBlock()
            model.fs.unit = CoagulationFlocculation(default={
                "property_package": model.fs.properties })

        # NOTE: package must also contain ('Liq','TSS') as a component
