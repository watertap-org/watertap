import pytest
from pyomo.environ import ConcreteModel, SolverFactory, TerminationCondition, \
    value, Constraint, Var, Objective, Expression, SolverStatus
from pyomo.environ import units as pyunits
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent, check_units_equivalent
import pyomo.util.infeasible as infeas
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom, number_variables
import idaes.core.util.model_statistics as stats
from idaes.core.util.constants import Constants
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.property_models.ion_DSPMDE_prop_pack import DSPMDEParameterBlock
from watertap.unit_models.electrodialysis_0D import Electrodialysis0D

from idaes.core.util import get_solver

# Get default solver for testing
solver = get_solver()

class TestElec0D():
    @pytest.fixture(scope="class")
    def ED_model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})
        # create dict to define ions (the prop pack of Adam requires this)
        ion_dict = {
            "solute_list": ["Na_+", "Cl_-"],
            "diffusivity_data": {("Liq", "Na_+"): 1.33e-9,
                                 ("Liq", "Cl_-"): 2.03e-9},
            "mw_data": {"H2O": 18e-3,
                        "Na_+": 23e-3,
                        "Cl_-": 35e-3},
            "stokes_radius_data": {"Na_+": 0.184e-9,
                                   "Cl_-": 0.121e-9},
            "charge": {"Na_+": 1,
                       "Cl_-": -1},
        }
        # attach prop pack to flowsheet
        m.fs.properties = DSPMDEParameterBlock(default=ion_dict)

        # build the unit model, pass prop pack to the model
        m.fs.unit = Electrodialysis0D(default={
            "property_package": m.fs.properties })

        return m

    @pytest.mark.unit
    def test_stats(self, ED_model):
        m = ED_model
        assert degrees_of_freedom(m.fs) == 8

        assert_units_consistent(m)

        m.fs.unit.inlet_dilute.pressure.fix(101325)
        m.fs.unit.inlet_dilute.temperature.fix(298.15)
        m.fs.unit.inlet_dilute.flow_mol_phase_comp[0, 'Liq', 'H2O'].fix(1)
        m.fs.unit.inlet_dilute.flow_mol_phase_comp[0, 'Liq', 'Na_+'].fix(0.01)
        m.fs.unit.inlet_dilute.flow_mol_phase_comp[0, 'Liq', 'Cl_-'].fix(0.01)

        m.fs.unit.inlet_concentrate.pressure.fix(101325)
        m.fs.unit.inlet_concentrate.temperature.fix(298.15)
        m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, 'Liq', 'H2O'].fix(1)
        m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, 'Liq', 'Na_+'].fix(0.01)
        m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, 'Liq', 'Cl_-'].fix(0.01)

        assert degrees_of_freedom(m.fs) == 0

    @pytest.mark.component
    def test_initialize(self, ED_model):
        m = ED_model

        m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1, index=('Liq', 'H2O'))
        m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1e2, index=('Liq', 'Na_+'))
        m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1e2, index=('Liq', 'Cl_-'))

        iscale.calculate_scaling_factors(m.fs)

        # Expected to fail here
        m.fs.unit.initialize(optarg=solver.options, outlvl=idaeslog.DEBUG)

        assert degrees_of_freedom(m.fs) == 0

    @pytest.mark.component
    def test_solve(self, ED_model):
        m = ED_model

        # Expected to fail here
        results = solver.solve(m, tee=True)

        assert results.solver.termination_condition == \
               TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok
