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

"""
   Tests to check reactions with water dissociation and calcium hydroxide (CaOH2)
   with equilibrium and stoichiometric reactor and return correct pH value.
"""
import pytest

# Importing the object for units from pyomo
from pyomo.environ import units as pyunits

from pyomo.environ import (ConcreteModel,
                           TerminationCondition,
                           SolverStatus,
                           value,
                           Var,
                           Suffix)

from pyomo.util.check_units import (assert_units_consistent,
                                    assert_units_equivalent)

from idaes.core import (FlowsheetBlock,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType)

import idaes.logger as idaeslog

from idaes.generic_models.unit_models.stoichiometric_reactor import \
    StoichiometricReactor

from idaes.core.util import scaling as iscale

from idaes.core.util import get_solver
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              fixed_variables_set,
                                              activated_constraints_set,
                                              number_variables,
                                              number_total_constraints,
                                              number_unused_variables)
from idaes.core.util.testing import (PhysicalParameterTestBlock,
                                     ReactionParameterTestBlock,
                                     initialization_tester)

from idaes.core import AqueousPhase, FlowsheetBlock, EnergyBalanceType
from idaes.core.components import Solvent, Solute, Cation, Anion
from idaes.core.phases import PhaseType as PT

# Import the idaes objects for Generic Properties and Reactions
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)
from idaes.generic_models.properties.core.generic.generic_reaction import (
        GenericReactionParameterBlock)

# Imports from idaes generic models
import idaes.generic_models.properties.core.pure.Perrys as Perrys
from idaes.generic_models.properties.core.pure.ConstantProperties import Constant
from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.ideal import Ideal
from idaes.generic_models.properties.core.reactions.rate_constant import arrhenius
from idaes.generic_models.properties.core.reactions.rate_forms import power_law_rate

from idaes.generic_models.unit_models.cstr import CSTR

# Importing the enum for concentration unit basis used in the 'get_concentration_term' function
from idaes.generic_models.properties.core.generic.generic_reaction import ConcentrationForm

# Import the object/function for heat of reaction
from idaes.generic_models.properties.core.reactions.dh_rxn import constant_dh_rxn

# Import safe log power law equation
from idaes.generic_models.properties.core.reactions.equilibrium_forms import log_power_law_equil

# Import k-value functions
from idaes.generic_models.properties.core.reactions.equilibrium_constant import (
    gibbs_energy,
    van_t_hoff)

# Import log10 function from pyomo
from pyomo.environ import log10

thermo_config = {
    "components": {
        'H2O': {"type": Solvent,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (18.0153, pyunits.g/pyunits.mol),
                    "pressure_crit": (220.64E5, pyunits.Pa),
                    "temperature_crit": (647, pyunits.K),
                    # Comes from Perry's Handbook:  p. 2-98
                    "dens_mol_liq_comp_coeff": {
                        '1': (5.459, pyunits.kmol*pyunits.m**-3),
                        '2': (0.30542, pyunits.dimensionless),
                        '3': (647.13, pyunits.K),
                        '4': (0.081, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-285.830, pyunits.kJ/pyunits.mol),
                    "enth_mol_form_vap_comp_ref": (0, pyunits.kJ/pyunits.mol),
                    # Comes from Perry's Handbook:  p. 2-174
                    "cp_mol_liq_comp_coeff": {
                        '1': (2.7637E5, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (-2.0901E3, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (8.125, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (-1.4116E-2, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (9.3701E-6, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "cp_mol_ig_comp_coeff": {
                        'A': (30.09200, pyunits.J/pyunits.mol/pyunits.K),
                        'B': (6.832514, pyunits.J*pyunits.mol**-1*pyunits.K**-1*pyunits.kiloK**-1),
                        'C': (6.793435, pyunits.J*pyunits.mol**-1*pyunits.K**-1*pyunits.kiloK**-2),
                        'D': (-2.534480, pyunits.J*pyunits.mol**-1*pyunits.K**-1*pyunits.kiloK**-3),
                        'E': (0.082139, pyunits.J*pyunits.mol**-1*pyunits.K**-1*pyunits.kiloK**2),
                        'F': (-250.8810, pyunits.kJ/pyunits.mol),
                        'G': (223.3967, pyunits.J/pyunits.mol/pyunits.K),
                        'H': (0, pyunits.kJ/pyunits.mol)},
                    "entr_mol_form_liq_comp_ref": (69.95, pyunits.J/pyunits.K/pyunits.mol),
                    "pressure_sat_comp_coeff": {
                        'A': (4.6543, None),  # [1], temperature range 255.9 K - 373 K
                        'B': (1435.264, pyunits.K),
                        'C': (-64.848, pyunits.K)}
                                },
                    # End parameter_data
                    },
        'H_+': {"type": Cation, "charge": 1,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (1.00784, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (5.459, pyunits.kmol*pyunits.m**-3),
                        '2': (0.30542, pyunits.dimensionless),
                        '3': (647.13, pyunits.K),
                        '4': (0.081, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-230.000, pyunits.kJ/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (2.7637E5, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (-2.0901E3, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (8.125, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (-1.4116E-2, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (9.3701E-6, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "entr_mol_form_liq_comp_ref": (-10.75, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'OH_-': {"type": Anion,
                "charge": -1,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (17.008, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (5.459, pyunits.kmol*pyunits.m**-3),
                        '2': (0.30542, pyunits.dimensionless),
                        '3': (647.13, pyunits.K),
                        '4': (0.081, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-230.000, pyunits.kJ/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (2.7637E5, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (-2.0901E3, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (8.125, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (-1.4116E-2, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (9.3701E-6, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "entr_mol_form_liq_comp_ref": (-10.75, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    },
        'Ca(OH)2': {  "type": Solute,  "valid_phase_types": PT.aqueousPhase,
                    # Define the methods used to calculate the following properties
                    "dens_mol_liq_comp": Constant,
                    "enth_mol_liq_comp": Constant,
                    "cp_mol_liq_comp": Constant,
                    "entr_mol_liq_comp": Constant,
                    # Parameter data is always associated with the methods defined above
                    "parameter_data": {
                        "mw": (74.093, pyunits.g/pyunits.mol),
                        "dens_mol_liq_comp_coeff": (55, pyunits.kmol*pyunits.m**-3),
                        "enth_mol_form_liq_comp_ref": (-945.53, pyunits.kJ/pyunits.mol),
                        "cp_mol_liq_comp_coeff": (167039, pyunits.J/pyunits.kmol/pyunits.K),
                        "entr_mol_form_liq_comp_ref": (100, pyunits.J/pyunits.K/pyunits.mol)
                            },
                    # End parameter_data
                    },
        'Ca_2+': {"type": Cation, "charge": 2,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (40.078, pyunits.g/pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        '1': (13.5, pyunits.kmol*pyunits.m**-3),
                        '2': (1, pyunits.dimensionless),
                        '3': (1, pyunits.K),
                        '4': (1, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-542.83, pyunits.J/pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        '1': (2.7637E5, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (-2.0901E3, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (8.125, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (-1.4116E-2, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (9.3701E-6, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "entr_mol_form_liq_comp_ref": (-53, pyunits.J/pyunits.K/pyunits.mol)
                                },
                    # End parameter_data
                    }
              },
              # End Component list
        "phases":  {'Liq': {"type": AqueousPhase,
                            "equation_of_state": Ideal},
                    },

        "state_definition": FTPx,
        "state_bounds": {"flow_mol": (0, 50, 100),
                         "temperature": (273.15, 300, 650),
                         "pressure": (5e4, 1e5, 1e6)
                     },
        "pressure_ref": 1e5,
        "temperature_ref": 300,
        "base_units": {"time": pyunits.s,
                       "length": pyunits.m,
                       "mass": pyunits.kg,
                       "amount": pyunits.mol,
                       "temperature": pyunits.K},
        # Inherent reactions
        "inherent_reactions": {
            "H2O_Kw": {
                    "stoichiometry": {("Liq", "H2O"): -1,
                                     ("Liq", "H_+"): 1,
                                     ("Liq", "OH_-"): 1},
                   "heat_of_reaction": constant_dh_rxn,
                   "equilibrium_constant": van_t_hoff,
                   "equilibrium_form": log_power_law_equil,
                   "concentration_form": ConcentrationForm.moleFraction,
                   "parameter_data": {
                       "dh_rxn_ref": (55.830, pyunits.J/pyunits.mol),
                       "k_eq_ref": (10**-14/55.2/55.2, pyunits.dimensionless),
                       "T_eq_ref": (298, pyunits.K),

                       # By default, reaction orders follow stoichiometry
                       #    manually set reaction order here to override
                       "reaction_order": {("Liq", "H2O"): 0,
                                        ("Liq", "H_+"): 1,
                                        ("Liq", "OH_-"): 1}
                        }
                        # End parameter_data
                   }
                   # End R1
             }
             # End equilibrium_reactions

}
    # End thermo_config definition


# This config is REQUIRED to use EquilibriumReactor even if we have no equilibrium reactions
reaction_config = {
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},
    "equilibrium_reactions": {
        "dummy": {
                "stoichiometry": {},
                "equilibrium_form": log_power_law_equil,
               }
         },
         # End equilibrium_reactions
    "rate_reactions": {
        "R1": {"stoichiometry": {("Liq", "Ca(OH)2"): -1,
                                 ("Liq", "Ca_2+"): 1,
                                 ("Liq", "OH_-"): 2},
               "heat_of_reaction": constant_dh_rxn,
               "rate_constant" : arrhenius,
               "rate_form" : power_law_rate,
               "concentration_form" : ConcentrationForm.moleFraction,
               "parameter_data": {
                   "arrhenius_const" : (1, pyunits.mol/pyunits.m**3/pyunits.s),
                   "energy_activation" : (0, pyunits.J/pyunits.mol),
                   "dh_rxn_ref": (0, pyunits.J/pyunits.mol)
              }
         }
    }
}
    # End reaction_config definition


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------
class TestWaterStoich(object):
    @pytest.fixture(scope="class")
    def water_stoich(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.thermo_params = GenericParameterBlock(default=thermo_config)
        m.fs.rxn_params = GenericReactionParameterBlock(
                default={"property_package": m.fs.thermo_params, **reaction_config})

        m.fs.unit = StoichiometricReactor(default={
                "property_package": m.fs.thermo_params,
                "reaction_package": m.fs.rxn_params,
                "has_heat_transfer": False,
                "has_heat_of_reaction": False,
                "energy_balance_type": EnergyBalanceType.none,
                "has_pressure_change": False})

        m.fs.unit.inlet.mole_frac_comp[0, "Ca_2+"].fix( 0. )
        m.fs.unit.inlet.mole_frac_comp[0, "Ca(OH)2"].fix( 0.0000018 )
        m.fs.unit.inlet.mole_frac_comp[0, "H_+"].fix( 0. )
        m.fs.unit.inlet.mole_frac_comp[0, "OH_-"].fix( 0. )
        m.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix( 0.9999982 )

        m.fs.unit.inlet.pressure.fix(101325.0)
        m.fs.unit.inlet.temperature.fix(298.)
        m.fs.unit.inlet.flow_mol.fix(10)

        m.fs.unit.outlet.temperature.fix(298.)
        m.fs.unit.outlet.mole_frac_comp[0, "Ca(OH)2"].fix( 0.0000018*0.95 )
        #m.fs.unit.rate_reaction_extent[0, 'R1'].fix(0)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, water_stoich):

        m = water_stoich 
        assert hasattr(m.fs.thermo_params, 'component_list')
        assert len(m.fs.thermo_params.component_list) == 5
        assert 'H2O' in m.fs.thermo_params.component_list
        assert isinstance(m.fs.thermo_params.H2O, Solvent)
        assert 'H_+' in m.fs.thermo_params.component_list
        assert isinstance(m.fs.thermo_params.component('H_+'), Cation)
        assert 'OH_-' in m.fs.thermo_params.component_list
        assert isinstance(m.fs.thermo_params.component('OH_-'), Anion)

        assert hasattr(m.fs.thermo_params, 'phase_list')
        assert len(m.fs.thermo_params.phase_list) == 1
        assert isinstance(m.fs.thermo_params.Liq, AqueousPhase)

    @pytest.mark.unit
    def test_units_stoich(self, water_stoich):
        m = water_stoich
        assert_units_consistent(m)

    @pytest.mark.unit
    def test_dof_stoich(self, water_stoich):
        m = water_stoich
        assert (degrees_of_freedom(m) == 0)

    @pytest.mark.unit
    def test_stats_stoich(self, water_stoich):
        m = water_stoich
        assert (number_variables(m) == 121)
        assert (number_total_constraints(m) == 40)
        assert (number_unused_variables(m) == 67)

    @pytest.mark.component
    def test_scaling_stoich(self, water_stoich):
        m = water_stoich

        # Iterate through the reactions to set appropriate eps values
        factor = 1e-4
        for rid in m.fs.thermo_params.inherent_reaction_idx:
            scale = value(m.fs.unit.control_volume.properties_out[0.0].k_eq[rid].expr)
            # Want to set eps in some fashion similar to this
            if scale < 1e-16:
                m.fs.thermo_params.component("reaction_"+rid).eps.value = scale*factor
            else:
                m.fs.thermo_params.component("reaction_"+rid).eps.value = 1e-16*factor

        for i in m.fs.unit.control_volume.inherent_reaction_extent_index:
            scale = value(m.fs.unit.control_volume.properties_out[0.0].k_eq[i[1]].expr)
            iscale.set_scaling_factor(m.fs.unit.control_volume.inherent_reaction_extent[0.0,i[1]], 10/scale)
            iscale.constraint_scaling_transform(m.fs.unit.control_volume.properties_out[0.0].
                    inherent_equilibrium_constraint[i[1]], 0.1)

        # Scaling for kinetic reactions
        for i in m.fs.rxn_params.rate_reaction_idx:
            scale = value(m.fs.unit.control_volume.reactions[0.0].reaction_rate[i].expr)
            iscale.set_scaling_factor(m.fs.unit.control_volume.rate_reaction_extent[0.0,i], 10/scale)

        # Next, try adding scaling for species
        min = 1e-6
        for i in m.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp:
            # i[0] = phase, i[1] = species
            if m.fs.unit.inlet.mole_frac_comp[0, i[1]].value > min:
                scale = m.fs.unit.inlet.mole_frac_comp[0, i[1]].value
            else:
                scale = min
            iscale.set_scaling_factor(m.fs.unit.control_volume.properties_out[0.0].mole_frac_comp[i[1]], 10/scale)
            iscale.set_scaling_factor(m.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[i], 10/scale)
            iscale.set_scaling_factor(m.fs.unit.control_volume.properties_out[0.0].flow_mol_phase_comp[i], 10/scale)
            iscale.constraint_scaling_transform(
                m.fs.unit.control_volume.properties_out[0.0].component_flow_balances[i[1]], 10/scale)
            iscale.constraint_scaling_transform(m.fs.unit.control_volume.material_balances[0.0,i[1]], 10/scale)

        iscale.set_scaling_factor(m.fs.unit.control_volume.rate_reaction_extent[0.0,'R1'], 1)
        iscale.calculate_scaling_factors(m.fs.unit)

        m.fs.unit.control_volume.pprint()
        assert False

        assert isinstance(m.fs.unit.control_volume.scaling_factor, Suffix)

        assert isinstance(m.fs.unit.control_volume.properties_out[0.0].scaling_factor, Suffix)

        assert isinstance(m.fs.unit.control_volume.properties_in[0.0].scaling_factor, Suffix)

    @pytest.mark.component
    def test_initialize_inherent(self, water_stoich):
        m = water_stoich
        state_args = {'mole_frac_comp':
                        {   'Ca(OH)2': 0.0000018,
                            'Ca_2+': 0.0000018,
                            'H2O': 1,
                            'H_+': 10**-7/55.6,
                            'OH_-': 10**-7/55.6
                        },
                        'pressure': 101325,
                        'temperature': 298,
                        'flow_mol': 10
                    }
        orig_fixed_vars = fixed_variables_set(m)
        orig_act_consts = activated_constraints_set(m)

        solver.options['bound_push'] = 1e-10
        solver.options['mu_init'] = 1e-6
        m.fs.unit.initialize(state_args=state_args, optarg=solver.options, outlvl=idaeslog.DEBUG)

        fin_fixed_vars = fixed_variables_set(m)
        fin_act_consts = activated_constraints_set(m)

        print(value(m.fs.unit.outlet.temperature[0]))
        assert degrees_of_freedom(m) == 0

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

    @pytest.mark.component
    def test_solve_inherent(self, water_stoich):
        m = water_stoich
        solver.options['max_iter'] = 4000
        results = solver.solve(m)
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_solution_inherent(self, water_stoich):
        m = water_stoich

        assert pytest.approx(10, rel=1e-5) == value(m.fs.unit.outlet.flow_mol[0])
        print(value(m.fs.unit.outlet.flow_mol[0]))
        assert pytest.approx(101325, rel=1e-5) == value(m.fs.unit.outlet.pressure[0])
        print(value(m.fs.unit.outlet.pressure[0]))

        total_molar_density = \
            value(m.fs.unit.control_volume.properties_out[0.0].dens_mol_phase['Liq'])/1000
        print(total_molar_density)
        assert pytest.approx(55.2336, rel=1e-5) == total_molar_density
        pH = -value(log10(m.fs.unit.outlet.mole_frac_comp[0, "H_+"]*total_molar_density))
        pOH = -value(log10(m.fs.unit.outlet.mole_frac_comp[0, "OH_-"]*total_molar_density))
        print(-value(log10(m.fs.unit.outlet.mole_frac_comp[0, "H_+"]*total_molar_density)))
        print(-value(log10(m.fs.unit.outlet.mole_frac_comp[0, "OH_-"]*total_molar_density)))
#        assert pytest.approx(6.9997414, rel=1e-5) == pH
#        assert pytest.approx(6.9997414, rel=1e-5) == pOH
#        assert pytest.approx(0.989999996, rel=1e-5) == value(m.fs.unit.outlet.mole_frac_comp[0.0, 'H2O'])
      
        print(value((m.fs.unit.outlet.mole_frac_comp[0, "Ca(OH)2"])))
        print(value((m.fs.unit.outlet.mole_frac_comp[0, "Ca_2+"])))
        assert False
