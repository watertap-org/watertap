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

"""
    This test is to establish that the core chemistry packages in IDAES solve
    a simple water dissociation problem and return the correct pH value.
"""
import enum

# Importing testing libraries
import pytest

# Importing the object for units from pyomo
from pyomo.environ import units as pyunits

# Imports from idaes core
from idaes.core import AqueousPhase
from idaes.core.base.components import Solvent, Solute, Cation, Anion
from idaes.core.base.phases import PhaseType as PT

# Imports from idaes generic models
import idaes.models.properties.modular_properties.pure.Perrys as Perrys
from idaes.models.properties.modular_properties.state_definitions import FTPx
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
    gibbs_energy,
    van_t_hoff,
)

# Import specific pyomo objects
from pyomo.environ import (
    ConcreteModel,
    SolverStatus,
    TerminationCondition,
    value,
    Suffix,
)

from idaes.core.util import scaling as iscale

# Import pyomo methods to check the system units
from pyomo.util.check_units import assert_units_consistent

# Import idaes methods to check the model during construction
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    fixed_variables_set,
    activated_constraints_set,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)

# Import the idaes objects for Generic Properties and Reactions
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock,
)

# Import the idaes object for the EquilibriumReactor unit model
from idaes.models.unit_models.equilibrium_reactor import EquilibriumReactor

# Import the core idaes objects for Flowsheets and types of balances
from idaes.core import FlowsheetBlock

# Import log10 function from pyomo
from pyomo.environ import log10


__authors__ = [
    "Austin Ladshaw",
    "Ludovico Bianchi",
    "Dan Gunter",
]
__author__ = __authors__[0]


class Variant(str, enum.Enum):
    equilibrium = "equilibrium"
    inherent = "inherent"

    def __str__(self):
        return f"{self.__class__.__name__}.{self.name}"

    @property
    def is_equilibrium(self):
        return self is Variant.equilibrium


@pytest.fixture(scope="module", params=[Variant.equilibrium, Variant.inherent], ids=str)
def variant(request) -> Variant:
    return request.param


@pytest.fixture(scope="module")
def base_units():
    return {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    }


@pytest.fixture(scope="module")
def solver():
    s = get_solver()
    s.options["max_iter"] = 200
    return s


@pytest.fixture(scope="module")
def thermo_config(base_units):
    # Configuration dictionary
    data = {
        "components": {
            "H2O": {
                "type": Solvent,
                # Define the methods used to calculate the following properties
                "dens_mol_liq_comp": Perrys,
                "enth_mol_liq_comp": Perrys,
                "cp_mol_liq_comp": Perrys,
                "entr_mol_liq_comp": Perrys,
                # Parameter data is always associated with the methods defined above
                "parameter_data": {
                    "mw": (18.0153, pyunits.g / pyunits.mol),
                    "pressure_crit": (220.64e5, pyunits.Pa),
                    "temperature_crit": (647, pyunits.K),
                    # Comes from Perry's Handbook:  p. 2-98
                    "dens_mol_liq_comp_coeff": {
                        "1": (5.459, pyunits.kmol * pyunits.m**-3),
                        "2": (0.30542, pyunits.dimensionless),
                        "3": (647.13, pyunits.K),
                        "4": (0.081, pyunits.dimensionless),
                    },
                    "enth_mol_form_liq_comp_ref": (-285.830, pyunits.kJ / pyunits.mol),
                    "enth_mol_form_vap_comp_ref": (0, pyunits.kJ / pyunits.mol),
                    # Comes from Perry's Handbook:  p. 2-174
                    "cp_mol_liq_comp_coeff": {
                        "1": (2.7637e5, pyunits.J / pyunits.kmol / pyunits.K),
                        "2": (-2.0901e3, pyunits.J / pyunits.kmol / pyunits.K**2),
                        "3": (8.125, pyunits.J / pyunits.kmol / pyunits.K**3),
                        "4": (-1.4116e-2, pyunits.J / pyunits.kmol / pyunits.K**4),
                        "5": (9.3701e-6, pyunits.J / pyunits.kmol / pyunits.K**5),
                    },
                    "cp_mol_ig_comp_coeff": {
                        "A": (30.09200, pyunits.J / pyunits.mol / pyunits.K),
                        "B": (
                            6.832514,
                            pyunits.J
                            * pyunits.mol**-1
                            * pyunits.K**-1
                            * pyunits.kiloK**-1,
                        ),
                        "C": (
                            6.793435,
                            pyunits.J
                            * pyunits.mol**-1
                            * pyunits.K**-1
                            * pyunits.kiloK**-2,
                        ),
                        "D": (
                            -2.534480,
                            pyunits.J
                            * pyunits.mol**-1
                            * pyunits.K**-1
                            * pyunits.kiloK**-3,
                        ),
                        "E": (
                            0.082139,
                            pyunits.J
                            * pyunits.mol**-1
                            * pyunits.K**-1
                            * pyunits.kiloK**2,
                        ),
                        "F": (-250.8810, pyunits.kJ / pyunits.mol),
                        "G": (223.3967, pyunits.J / pyunits.mol / pyunits.K),
                        "H": (0, pyunits.kJ / pyunits.mol),
                    },
                    "entr_mol_form_liq_comp_ref": (
                        69.95,
                        pyunits.J / pyunits.K / pyunits.mol,
                    ),
                    "pressure_sat_comp_coeff": {
                        "A": (4.6543, None),  # [1], temperature range 255.9 K - 373 K
                        "B": (1435.264, pyunits.K),
                        "C": (-64.848, pyunits.K),
                    },
                },
                # End parameter_data
            },
            "H_+": {
                "type": Cation,
                "charge": 1,
                # Define the methods used to calculate the following properties
                "dens_mol_liq_comp": Perrys,
                "enth_mol_liq_comp": Perrys,
                "cp_mol_liq_comp": Perrys,
                "entr_mol_liq_comp": Perrys,
                # Parameter data is always associated with the methods defined above
                "parameter_data": {
                    "mw": (1.00784, pyunits.g / pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        "1": (5.459, pyunits.kmol * pyunits.m**-3),
                        "2": (0.30542, pyunits.dimensionless),
                        "3": (647.13, pyunits.K),
                        "4": (0.081, pyunits.dimensionless),
                    },
                    "enth_mol_form_liq_comp_ref": (-230.000, pyunits.kJ / pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        "1": (2.7637e5, pyunits.J / pyunits.kmol / pyunits.K),
                        "2": (-2.0901e3, pyunits.J / pyunits.kmol / pyunits.K**2),
                        "3": (8.125, pyunits.J / pyunits.kmol / pyunits.K**3),
                        "4": (-1.4116e-2, pyunits.J / pyunits.kmol / pyunits.K**4),
                        "5": (9.3701e-6, pyunits.J / pyunits.kmol / pyunits.K**5),
                    },
                    "entr_mol_form_liq_comp_ref": (
                        -10.75,
                        pyunits.J / pyunits.K / pyunits.mol,
                    ),
                },
                # End parameter_data
            },
            "OH_-": {
                "type": Anion,
                "charge": -1,
                # Define the methods used to calculate the following properties
                "dens_mol_liq_comp": Perrys,
                "enth_mol_liq_comp": Perrys,
                "cp_mol_liq_comp": Perrys,
                "entr_mol_liq_comp": Perrys,
                # Parameter data is always associated with the methods defined above
                "parameter_data": {
                    "mw": (17.008, pyunits.g / pyunits.mol),
                    "dens_mol_liq_comp_coeff": {
                        "1": (5.459, pyunits.kmol * pyunits.m**-3),
                        "2": (0.30542, pyunits.dimensionless),
                        "3": (647.13, pyunits.K),
                        "4": (0.081, pyunits.dimensionless),
                    },
                    "enth_mol_form_liq_comp_ref": (-230.000, pyunits.kJ / pyunits.mol),
                    "cp_mol_liq_comp_coeff": {
                        "1": (2.7637e5, pyunits.J / pyunits.kmol / pyunits.K),
                        "2": (-2.0901e3, pyunits.J / pyunits.kmol / pyunits.K**2),
                        "3": (8.125, pyunits.J / pyunits.kmol / pyunits.K**3),
                        "4": (-1.4116e-2, pyunits.J / pyunits.kmol / pyunits.K**4),
                        "5": (9.3701e-6, pyunits.J / pyunits.kmol / pyunits.K**5),
                    },
                    "entr_mol_form_liq_comp_ref": (
                        -10.75,
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
        "state_definition": FTPx,
        "state_bounds": {
            "flow_mol": (0, 50, 100),
            "temperature": (273.15, 300, 650),
            "pressure": (5e4, 1e5, 1e6),
        },
        "pressure_ref": 1e5,
        "temperature_ref": 300,
        "base_units": base_units,
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
    return data


# Define the reaction_config for water dissociation
@pytest.fixture(scope="module")
def water_reaction_config(base_units):
    return {
        "base_units": base_units,
        "equilibrium_reactions": {
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
        }
        # End equilibrium_reactions
    }
    # End reaction_config definition


def _get_without_inherent_reactions(
    config_dict: dict, key="inherent_reactions"
) -> dict:
    data = dict(config_dict)
    del data[key]
    return data


class TestPureWater:
    @pytest.fixture
    def model(self, thermo_config, variant: Variant, water_reaction_config):
        if variant.is_equilibrium:
            thermo_config = _get_without_inherent_reactions(thermo_config)
        model = ConcreteModel()
        model.fs = FlowsheetBlock(dynamic=False)
        model.fs.thermo_params = GenericParameterBlock(**thermo_config)
        print(water_reaction_config)
        model.fs.rxn_params = GenericReactionParameterBlock(
            property_package=model.fs.thermo_params, **water_reaction_config
        )
        model.fs.unit = EquilibriumReactor(
            property_package=model.fs.thermo_params,
            reaction_package=model.fs.rxn_params,
            has_rate_reactions=False,
            has_equilibrium_reactions=variant.is_equilibrium,
            has_heat_transfer=False,
            has_heat_of_reaction=False,
            has_pressure_change=False,
        )

        model.fs.unit.inlet.mole_frac_comp[0, "H_+"].fix(0.0)
        model.fs.unit.inlet.mole_frac_comp[0, "OH_-"].fix(0.0)
        model.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix(1.0)
        model.fs.unit.inlet.pressure.fix(101325.0)
        model.fs.unit.inlet.temperature.fix(298.0)
        model.fs.unit.inlet.flow_mol.fix(10)

        return model

    @pytest.mark.unit
    def test_build_model(self, model):

        assert hasattr(model.fs.thermo_params, "component_list")
        assert len(model.fs.thermo_params.component_list) == 3
        assert "H2O" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.H2O, Solvent)
        assert "H_+" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("H_+"), Cation)
        assert "OH_-" in model.fs.thermo_params.component_list
        assert isinstance(model.fs.thermo_params.component("OH_-"), Anion)

        assert hasattr(model.fs.thermo_params, "phase_list")
        assert len(model.fs.thermo_params.phase_list) == 1
        assert isinstance(model.fs.thermo_params.Liq, AqueousPhase)

    @pytest.mark.unit
    def test_units_inherent(self, model):
        assert_units_consistent(model)

    @pytest.mark.unit
    def test_dof(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.fixture
    def expected_stats(self, variant):
        expected = {
            Variant.inherent: {
                number_variables: 81,
                number_total_constraints: 28,
                number_unused_variables: 12,
            },
            Variant.equilibrium: {
                number_variables: 75,
                number_total_constraints: 28,
                number_unused_variables: 6,
            },
        }
        return expected[variant]

    @pytest.mark.unit
    def test_stats(self, model, expected_stats):
        for check_func, expected_result in expected_stats.items():
            check_result = check_func(model)
            assert check_result == expected_result

    def _do_inherent_reaction_scaling(self, model):
        for i in model.fs.unit.control_volume.inherent_reaction_extent_index:
            scale = value(
                model.fs.unit.control_volume.properties_out[0.0].k_eq[i[1]].expr
            )
            iscale.set_scaling_factor(
                model.fs.unit.control_volume.inherent_reaction_extent[0.0, i[1]],
                10 / scale,
            )
            iscale.constraint_scaling_transform(
                model.fs.unit.control_volume.properties_out[
                    0.0
                ].inherent_equilibrium_constraint[i[1]],
                0.1,
            )
        return model

    def _do_equilibrium_reaction_scaling(self, model):
        for i in model.fs.unit.control_volume.equilibrium_reaction_extent_index:
            scale = value(model.fs.unit.control_volume.reactions[0.0].k_eq[i[1]].expr)
            iscale.set_scaling_factor(
                model.fs.unit.control_volume.equilibrium_reaction_extent[0.0, i[1]],
                10 / scale,
            )
            iscale.constraint_scaling_transform(
                model.fs.unit.control_volume.reactions[0.0].equilibrium_constraint[
                    i[1]
                ],
                0.1,
            )
        return model

    @pytest.fixture
    def model_scaling(self, model, variant: Variant):
        map_variant_reaction_scale_func = {
            Variant.equilibrium: self._do_equilibrium_reaction_scaling,
            Variant.inherent: self._do_inherent_reaction_scaling,
        }

        scale_func = map_variant_reaction_scale_func[variant]
        scale_func(model)

        # Next, try adding scaling for species
        min = 1e-3
        for i in model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp:
            # i[0] = phase, i[1] = species
            if model.fs.unit.inlet.mole_frac_comp[0, i[1]].value > min:
                scale = model.fs.unit.inlet.mole_frac_comp[0, i[1]].value
            else:
                scale = min
            iscale.set_scaling_factor(
                model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp[i[1]],
                10 / scale,
            )
            iscale.set_scaling_factor(
                model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[
                    i
                ],
                10 / scale,
            )
            iscale.set_scaling_factor(
                model.fs.unit.control_volume.properties_out[0.0].flow_mol_phase_comp[i],
                10 / scale,
            )
            iscale.constraint_scaling_transform(
                model.fs.unit.control_volume.properties_out[
                    0.0
                ].component_flow_balances[i[1]],
                10 / scale,
            )
            iscale.constraint_scaling_transform(
                model.fs.unit.control_volume.material_balances[0.0, i[1]], 10 / scale
            )

        iscale.calculate_scaling_factors(model.fs.unit)
        return model

    @pytest.mark.component
    def test_scaling(self, model_scaling, variant: Variant):
        model = model_scaling

        assert isinstance(model.fs.unit.control_volume.scaling_factor, Suffix)

        assert isinstance(
            model.fs.unit.control_volume.properties_out[0.0].scaling_factor, Suffix
        )

        assert isinstance(
            model.fs.unit.control_volume.properties_in[0.0].scaling_factor, Suffix
        )

        if variant.is_equilibrium:
            # When using equilibrium reactions, there is another set of scaling factors calculated
            assert isinstance(
                model.fs.unit.control_volume.reactions[0.0].scaling_factor, Suffix
            )

    @pytest.fixture
    def state_args(self):
        return {
            "mole_frac_comp": {
                "H2O": 1,
                "H_+": 10**-7 / 55.6,
                "OH_-": 10**-7 / 55.6,
            },
            "pressure": 101325,
            "temperature": 298,
            "flow_mol": 10,
        }

    @pytest.fixture
    def model_initialize(self, model, state_args, solver):
        def _collect_data_to_check(m):
            return {
                fixed_variables_set: fixed_variables_set(m),
                activated_constraints_set: activated_constraints_set(m),
            }

        data_before = _collect_data_to_check(model)
        model.fs.unit.initialize(state_args=state_args, optarg=solver.options)

        data_after = _collect_data_to_check(model)
        return model, data_before, data_after

    @pytest.mark.component
    def test_initialize(self, model_initialize):
        model, data_before, data_after = model_initialize

        assert degrees_of_freedom(model) == 0

        for key, value_before in data_before.items():
            value_after = data_after[key]
            assert len(value_before) == len(value_after)

    @pytest.fixture
    def model_solve(self, model, solver):
        results = solver.solve(model, tee=True)
        return model, results

    @pytest.mark.component
    def test_solve(self, model_solve):
        _, results = model_solve
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_solution(self, model_solve):
        model, _ = model_solve

        assert pytest.approx(298, rel=1e-5) == value(
            model.fs.unit.outlet.temperature[0]
        )
        assert pytest.approx(10, rel=1e-5) == value(model.fs.unit.outlet.flow_mol[0])
        assert pytest.approx(101325, rel=1e-5) == value(
            model.fs.unit.outlet.pressure[0]
        )

        total_molar_density = (
            value(
                model.fs.unit.control_volume.properties_out[0.0].dens_mol_phase["Liq"]
            )
            / 1000
        )
        assert pytest.approx(55.2336, rel=1e-5) == total_molar_density

        pH = -value(
            log10(model.fs.unit.outlet.mole_frac_comp[0, "H_+"] * total_molar_density)
        )
        pOH = -value(
            log10(model.fs.unit.outlet.mole_frac_comp[0, "OH_-"] * total_molar_density)
        )
        assert pytest.approx(6.9997414, rel=1e-5) == pH
        assert pytest.approx(6.9997414, rel=1e-5) == pOH
        assert pytest.approx(0.99999, rel=1e-5) == value(
            model.fs.unit.outlet.mole_frac_comp[0.0, "H2O"]
        )


class TestPureWaterEDB(TestPureWater):
    @pytest.fixture
    def thermo_config(self, edb):
        base = edb.get_base("default_thermo")
        elements = ["H", "O"]
        components = []
        # Add the components
        for c in edb.get_components(element_names=elements):
            # Should remove these since we don't have a vapor phase
            # NOTE: The model still runs without removing these, but
            #       the 'stats' won't pass their checks due to the addition
            #       of new system variables that they inherently bring
            c.remove("enth_mol_ig_comp")
            c.remove("pressure_sat_comp")
            base.add(c)
            components.append(c.name)
        # Add the reactions
        for r in edb.get_reactions(component_names=components):
            # Set a custom reaction order
            r.set_reaction_order("Liq", {"H2O": 0, "H_+": 1, "OH_-": 1})
            r._data["type"] = "inherent"
            base.add(r)
        return base.idaes_config

    @pytest.fixture
    def water_reaction_config(self, edb):
        elements = ["H", "O"]
        components = [c.name for c in edb.get_components(element_names=elements)]
        base = edb.get_base("reaction")
        # Add the reactions
        for r in edb.get_reactions(component_names=components):
            # Set a custom reaction order
            r.set_reaction_order("Liq", {"H2O": 0, "H_+": 1, "OH_-": 1})
            # Need to remove this to avoid errors when using the generated config
            r.remove_parameter("ds_rxn_ref")
            base.add(r)
        return base.idaes_config
