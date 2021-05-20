"""
Shared input data for tests
"""
from pyomo.environ import units as pyunits
from idaes.core.phases import PhaseType as PT
import idaes.generic_models.properties.core.pure.Perrys as Perrys

reaction_data = [
    {
        "stoichiometry": {"Liq/NH4 +": -1, "Liq/H +": 1, "Liq/NH3": 1},
        "heat_of_reaction": "constant_dh_rxn",
        "equilibrium_constant": "van_t_hoff_aqueous",
        "equilibrium_form": "log_power_law",
        "concentration_form": "ConcentrationForm.molarity",
        "parameter_data": {
            "dh_rxn_ref": [52.21, "U.kJ/U.mol"],
            "ds_rxn_ref": [-2.4, "U.J/U.mol/U.K"],
            "reaction_order": {"Liq/NH4 +": -1, "Liq/H +": 1, "Liq/NH3": 1},
        },
        "type": "reaction",
        "reaction_type": "equilibrium",
        "name": "NH4_Ka",
        "components": ["NH4", "Ka"],
        "base_units": {
            "time": "U.s",
            "length": "U.m",
            "mass": "U.kg",
            "amount": "U.mol",
            "temperature": "U.K",
        },
    },
    {
        "stoichiometry": {"Liq/HCO3 -": -1, "Liq/H +": 1, "Liq/CO3 2-": 1},
        "heat_of_reaction": "constant_dh_rxn",
        "equilibrium_constant": "van_t_hoff_aqueous",
        "equilibrium_form": "log_power_law",
        "concentration_form": "ConcentrationForm.molarity",
        "parameter_data": {
            "dh_rxn_ref": [14.9, "U.kJ/U.mol"],
            "ds_rxn_ref": [-148.1, "U.J/U.mol/U.K"],
            "reaction_order": {"Liq/HCO3 -": -1, "Liq/H +": 1, "Liq/CO3 2-": 1},
        },
        "type": "reaction",
        "reaction_type": "equilibrium",
        "name": "H2CO3_Ka2",
        "components": ["H2CO3", "Ka2"],
        "base_units": {
            "time": "U.s",
            "length": "U.m",
            "mass": "U.kg",
            "amount": "U.mol",
            "temperature": "U.K",
        },
    },
]
component_data = [
    {
        "type": "component",
        "valid_phase_types": "PT.liquidPhase",
        "dens_mol_liq_comp": "Perrys",
        "enth_mol_liq_comp": "Perrys",
        "cp_mol_liq_comp": "Perrys",
        "entr_mol_liq_comp": "Perrys",
        "parameter_data": {
            "mw": [74.09, "U.g/U.mol"],
            "dens_mol_liq_comp_coeff": {
                "1": [13.5, "U.kmol*U.m**-3"],
                "2": [1, "U.dimensionless"],
                "3": [1, "U.K"],
                "4": [1, "U.dimensionless"],
            },
            "enth_mol_form_liq_comp_ref": [-1003, "U.kJ/U.mol"],
            "cp_mol_liq_comp_coeff": {
                "1": [276370.0, "U.J/U.kmol/U.K"],
                "2": [-2090.1, "U.J/U.kmol/U.K**2"],
                "3": [8.125, "U.J/U.kmol/U.K**3"],
                "4": [-0.014116, "U.J/U.kmol/U.K**4"],
                "5": [9.3701e-06, "U.J/U.kmol/U.K**5"],
            },
            "entr_mol_form_liq_comp_ref": [-74.5, "U.J/U.K/U.mol"],
        },
        "name": "Ca[OH]2",
    },
    {
        "type": "component",
        "valid_phase_types": "PT.liquidPhase",
        "dens_mol_liq_comp": "Perrys",
        "enth_mol_liq_comp": "Perrys",
        "cp_mol_liq_comp": "Perrys",
        "entr_mol_liq_comp": "Perrys",
        "parameter_data": {
            "mw": [57.082, "U.g/U.mol"],
            "dens_mol_liq_comp_coeff": {
                "1": [13.5, "U.kmol*U.m**-3"],
                "2": [1, "U.dimensionless"],
                "3": [1, "U.K"],
                "4": [1, "U.dimensionless"],
            },
            "enth_mol_form_liq_comp_ref": [-756.18, "U.kJ/U.mol"],
            "cp_mol_liq_comp_coeff": {
                "1": [276370.0, "U.J/U.kmol/U.K"],
                "2": [-2090.1, "U.J/U.kmol/U.K**2"],
                "3": [8.125, "U.J/U.kmol/U.K**3"],
                "4": [-0.014116, "U.J/U.kmol/U.K**4"],
                "5": [9.3701e-06, "U.J/U.kmol/U.K**5"],
            },
            "entr_mol_form_liq_comp_ref": [16.95, "U.J/U.K/U.mol"],
        },
        "name": "CaOH +",
    },
    {
        "type": "component",
        "valid_phase_types": "PT.liquidPhase",
        "dens_mol_liq_comp": "Perrys",
        "enth_mol_liq_comp": "Perrys",
        "cp_mol_liq_comp": "Perrys",
        "entr_mol_liq_comp": "Perrys",
        "parameter_data": {
            "mw": [40.078, "U.g/U.mol"],
            "dens_mol_liq_comp_coeff": {
                "1": [13.5, "U.kmol*U.m**-3"],
                "2": [1, "U.dimensionless"],
                "3": [1, "U.K"],
                "4": [1, "U.dimensionless"],
            },
            "enth_mol_form_liq_comp_ref": [-542.83, "U.J/U.mol"],
            "cp_mol_liq_comp_coeff": {
                "1": [276370.0, "U.J/U.kmol/U.K"],
                "2": [-2090.1, "U.J/U.kmol/U.K**2"],
                "3": [8.125, "U.J/U.kmol/U.K**3"],
                "4": [-0.014116, "U.J/U.kmol/U.K**4"],
                "5": [9.3701e-06, "U.J/U.kmol/U.K**5"],
            },
            "entr_mol_form_liq_comp_ref": [-53, "U.J/U.K/U.mol"],
        },
        "name": "Ca 2+",
    },
]

Ca_thermo_config = {
    "components": {
        'Ca 2+': {"type": "Component",
                  "valid_phase_types": PT.liquidPhase,
                  # Define the methods used to calculate the following properties
                  "dens_mol_liq_comp": Perrys,
                  "enth_mol_liq_comp": Perrys,
                  "cp_mol_liq_comp": Perrys,
                  "entr_mol_liq_comp": Perrys,
                  # Parameter data is always associated with the methods defined above
                  "parameter_data": {
                      "mw": (40.078, pyunits.g / pyunits.mol),
                      "dens_mol_liq_comp_coeff": {
                          '1': (13.5, pyunits.kmol * pyunits.m ** -3),
                          '2': (1, pyunits.dimensionless),
                          '3': (1, pyunits.K),
                          '4': (1, pyunits.dimensionless)},
                      "enth_mol_form_liq_comp_ref": (-542.83, pyunits.J / pyunits.mol),
                      "cp_mol_liq_comp_coeff": {
                          '1': (2.7637E5, pyunits.J / pyunits.kmol / pyunits.K),
                          '2': (-2.0901E3, pyunits.J / pyunits.kmol / pyunits.K ** 2),
                          '3': (8.125, pyunits.J / pyunits.kmol / pyunits.K ** 3),
                          '4': (-1.4116E-2, pyunits.J / pyunits.kmol / pyunits.K ** 4),
                          '5': (9.3701E-6, pyunits.J / pyunits.kmol / pyunits.K ** 5)},
                      "entr_mol_form_liq_comp_ref": (-53, pyunits.J / pyunits.K / pyunits.mol)
                  },
                  # End parameter_data
                  }
    }
    # End Component list
}
# End thermo_config definition

Ca_thermo_data = {"_id": {"$oid": "6099cd12507b1e55f1707555"}, "type": "component",
                  "valid_phase_types": "PT.liquidPhase", "dens_mol_liq_comp": "Perrys", "enth_mol_liq_comp": "Perrys",
                  "cp_mol_liq_comp": "Perrys", "entr_mol_liq_comp": "Perrys",
                  "parameter_data": {"mw": [40.078, "U.g/U.mol"],
                                     "dens_mol_liq_comp_coeff": [[13.5, "U.kmol*U.m**-3"], [1, "U.dimensionless"],
                                                                 [1, "U.K"], [1, "U.dimensionless"]],
                                     "enth_mol_form_liq_comp_ref": [-542.83, "U.J/U.mol"],
                                     "cp_mol_liq_comp_coeff": [[276370, "U.J/U.kmol/U.K"],
                                                               [-2090.1, "U.J/U.kmol/U.K**2"],
                                                               [8.125, "U.J/U.kmol/U.K**3"],
                                                               [-0.014116, "U.J/U.kmol/U.K**4"],
                                                               [0.0000093701, "U.J/U.kmol/U.K**5"]],
                                     "entr_mol_form_liq_comp_ref": [-53, "U.J/U.K/U.mol"]}, "name": "Ca 2+",
                  "elements": ["Ca"]}
