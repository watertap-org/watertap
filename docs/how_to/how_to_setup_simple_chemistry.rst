How to setup simple chemistry
=============================

.. _GenericProperties: https://idaes-pse.readthedocs.io/en/stable/user_guide/components/property_package/general/index.html#generic-property-package-framework
.. _GenericReactions: https://idaes-pse.readthedocs.io/en/stable/user_guide/components/property_package/general_reactions/index.html

In ProteusLib, chemistry modules leverage the Generic Properties
(`GenericProperties`_)
and Generic Reactions
(`GenericReactions`_)
objects in IDAES. These objects can be used in conjunction with any unit process
where chemical reactions need to be considered. In this guide, we will cover how
to use these built-in objects within your own unit process.

What you will need
------------------

1. [Required] A **thermo-properties** configuration dictionary
2. [Optional] A **reaction-properties** configuration dictionary
3. [Required] A **unit model** upon which to build the chemistry module

.. note::

    The **reaction-properties** dictionary is optional for certain unit models and
    required for others. Additionally, certain reactions can actually be defined
    in the **thermo-properties** dictionary and would therefore NOT be included in
    the **reaction-properties** dictionary. The differences will be covered in another
    how-to guide.

The **thermo-properties** configuration dictionary
--------------------------------------------------

You will always be required to provide a configuration dictionary for thermodynamic
properties of the chemical species of interest in your process. At a minimum, the
**thermo-properties** dictionary MUST contain the following keys...

+----------------------+-------------------------------------------------------------------------------------------+
|     Key              |  Description                                                                              |
+======================+===========================================================================================+
| components           | dictionary containing the chemical species of interest and their properties               |
+----------------------+-------------------------------------------------------------------------------------------+
| phases               | dictionary containing the phases (gas, liquid, etc) for the system and its species        |
+----------------------+-------------------------------------------------------------------------------------------+
| state_definition     | IDAES object which defines how inlet/outlet ports define the initial state of the system  |
+----------------------+-------------------------------------------------------------------------------------------+
| base_units           | dictionary containing the base units the model will use for all equations in the system   |
+----------------------+-------------------------------------------------------------------------------------------+

.. note::

    Depending on your system, there may be additional keys in the **thermo-properties**
    configuration dictionary that may be required. For instance, multi-phase problems
    also require you to define how the model should represent the phase equilibria
    between the system's phases. An in depth discussion of all options and methods
    is beyond the scope of this guide. For additional information, refer to the IDAES
    documentation (`GenericProperties`_).

Example thermo-properties configuration dictionary
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this example, we will define a simple thermo-properties configuration dictionary
for a chemical system that contains only water.

.. code-block:: python

    # Configuration dictionary
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
                       "equilibrium_constant": gibbs_energy,
                       "equilibrium_form": log_power_law_equil,
                       "concentration_form": ConcentrationForm.molarity,
                       "parameter_data": {
                           "dh_rxn_ref": (55.830, pyunits.kJ/pyunits.mol),
                           "ds_rxn_ref": (-80.7, pyunits.J/pyunits.mol/pyunits.K),
                           "T_eq_ref": (300, pyunits.K),

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
