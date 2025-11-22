.. _how_to_use_inherent_reactions:

How to use inherent reactions
=============================

In :ref:`How to setup simple chemistry<how_to_setup_simple_chemistry>`, we showed
an example of setting up a **thermo_config** and **reaction_config** and put all
reactions in that **reaction_config**. However, it is possible to place reactions
into the **thermo_config** itself using **Inherent Reactions**. This guide will
perform the same walk through and setup from :ref:`How to setup simple chemistry<how_to_setup_simple_chemistry>`,
but will place reactions into the **thermo_config** itself.


Inherent Reaction vs Other Reactions
------------------------------------

The **Inherent Reactions** are reactions that will be common to all unit processes
within a flowsheet. Thus, it is convenient to put those common reactions into
the **thermo_config** as inherent. Then, the non-inherent reactions in the **reaction_config**
will be unique to a specific unit process.

.. note::

    It is possible to mix inherent reactions and non-inherent reactions within the
    same unit process. However, user's need to take care and make sure that if a
    reaction is in the **thermo_config** as inherent, then it SHOULD NOT also
    show up in the **reaction_config**. This would create a degeneracy in the model.


Example **thermo_config**
^^^^^^^^^^^^^^^^^^^^^^^^^

In this example, we will define a simple thermo-properties configuration dictionary
for a chemical system that contains only water. The water dissociation reaction will
be declared as **inherent** and, thus, be a part of this configuration dictionary.

.. testsetup::

   # quiet idaes logs
   import idaes.logger as idaeslogger
   idaeslogger.getLogger('ideas.core').setLevel('CRITICAL')
   idaeslogger.getLogger('idaes.init').setLevel('CRITICAL')

.. testcode::

    # Importing the object for units from pyomo
    from pyomo.environ import units as pyunits

    # Imports from idaes core
    from idaes.core import AqueousPhase
    from idaes.core.base.components import Solvent, Cation, Anion

    # Imports from idaes generic models
    import idaes.models.properties.modular_properties.pure.Perrys as Perrys
    from idaes.models.properties.modular_properties.pure.ConstantProperties import Constant
    from idaes.models.properties.modular_properties.state_definitions import FTPx
    from idaes.models.properties.modular_properties.eos.ideal import Ideal

    # Import the object/function for heat of reaction
    from idaes.models.properties.modular_properties.reactions.dh_rxn import constant_dh_rxn

    # Import built-in Gibb's Energy function
    from idaes.models.properties.modular_properties.reactions.equilibrium_constant import van_t_hoff

    # Import safe log power law equation
    from idaes.models.properties.modular_properties.reactions.equilibrium_forms import log_power_law_equil

    # Importing the enum for concentration unit basis used in the 'get_concentration_term' function
    from idaes.models.properties.modular_properties.base.generic_reaction import ConcentrationForm

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
                        # Parameters here come from Perry's Handbook:  p. 2-98
                        "dens_mol_liq_comp_coeff": {
                            'eqn_type': 1,
                            '1': (5.459, pyunits.kmol*pyunits.m**-3),
                            '2': (0.30542, pyunits.dimensionless),
                            '3': (647.13, pyunits.K),
                            '4': (0.081, pyunits.dimensionless)},
                        "enth_mol_form_liq_comp_ref": (-285.830, pyunits.kJ/pyunits.mol),
                        "enth_mol_form_vap_comp_ref": (0, pyunits.kJ/pyunits.mol),
                        # Parameters here come Perry's Handbook:  p. 2-174
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
                        "entr_mol_form_liq_comp_ref": (69.95, pyunits.J/pyunits.K/pyunits.mol)
                        # End parameter_data
                        }},
            'H_+': {"type": Cation, "charge": 1,
                  # Define the methods used to calculate the following properties
                  "dens_mol_liq_comp": Constant,
                  "enth_mol_liq_comp": Constant,
                  "cp_mol_liq_comp": Constant,
                  "entr_mol_liq_comp": Constant,
                  # Parameter data is always associated with the methods defined above
                  "parameter_data": {
                        "mw": (1.00784, pyunits.g/pyunits.mol),
                        "dens_mol_liq_comp_coeff": (55, pyunits.kmol*pyunits.m**-3),
                        "enth_mol_form_liq_comp_ref": (0, pyunits.kJ/pyunits.mol),
                        "cp_mol_liq_comp_coeff": (75000, pyunits.J/pyunits.kmol/pyunits.K),
                        "entr_mol_form_liq_comp_ref": (0, pyunits.J/pyunits.K/pyunits.mol)
                                    },
                        # End parameter_data
                        },
            'OH_-': {"type": Anion, "charge": -1,
                  # Define the methods used to calculate the following properties
                  "dens_mol_liq_comp": Constant,
                  "enth_mol_liq_comp": Constant,
                  "cp_mol_liq_comp": Constant,
                  "entr_mol_liq_comp": Constant,
                  # Parameter data is always associated with the methods defined above
                  "parameter_data": {
                        "mw": (17.008, pyunits.g/pyunits.mol),
                        "dens_mol_liq_comp_coeff": (55, pyunits.kmol*pyunits.m**-3),
                        "enth_mol_form_liq_comp_ref": (-230.000, pyunits.kJ/pyunits.mol),
                        "cp_mol_liq_comp_coeff": (75000, pyunits.J/pyunits.kmol/pyunits.K),
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

            # This is an optional dictionary to setup bounds on
            #   the state variables. Names below MUST correspond
            #   to the 'FTPx' type state definition
            "state_bounds": {"flow_mol": (0, 50, 100),
                             "temperature": (273.15, 300, 650),
                             "pressure": (5e4, 1e5, 1e6)
                         },

            # These are generally optional parameters, however, because we
            #   are using the Perry's model to calculate temperature dependent
            #   properties, we MUST provide these here.
            "pressure_ref": 1e5,
            "temperature_ref": 300,

            # Our dictionary for base units MUST define the following
            "base_units": {"time": pyunits.s,
                           "length": pyunits.m,
                           "mass": pyunits.kg,
                           "amount": pyunits.mol,
                           "temperature": pyunits.K},

             # Inherent reactions
             #    These are added just like any other equilibrium reaction
             #    would be defined in a reaction config
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
                  }
                  # End inherent reactions
        }
        # End thermo_config definition

For a detailed analysis of everything from above, see
:ref:`How to setup simple chemistry<how_to_setup_simple_chemistry>`.


.. note::

    Even if your system only involves inherent reactions, you may still be required to provide
    a reaction configuration dictionary to construct certain unit models (such as EquilibriumReactor).
    However, your reaction configuration may be a blank or a **dummy** config (see below)

Example of a dummy **reaction_config**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. testcode::

    # Importing the object for units from pyomo
    from pyomo.environ import units as pyunits

    # Import safe log power law equation
    from idaes.models.properties.modular_properties.reactions.equilibrium_forms import log_power_law_equil

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
                   # End reaction
             }
             # End equilibrium_reactions
        }
        # End reaction_config definition


Example: Using our configuration dictionaries in an EquilibriumReactor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Recall, we had named our configuration dictionaries as ``thermo_config`` and
``reaction_config``. We will reference those dictionary names in the example
code below.

.. testcode::
    
    # Import specific pyomo objects
    from pyomo.environ import ConcreteModel

    # Import the core idaes objects for Flowsheets and types of balances
    from idaes.core import FlowsheetBlock

    # Import the idaes objects for Generic Properties and Reactions
    from idaes.models.properties.modular_properties.base.generic_property import (
        GenericParameterBlock,
    )
    from idaes.models.properties.modular_properties.base.generic_reaction import (
        GenericReactionParameterBlock,
    )

    # Import the idaes object for the EquilibriumReactor unit model
    from idaes.models.unit_models.equilibrium_reactor import EquilibriumReactor

    # Create an instance of a pyomo model
    model = ConcreteModel()

    # Add an IDAES flowsheet to that model
    model.fs = FlowsheetBlock(dynamic=False)

    # Add a thermo parameter block to that flowsheet
    #   Here, we are passing our 'thermo_config' dictionary we created earlier
    model.fs.thermo_params = GenericParameterBlock(**thermo_config)

    # Add a reaction parameter block to that flowsheet
    #   Here, we are passing our thermo block created above as the property package
    #   and then giving our 'reaction_config' as the instructions for how the
    #   reactions will be constructed from the thermo package.
    model.fs.rxn_params = GenericReactionParameterBlock(
        property_package=model.fs.thermo_params, **reaction_config
    )

    # Add an EquilibriumReactor object as the unit model
    #   Here, we pass both the thermo package and reaction package, as well
    #   as a number of other arguments to help define how this unit process
    #   will behave.
    #
    # NOTE: What is different here is now we state that there are no
    #       equilibrium reactions in this unit model because we defined
    #       those reactions as inherent.
    model.fs.unit = EquilibriumReactor(
        property_package=model.fs.thermo_params,
        reaction_package=model.fs.rxn_params,
        has_rate_reactions=False,
        has_equilibrium_reactions=False,
        has_heat_transfer=False,
        has_heat_of_reaction=False,
        has_pressure_change=False,
    )

    # At this point, you can 'fix' your inlet/outlet state conditions,
    #     setup scaling factors, initialize the model, then solve the model
    #     just as you would with any other IDAES flowsheet
