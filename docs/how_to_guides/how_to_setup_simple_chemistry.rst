.. _how_to_setup_simple_chemistry:

How to setup simple chemistry
=============================

.. note:: This page provides a manual approach to building an IDAES configuration dictionary.
    The same result can be achieved by using the :doc:`Electrolyte Database (EDB)</technical_reference/edb/index>`. 
    Examples of this are provided under :ref:`How to use EDB <how_to_use_edb>`.
.. _GenericProperties: https://idaes-pse.readthedocs.io/en/stable/user_guide/components/property_package/general/index.html#generic-property-package-framework
.. _GenericReactions: https://idaes-pse.readthedocs.io/en/stable/user_guide/components/property_package/general_reactions/index.html
.. _Perrys: https://idaes-pse.readthedocs.io/en/stable/user_guide/components/property_package/general/pure/Perrys.html
.. _Constant: https://idaes-pse.readthedocs.io/en/stable/user_guide/components/property_package/general/pure/ConstantProperties.html
.. _StateDefinition: https://idaes-pse.readthedocs.io/en/stable/user_guide/components/property_package/general/state_definition.html
.. _EquationOfState: https://idaes-pse.readthedocs.io/en/stable/user_guide/components/property_package/general/eos/ideal.html
.. _Components: https://idaes-pse.readthedocs.io/en/stable/user_guide/components/property_package/general/component_def.html
.. _Phases: https://idaes-pse.readthedocs.io/en/stable/user_guide/components/property_package/general/phase_def.html
.. _RateReactions: https://idaes-pse.readthedocs.io/en/stable/user_guide/components/property_package/general_reactions/rate_rxns.html
.. _EquilibriumReactions: https://idaes-pse.readthedocs.io/en/stable/user_guide/components/property_package/general_reactions/equil_rxns.html
.. _ReactionMethods: https://idaes-pse.readthedocs.io/en/stable/user_guide/components/property_package/general_reactions/method_libraries.html#reaction-module-libraries
.. _ConcentrationForm: https://idaes-pse.readthedocs.io/en/stable/user_guide/components/property_package/general_reactions/rate_rxns.html#concentration-form
.. _UnitModels: https://idaes-pse.readthedocs.io/en/stable/technical_specs/model_libraries/generic/unit_models/index.html
.. _EquilibriumReactor: https://idaes-pse.readthedocs.io/en/stable/technical_specs/model_libraries/generic/unit_models/equilibrium.html
.. _IDAESWorkflow: https://idaes-pse.readthedocs.io/en/stable/user_guide/workflow/general.html

In WaterTAP, chemistry modules leverage the Generic Properties
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
    how-to guide (see :ref:`How to use inherent reactions<how_to_use_inherent_reactions>`)


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
    from idaes.core.components import Solvent, Cation, Anion

    # Imports from idaes generic models
    import idaes.generic_models.properties.core.pure.Perrys as Perrys
    from idaes.generic_models.properties.core.pure.ConstantProperties import Constant
    from idaes.generic_models.properties.core.state_definitions import FTPx
    from idaes.generic_models.properties.core.eos.ideal import Ideal

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
        }
        # End thermo_config definition

There is a significant amount to breakdown here, so let's discuss a couple of things
step by step...

**(1)** All components need a ``"type"``. For this, you have a number of ``"type"`` options within IDAES.
Generally, the ``"H2O"`` component should always be a ``Solvent`` within WaterTAP. Charged species
will always be either a ``Cation`` or ``Anion`` depending on the sign of their actual ``"charge"``.
More information on ``"components"`` can be found at `Components`_.

**(2)** All components need to have methods defined for calculating thermodynamic properties such as
``"dens_mol_liq_comp"``, ``"cp_mol_liq_comp"``, ``"enth_mol_liq_comp"``, and ``"entr_mol_liq_comp"``.
In this example, we used the ``Perrys`` method for ``"H2O"`` and the ``Constant`` method for
both of our ions. When we declare a specific method to calculate these properties, we are then
REQUIRED to include specific parameter information in the ``"parameter_data"`` dictionary
defined within each component dictionary. For additional information regarding those parameter
needs, have a look at `Perrys`_ and `Constant`_ methods in IDAES.

**(3)** In this example, we are just setting up a configuration for water only. Thus, we are
not particularly interested in any other phases. In this case, we define the ``"phases"``
dictionary to contain a single phase we named ``'Liq'`` and declared this to be an ``AqueousPhase``.
In WaterTAP, most of our models will be using ``AqueousPhase``, but may add additional phases
for effects such as precipitation and/or gas-absorbtion. Also, it should be noted that each phase
must also define a method for the ``"equation_of_state"`` argument. In this case, we are assuming
that the phase behaves under the ``Ideal`` assumption. For more information on phases and equations
of state, see `Phases`_ and `EquationOfState`_.

**(4)** We chose to define the ``"state_definition"`` as ``FTPx``, however, there are many more
options available. More information can be found in `StateDefinition`_.

.. note::

    Much of the difficulties and complications with setting up a proper **thermo-properties**
    configuration dictionary can be handled by the **Electrolyte Database** system in
    WaterTAP (Documentation pending)



The **reaction-properties** configuration dictionary
----------------------------------------------------

If you did not include reactions in the **thermo-properties** dictionary
(see :ref:`How to use inherent reactions<how_to_use_inherent_reactions>`)
and your system involves reactions, then you MUST also create and
provide a **reaction-properties** configuration dictionary. Unlike the **thermo-properties**
configuration dictionary, most of the keys within the **reaction-properties** dictionary
are optional and depend on your system. The major keys to be aware of are as follows...

+-----------------------+-------------------------------------------------------------------------------------------+
|     Key               |  Description                                                                              |
+=======================+===========================================================================================+
| base_units            | dictionary containing the base units the model uses (same as the **thermo-properties**)   |
+-----------------------+-------------------------------------------------------------------------------------------+
| equilibrium_reactions | dictionary containing the full set of equilibrium reactions in the system                 |
+-----------------------+-------------------------------------------------------------------------------------------+
| rate_reactions        | dictionary containing the full set of rate reactions in the system                        |
+-----------------------+-------------------------------------------------------------------------------------------+

.. note::

    Each type of reaction (``equilibrium_reactions`` and ``rate_reactions``) have
    their own sets of parameters and methods to be declared. More information on
    how to set up these arguments and the methods available can be found at
    `GenericReactions`_. You can go directly to either methods by following
    the following links (`EquilibriumReactions`_ and `RateReactions`_).


Example reaction-properties configuration dictionary
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Following from our previous example for the **thermo-properties** configuration
dictionary, here we will show how you setup a **reaction-properties** configuration
dictionary for the dissociation of water. Since water dissociation is a fast acid
reaction, we will model it as an equilibrium reaction.

.. testcode::

    # Importing the object for units from pyomo
    from pyomo.environ import units as pyunits

    # Import the object/function for heat of reaction
    from idaes.generic_models.properties.core.reactions.dh_rxn import constant_dh_rxn

    # Import built-in Gibb's Energy function
    from idaes.generic_models.properties.core.reactions.equilibrium_constant import van_t_hoff

    # Import safe log power law equation
    from idaes.generic_models.properties.core.reactions.equilibrium_forms import log_power_law_equil

    # Importing the enum for concentration unit basis used in the 'get_concentration_term' function
    from idaes.generic_models.properties.core.generic.generic_reaction import ConcentrationForm

    reaction_config = {
        "base_units": {"time": pyunits.s,
                       "length": pyunits.m,
                       "mass": pyunits.kg,
                       "amount": pyunits.mol,
                       "temperature": pyunits.K},
        "equilibrium_reactions": {
            "H2O_Kw": {
                    "stoichiometry": {("Liq", "H2O"): -1,
                                     ("Liq", "H_+"): 1,
                                     ("Liq", "OH_-"): 1},
                   "heat_of_reaction": constant_dh_rxn,
                   "equilibrium_constant": van_t_hoff,
                   "equilibrium_form": log_power_law_equil,
                   "concentration_form": ConcentrationForm.moleFraction,
                   "parameter_data": {
                       "dh_rxn_ref": (55.830, pyunits.kJ/pyunits.mol),
                       "k_eq_ref": (10**-14/55.2/55.2, pyunits.dimensionless),
                       "T_eq_ref": (298, pyunits.K),

                       # By default, reaction orders follow stoichiometry, so
                       #    we manually set reaction order here to override.
                       #    In our case, the water dissociation reaction is
                       #    mathematically represented by Kw = [H_+]*[OH_-]
                       #    thus, this reaction is of order 0 with respect

                       #    to the [H2O] concentration.
                       "reaction_order": {("Liq", "H2O"): 0,
                                        ("Liq", "H_+"): 1,
                                        ("Liq", "OH_-"): 1}
                        }
                        # End parameter_data
                   }
                   # End reaction H2O_Kw
             }
             # End equilibrium_reactions
        }
        # End reaction_config definition

There is a significant amount of information and options available, so we will
just go through some things of note here.

**(1)** Each reaction you add to the model will have its own dictionary with
essentially the same format as the ``"H2O_Kw"`` dictionary shown above. Make sure
that each reaction you add has a unique key within the ``"equilibrium_reactions"``
parent dictionary.

**(2)** The first thing you need to define about the reaction is the stoichiometry.
In IDAES, we follow the convention that **products** of a reaction should have
positive stoichiometric values and **reactants** of a reaction should have negative
stoichiometric values. This is true for both ``equilibrium_reactions`` and
``rate_reactions``.

**(3)** The ``"stoichiometry"`` dictionary under the reaction has tuple keys. In
this format, the first item in the tuple is the ``phase`` of the species involved in
the reaction and the second item in the tuple is the ``name`` of the species. Recall
that in the **thermo-properties** configuration dictionary, we named the ``AqueousPhase``
as ``"Liq"``, thus we must reference that same name here in the **reaction-properties**
configuration dictionary. The specific species must also be referenced by the names
they were given in the **thermo-properties** configuration dictionary.

**(4)** You must provide methods/options for each of the following: ``"heat_of_reaction"``,
``"equilibrium_constant"``, ``"equilibrium_form"``, and ``"concentration_form"``. These
methods will define how IDAES computes the heat of reaction in the energy balance, the
equilibrium constant or K value for this reaction constraint, the mathematical representation
of the equilibrium constraint, and what the concentration form is for the species involved
in this reaction, respectively. Many options are available for all of these and more
information on each can be found at `ReactionMethods`_ and `ConcentrationForm`_.

**(5)** The ``"parameter_data"`` dictionary must contain the parameter information
required by the chosen methods from **(4)** above. See `ReactionMethods`_ for more
details.

**(6)** Within the ``"parameter_data"`` dictionary is an optional dictionary for
``"reaction_order"``. If this dictionary is not provided, then it is assumed that
the order of the reaction form with respect to each species just follows the
``"stoichiometry"`` dictionary from above. However, in certain cases you may need
to override that assumption. In this particular case, we override the reaction
order to zero out the order with respect to the water concentration. This is
standard practice for aqueous acid-base chemistry.

.. note::

    The ``"reaction_order"`` dictionary follows the same sign convention for products
    and reactants as the ``"stoichiometry"`` dictionary. Positive signs for products
    and negative signs for reactants.


Defining a **unit model**
-------------------------

Once you have your **thermo-properties** and (optionally) your **reaction-properties**
configuration dictionaries setup, you will want to put them into a **unit model** so
that you can simulate that particular unit process with the chemistry you have
specified. Within IDAES, their are numerous **unit models** to chose from that
will support the inclusion of these chemistry configurations. A list of the
**unit models** available, and how to use them, are provided here (`UnitModels`_).

In this guide, we will not cover all the **unit models**, but will give one basic
example of how to use the configuration dictionaries defined above with the
`EquilibriumReactor`_ model.


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
    from idaes.generic_models.properties.core.generic.generic_property import GenericParameterBlock
    from idaes.generic_models.properties.core.generic.generic_reaction import GenericReactionParameterBlock

    # Import the idaes object for the EquilibriumReactor unit model
    from idaes.generic_models.unit_models.equilibrium_reactor import EquilibriumReactor

    # Create an instance of a pyomo model
    model = ConcreteModel()

    # Add an IDAES flowsheet to that model
    model.fs = FlowsheetBlock(default={"dynamic": False})

    # Add a thermo parameter block to that flowsheet
    #   Here, we are passing our 'thermo_config' dictionary we created earlier
    model.fs.thermo_params = GenericParameterBlock(default=thermo_config)

    # Add a reaction parameter block to that flowsheet
    #   Here, we are passing our thermo block created above as the property package
    #   and then giving our 'reaction_config' as the instructions for how the
    #   reactions will be constructed from the thermo package.
    model.fs.rxn_params = GenericReactionParameterBlock(
                default={"property_package": model.fs.thermo_params, **reaction_config})

    # Add an EquilibriumReactor object as the unit model
    #   Here, we pass both the thermo package and reaction package, as well
    #   as a number of other arguments to help define how this unit process
    #   will behave.
    model.fs.unit = EquilibriumReactor(default={
                "property_package": model.fs.thermo_params,
                "reaction_package": model.fs.rxn_params,
                "has_rate_reactions": False,
                "has_equilibrium_reactions": True,
                "has_heat_transfer": False,
                "has_heat_of_reaction": False,
                "has_pressure_change": False})

    # At this point, you can 'fix' your inlet/outlet state conditions,
    #     setup scaling factors, initialize the model, then solve the model
    #     just as you would with any other IDAES flowsheet


In the example code above, we show how to setup the thermo and reaction packages
and place them into the `EquilibriumReactor` unit model, but do not go further.
Additional instructions for setting up and solving unit models can be found at
`IDAESWorkflow`_.
