How to setup simple chemistry
=============================

.. _GenericProperties: https://idaes-pse.readthedocs.io/en/stable/user_guide/components/property_package/general/index.html#generic-property-package-framework
.. _GenericReactions: https://idaes-pse.readthedocs.io/en/stable/user_guide/components/property_package/general_reactions/index.html
.. _Perrys: https://idaes-pse.readthedocs.io/en/stable/user_guide/components/property_package/general/pure/Perrys.html
.. _Constant: https://idaes-pse.readthedocs.io/en/stable/user_guide/components/property_package/general/pure/ConstantProperties.html
.. _StateDefinition: https://idaes-pse.readthedocs.io/en/stable/user_guide/components/property_package/general/state_definition.html
.. _EquationOfState: https://idaes-pse.readthedocs.io/en/stable/user_guide/components/property_package/general/eos/ideal.html
.. _Components: https://idaes-pse.readthedocs.io/en/stable/user_guide/components/property_package/general/component_def.html
.. _Phases: https://idaes-pse.readthedocs.io/en/stable/user_guide/components/property_package/general/phase_def.html

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
                        },
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
Generally, the ``"H2O"`` component should always be a ``Solvent`` within ProteusLib. Charged species
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
In ProteusLib, most of our models will be using ``AqueousPhase``, but may add additional phases
for effects such as precipitation and/or gas-absorbtion. Also, it should be noted that each phase
must also define a method for the ``"equation_of_state"`` argument. In this case, we are assuming
that the phase behaves under the ``Ideal`` assumption. For more information on phases and equations
of state, see `Phases`_ and `EquationOfState`_.

**(4)** We chose to define the ``"state_definition"`` as ``FTPx``, however, there are many more
options available. More information can be found in `StateDefinition`_.


The **reaction-properties** configuration dictionary
----------------------------------------------------

If you did not include reactions in the **thermo-properties** dictionary (a topic
discussed later) and your system involves reactions, then you MUST also create and
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
