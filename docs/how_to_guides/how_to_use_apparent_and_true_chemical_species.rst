How to use apparent and true chemical species
=============================================

.. note:: This page provides a manual approach to building an IDAES configuration dictionary.
    The same result can be achieved by using the :doc:`Electrolyte Database (EDB)</technical_reference/edb/index>`.

In WaterTAP, most all chemical processes simulated will be considered "true"
species, i.e., species that actually exist in an aqueous solution (e.g., Na+ and Cl- for NaCl). However, there
may be times when you want to have your system report the "apparent" species
between your inlet and/or outlet ports of a flowsheet. To do this, you need
to understand how to define a species in a given configuration file as "apparent",
which is a special component type in the **thermo-properties** configuration dictionary.

For more information on creating a **thermo-properties** configuration dictionary,
see :ref:`How to setup simple chemistry<how_to_setup_simple_chemistry>`.


What you need to update in the thermo-properties configuration dictionary
-------------------------------------------------------------------------

1. Add the ``"Apparent"`` and ``"StateIndex"`` objects under your import statements
2. Specify certain species as ``"Apparent"`` within its component definition
3. Give a list of "true" species it will dissociate into
4. Define a "state_components" argument in the config and specify ``"StateIndex.true"``


Example Configuration Dictionary
--------------------------------

This is best understood by going through an example. In this case, we will consider
building off of the :ref:`How to setup simple chemistry<how_to_setup_simple_chemistry>` guide and add ``"NaCl"`` as
an ``"Apparent"`` species.

.. testsetup::

   # quiet idaes logs
   import idaes.logger as idaeslogger
   idaeslogger.getLogger('ideas.core').setLevel('CRITICAL')
   idaeslogger.getLogger('idaes.init').setLevel('CRITICAL')

.. testcode::

  from pyomo.environ import units as pyunits
  from idaes.core import AqueousPhase

  # Add Apparent in this import statement
  from idaes.core.components import Solvent, Cation, Anion, Apparent

  import idaes.generic_models.properties.core.pure.Perrys as Perrys
  from idaes.generic_models.properties.core.pure.ConstantProperties import Constant
  from idaes.generic_models.properties.core.state_definitions import FTPx
  from idaes.generic_models.properties.core.eos.ideal import Ideal

  # Add the import for StateIndex
  from idaes.generic_models.properties.core.generic.generic_property import StateIndex

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
                      },
            'Na_+': {"type": Cation, "charge": 1,
                  # Define the methods used to calculate the following properties
                  "dens_mol_liq_comp": Constant,
                  "enth_mol_liq_comp": Constant,
                  "cp_mol_liq_comp": Constant,
                  "entr_mol_liq_comp": Constant,
                  # Parameter data is always associated with the methods defined above
                  "parameter_data": {
                        "mw": (22.989769, pyunits.g/pyunits.mol),
                        "dens_mol_liq_comp_coeff": (55, pyunits.kmol*pyunits.m**-3),
                        "enth_mol_form_liq_comp_ref": (-240.1, pyunits.kJ/pyunits.mol),
                        "cp_mol_liq_comp_coeff": (75000, pyunits.J/pyunits.kmol/pyunits.K),
                        "entr_mol_form_liq_comp_ref": (59, pyunits.J/pyunits.K/pyunits.mol)
                                    },
                        # End parameter_data
                        },
            'Cl_-': {"type": Anion, "charge": -1,
                  # Define the methods used to calculate the following properties
                  "dens_mol_liq_comp": Constant,
                  "enth_mol_liq_comp": Constant,
                  "cp_mol_liq_comp": Constant,
                  "entr_mol_liq_comp": Constant,
                  # Parameter data is always associated with the methods defined above
                  "parameter_data": {
                        "mw": (35.453, pyunits.g/pyunits.mol),
                        "dens_mol_liq_comp_coeff": (55, pyunits.kmol*pyunits.m**-3),
                        "enth_mol_form_liq_comp_ref": (-167.2, pyunits.kJ/pyunits.mol),
                        "cp_mol_liq_comp_coeff": (75000, pyunits.J/pyunits.kmol/pyunits.K),
                        "entr_mol_form_liq_comp_ref": (56.5, pyunits.J/pyunits.K/pyunits.mol)
                                    },
                        # End parameter_data
                        },

            # This is how an Apparent species is defined in the configuration dictionary
            #   it requires the same parameter arguments as True species, but also needs
            #   a dictionary for "dissociation_species" that tells how much of each
            #   true species this Apparent species is formed from.
            'NaCl': {"type": Apparent,
                  "dissociation_species": {"Na_+":1, "Cl_-":1},
                  # Define the methods used to calculate the following properties
                  "dens_mol_liq_comp": Constant,
                  "enth_mol_liq_comp": Constant,
                  "cp_mol_liq_comp": Constant,
                  "entr_mol_liq_comp": Constant,
                  # Parameter data is always associated with the methods defined above
                  "parameter_data": {
                      "dens_mol_liq_comp_coeff": (55, pyunits.kmol*pyunits.m**-3),
                      "enth_mol_form_liq_comp_ref": (-945.53, pyunits.kJ/pyunits.mol),
                      "cp_mol_liq_comp_coeff": (167039, pyunits.J/pyunits.kmol/pyunits.K),
                      "entr_mol_form_liq_comp_ref": (100, pyunits.J/pyunits.K/pyunits.mol)
                                    },
                            # End parameter_data
                            },
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

          # We must define the 'StateIndex' as "true". This is because in WaterTAP,
          #   all speciation reactions are defined on the true species, not the
          #   apparent species.
          "state_components": StateIndex.true,

          "pressure_ref": 1e5,
          "temperature_ref": 300,
          "base_units": {"time": pyunits.s,
                         "length": pyunits.m,
                         "mass": pyunits.kg,
                         "amount": pyunits.mol,
                         "temperature": pyunits.K},
      }
      # End thermo_config definition


.. note::

    When you define a species as ``"Apparent"`` and specify ``"state_components": StateIndex.true``,
    you cannot reference that species as part of your inlet variables or in any
    reactions in the system. The ``"StateIndex"`` is used to define what species
    can be used in reactions or in the inlet ports to set initial states. For WaterTAP,
    we will always define reactions on a true species basis.
