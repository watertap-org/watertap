.. _how_to_scale_chemical_species:

How to scale chemical species
=============================

In :ref:`How to scale chemical reactions<how_to_scale_chemical_reactions>`, we
discussed at length how you need to apply scaling factors to variables and constraints
associated with all the different types of reactions your system may have. Here,
we will add to that discussion how the variables for the chemical species in the
system should also be scaled.

.. note::

    IMPORTANT: Scaling of just the chemical species is insufficient for solving
    a chemistry module. User's MUST also scale the chemical reactions in the system
    due to the dilute nature of aqueous chemistry. See
    :ref:`How to scale chemical reactions<how_to_scale_chemical_reactions>`. Additionally,
    user's MUST also scale the energy balance equations. See also
    :ref:`How to scale energy balance for chemistry<how_to_scale_chemical_process_energy_balance>`.

After you have set all scaling factors for reactions, species, energy balances, etc.,
you need to call the 'calculate_scaling_factors' function on the entire model before
attempting to solve. See below:

.. code-block::

    # After all scaling factors and constraint transformations are complete,
    # call the following function on the model.
    iscale.calculate_scaling_factors(model.fs.unit)


Scaling of the chemical species depends on what your declared 'state_definition' was
in your **thermo_config**. See :ref:`How to setup simple chemistry<how_to_setup_simple_chemistry>`
for the detailed discussion on setting up a **thermo_config**. In this document, we will
discuss scaling chemical species for:

1. FTPx State Variables
2. FpcTP State Variables


Scaling for FTPx State Variables
--------------------------------

.. note::

    The following example assumes that your unit model is named **unit** and
    your **FlowsheetBlock** is named **fs**. For information on setup of the
    basic chemistry modules,
    see :ref:`How to setup simple chemistry<how_to_setup_simple_chemistry>`.

For this example, our **GenericParameterBlock** (named **thermo_params**) is using
the "state_definition" of **FTPx**.

.. code-block::

    # Specify a minimum division factor for scaling
    min = 1e-3
    # For species
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
            model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[i],
            10 / scale,
        )
        iscale.set_scaling_factor(
            model.fs.unit.control_volume.properties_out[0.0].flow_mol_phase_comp[i],
            10 / scale,
        )
        iscale.constraint_scaling_transform(
            model.fs.unit.control_volume.properties_out[0.0].component_flow_balances[i[1]],
            10 / scale,
        )
        iscale.constraint_scaling_transform(
            model.fs.unit.control_volume.material_balances[0.0, i[1]], 10 / scale
        )

    if hasattr(model.fs.unit.control_volume, "volume"):
        iscale.set_scaling_factor(
            model.fs.unit.control_volume.volume, 10 / model.fs.unit.volume[0.0].value
        )



Scaling for FpcTP State Variables
---------------------------------

.. note::

    The following example assumes that your unit model is named **unit** and
    your **FlowsheetBlock** is named **fs**. For information on setup of the
    basic chemistry modules,
    see :ref:`How to setup simple chemistry<how_to_setup_simple_chemistry>`.

For this example, our **GenericParameterBlock** (named **thermo_params**) is using
the "state_definition" of **FpcTP**.

.. code-block::
    
    # Specify a minimum division factor for scaling
    min = 1e-3
    # For species
    for i in model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp:
        # i[0] = phase, i[1] = species
        if model.fs.unit.inlet.flow_mol_phase_comp[0, i[0], i[1]].value > min:
            scale = model.fs.unit.inlet.flow_mol_phase_comp[0, i[0], i[1]].value
        else:
            scale = min

        iscale.set_scaling_factor(
            model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp[i[1]],
            10 / scale,
        )
        iscale.set_scaling_factor(
            model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[i],
            10 / scale,
        )
        iscale.set_scaling_factor(
            model.fs.unit.control_volume.properties_out[0.0].flow_mol_phase_comp[i],
            10 / scale,
        )
        iscale.constraint_scaling_transform(
            model.fs.unit.control_volume.material_balances[0.0, i[1]], 10 / scale
        )

    if hasattr(model.fs.unit.control_volume, "volume"):
        iscale.set_scaling_factor(
            model.fs.unit.control_volume.volume, 10 / model.fs.unit.volume[0.0].value
        )
