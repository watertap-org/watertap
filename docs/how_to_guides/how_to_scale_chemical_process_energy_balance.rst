.. _how_to_scale_chemical_process_energy_balance:

How to scale energy balance for chemical process
================================================

In :ref:`How to scale chemical reactions<how_to_scale_chemical_reactions>` and
:ref:`How to scale chemical species<how_to_scale_chemical_species>` we covered
the most complex parts associated with scaling aqueous chemistry models in
WaterTAP. The last piece to mention is that we also must scale the energy balances
for the model.

.. note::

    IMPORTANT: Scaling of just the energy balances is insufficient for solving
    a chemistry module. User's MUST also scale the chemical reactions in the system
    due to the dilute nature of aqueous chemistry. See
    :ref:`How to scale chemical reactions<how_to_scale_chemical_reactions>`. Additionally,
    user's MUST also scale the chemical species. See also
    :ref:`How to scale chemical species<how_to_scale_chemical_species>`.

After you have set all scaling factors for reactions, species, energy balances, etc.,
you need to call the 'calculate_scaling_factors' function on the entire model before
attempting to solve. See below:

.. code-block::

    # After all scaling factors and constraint transformations are complete,
    # call the following function on the model.
    iscale.calculate_scaling_factors(model.fs.unit)


Scaling the energy balances
---------------------------

.. note::

    The following example assumes that your unit model is named **unit** and
    your **FlowsheetBlock** is named **fs**. For information on setup of the
    basic chemistry modules,
    see :ref:`How to setup simple chemistry<how_to_setup_simple_chemistry>`.

The scaling of the energy balances can be done based solely on how the enthalpy
flow term expressions are formulated in the model. This can be accomplished
simply by doing the following:

.. code-block::

    # Setup default min and max for scaling
    max = 1
    min = 1
    for phase in model.fs.unit.control_volume.properties_in[0.0].enth_mol_phase:
        val = abs(
            value(
                model.fs.unit.control_volume.properties_in[0.0].enth_mol_phase[phase].expr
            )
        )
        if val >= max:
            max = val
        if val <= min:
            val = min
        iscale.set_scaling_factor(
            model.fs.unit.control_volume.properties_in[0.0]._enthalpy_flow_term[phase],
            10 / val,
        )
        iscale.set_scaling_factor(
            model.fs.unit.control_volume.properties_out[0.0]._enthalpy_flow_term[phase],
            10 / val,
        )

    iscale.constraint_scaling_transform(
        model.fs.unit.control_volume.enthalpy_balances[0.0], 10 / max
    )

