.. _how_to_scale_chemical_species:

How to scale chemical species
=============================

.. warning::
    This scaling methods are meant to be used with the 'gradient-based' scaling
    functionality built into ipopt. It is not advised to use 'user-scaling' with
    the chemical modules, as this tends to have poor convergence.

In :ref:`How to scale chemical reactions<how_to_scale_chemical_reactions>`, we
discussed at length how you need to apply scaling factors to variables and constraints
associated with all the different types of reactions your system may have. Here,
we will add to that discussion how the variables for the chemical species in the
system should also be scaled.

.. note::

    The scaling factor constants used in this demonstration are the best found
    for these types of problems based on numerical testing of a variety of problems.
    These values can and may need to change for your particular problem.


Scaling Chemical Species Example
--------------------------------

.. note::

    The following example assumes that your unit model is named **unit** and
    your **FlowsheetBlock** is named **fs**. For information on setup of the
    basic chemistry modules,
    see :ref:`How to setup simple chemistry<how_to_setup_simple_chemistry>`.

For this example, our **GenericParameterBlock** (named **thermo_params**) is using
the "state_definition" of **FTPx**. This may change what variables the scaling
is applied to. See :ref:`How to setup simple chemistry<how_to_setup_simple_chemistry>`
for more information on "state_definition" options.

.. code-block::

    # Try adding scaling for species
    min = 1e-6
    for i in model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp:
        # i[0] = phase, i[1] = species
        if model.fs.unit.inlet.mole_frac_comp[0, i[1]].value > min:
            scale = model.fs.unit.inlet.mole_frac_comp[0, i[1]].value
        else:
            scale = min
        iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp[i[1]], 10/scale)
        iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].mole_frac_phase_comp[i], 10/scale)
        iscale.set_scaling_factor(model.fs.unit.control_volume.properties_out[0.0].flow_mol_phase_comp[i], 10/scale)
        iscale.constraint_scaling_transform(
            model.fs.unit.control_volume.properties_out[0.0].component_flow_balances[i[1]], 10/scale)
        iscale.constraint_scaling_transform(model.fs.unit.control_volume.material_balances[0.0,i[1]], 10/scale)
