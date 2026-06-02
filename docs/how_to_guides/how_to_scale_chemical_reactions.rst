.. _how_to_scale_chemical_reactions:

How to scale chemical reactions
===============================

.. note::

    The following guide and examples assume that your **GenericParameterBlock**
    is named **thermo_params**, your **GenericReactionParameterBlock** is named
    **rxn_params**, your unit model is named **unit**, and your **FlowsheetBlock**
    is named **fs**. For information on setup of the basic chemistry modules,
    see :ref:`How to setup simple chemistry<how_to_setup_simple_chemistry>`.

No model in WaterTAP can solve without proper scaling of the constraints
and variables within that model. This can be a difficult task for aqueous
chemistry systems as these systems must deal with extremely dilute or trace
chemical species and reaction coefficients that may vary from each other
by several orders of magnitude.

To scale a set of chemical reactions properly requires different function
arguments to be invoked depending on: (i) where the reactions are located
in the model and (ii) what types of reactions are called for. Each of these
different scenarios will be discussed and demonstrated below.

.. note::

    IMPORTANT: Scaling of just the chemical reactions is insufficient for solving
    a chemistry module. User's MUST also scale the chemical species in the system
    due to the dilute nature of aqueous chemistry. See
    :ref:`How to scale chemical species<how_to_scale_chemical_species>`. Additionally,
    user's MUST also scale the energy balance equations. See also
    :ref:`How to scale energy balance for chemistry<how_to_scale_chemical_process_energy_balance>`.

After you have set all scaling factors for reactions, species, energy balances, etc.,
you need to call the 'calculate_scaling_factors' function on the entire model before
attempting to solve. See below:

.. code-block::

    # After all scaling factors and constraint transformations are complete,
    # call the following function on the model.
    iscale.calculate_scaling_factors(model.fs.unit)


Types of Reactions
------------------

1. Inherent Reactions using **log_power_law_equil** form
2. Equilibrium Reactions using **log_power_law_equil** form
3. Rate Reactions using **power_law_rate** form
4. Solubility/Precipitation Reaction using **log_solubility_product** form
5. Stoichiometric Reactions



Inherent Reactions using **log_power_law_equil** form
-----------------------------------------------------

Inherent reactions (see :ref:`How to use inherent reactions<how_to_use_inherent_reactions>`)
are placed into the **GenericParameterBlock**
rather than the **GenericReactionParameterBlock**. As such, the scaling methods employed
must be based on the **thermo_config** dictionary for this set of reactions. The code segment below
shows a demonstration of scaling factors that generally work well for these types of reactions.


Inherent Scaling Demonstration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block::

    # Specify a minimum allowable scaling division factor
    min_scale = 1e-3
    # Add scaling factors for reactions by looping through the reactions
    for i in model.fs.unit.control_volume.inherent_reaction_extent_index:
        # i[0] = time, i[1] = reaction

        # Grab the 'k_eq_ref' value directly from the thermo config dictionary
        scale = max(
            min_scale,
            thermo_config["inherent_reactions"][i[1]]["parameter_data"]["k_eq_ref"][0],
        )
        iscale.set_scaling_factor(
            model.fs.unit.control_volume.inherent_reaction_extent[0.0, i[1]], 10 / scale
        )
        iscale.constraint_scaling_transform(
            model.fs.unit.control_volume.properties_out[
                0.0
            ].inherent_equilibrium_constraint[i[1]],
            0.1,
        )


What this code segment is doing is to first create scaling factors for the
**inherent_reaction_extent** variable based on the inverse
of the **k_eq_ref** value (or based on some minimum inverse value).

Lastly, we perform a scaling transformation of the **inherent_equilibrium_constraint**. In our case,
we scale this constraint down by a factor of 0.1 because the **log_power_law_equil** constraint
form is already scaled up reasonably well.


Equilibrium Reactions using **log_power_law_equil** form
--------------------------------------------------------

Equilibrium reactions (see :ref:`How to setup simple chemistry<how_to_setup_simple_chemistry>`)
are placed into the **GenericReactionParameterBlock**. As such, the scaling methods employed
must be applied to the **rxn_config** dictionary for this set of reactions. The code segment below
shows a demonstration of scaling factors that generally work well for these types of reactions.


Equilibrium Scaling Demonstration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block::

    # Specify a minimum allowable scaling division factor
    min_scale = 1e-3
    # Add scaling factors for reactions
    for i in model.fs.unit.control_volume.equilibrium_reaction_extent_index:
        # i[0] = time, i[1] = reaction

        # Grab the 'k_eq_ref' value from the reaction config
        scale = max(
            min_scale,
            rxn_config["equilibrium_reactions"][i[1]]["parameter_data"]["k_eq_ref"][0],
        )
        iscale.set_scaling_factor(
            model.fs.unit.control_volume.equilibrium_reaction_extent[0.0, i[1]], 10 / scale
        )
        iscale.constraint_scaling_transform(
            model.fs.unit.control_volume.reactions[0.0].equilibrium_constraint[i[1]], 0.1
        )


.. note::

    These scaling arguments are identical to the **Inherent Reaction** scaling methods,
    however, because these reactions exist in a different location of the model, we
    showed this here for completeness. All reactions, regardless of location, need scaling.


Rate Reactions using **power_law_rate** form
--------------------------------------------

Rate reactions only exist in the **GenericReactionParameterBlock** and so these scaling
arguments apply to **rxn_params** for these types of reactions. These are much simpler to
scale than both the **Inherent** and **Equilibrium** reactions, but are just as important
to apply scaling for. Below is a demonstration of applying scaling.

Rate Reaction Scaling Demonstration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block::

    # Scaling for kinetic reactions
    for i in model.fs.rxn_params.rate_reaction_idx:
        scale = value(model.fs.unit.control_volume.reactions[0.0].reaction_rate[i].expr)
        iscale.set_scaling_factor(
            model.fs.unit.control_volume.rate_reaction_extent[0.0, i], 1000 / scale
        )



.. note::

    We only need to call 'set_scaling_factor' here and NOT call 'constraint_scaling_transform'
    because this scaling factor will cascade into the constraints automatically once you call
    'calculate_scaling_factors' on the model. This is different from other reactions because
    there is no 'log form' for rate reactions. The 'log form' always requires some additional
    treatment.



Solubility/Precipitation Reaction using **log_solubility_product** form
-----------------------------------------------------------------------

To scale these reactions, you will use the same methods outlined above for **Equilibrium**
and **Inherent** reactions. However, there is an additional step. That additional step involves
setting a smoothing parameter **eps** (which is a factor unique to the 'log_solubility_product'
function). Below is a demonstration of setting up that smoothing parameter assuming your
solubility reactions are in the **rxn_params** object and the **rxn_config** dictionary.

Setting **eps** Smoothing Factor for Solubility Products
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block::

    # Specify a minimum allowable scaling factor for the eps
    factor = 1e-2
    for rid in model.fs.rxn_params.equilibrium_reaction_idx:
        # Grab the 'k_eq_ref' value from the reaction config
        scale = rxn_config["equilibrium_reactions"][rid]["parameter_data"]["k_eq_ref"][0]

        # NOTE: ONLY certain functions have an eps value that we need to set
        if hasattr(rxn_params.component("reaction_" + rid), "eps"):
            # highest allowable value for setting eps based on k_eq_ref
            if scale < 1e-16:
                model.fs.rxn_params.component("reaction_" + rid).eps.value = scale * factor
            else:
                model.fs.rxn_params.component("reaction_" + rid).eps.value = 1e-16 * factor



Stoichiometric Reactions
------------------------

Stoichiometric reactions are generally the simplest to scale. However, determining
how much to scale them by is not always clear. It depends on what are the expected
changes in molar flows due to the reaction. Since there is not always a clear way
to determine this, the demonstration below simply shows you where the scaling is
applied to within the framework.


Setting Scaling Factor for Stoichiometric Reaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The sample below just gives a brief demo of how to scale a stoichiometric reaction
named 'R1' by a given factor. All stoichiometric reactions are located in the
'control_volume' of the model and named 'rate_reaction_extent'. This is because
**Rate Reactions** and **Stoichiometric Reactions** have very similar implementations
in the IDAES framework.

.. code-block::
    
    # Specify a factor to scale by
    factor = 1
    iscale.set_scaling_factor(
        m.fs.unit.control_volume.rate_reaction_extent[0.0, "R1"], factor
    )

