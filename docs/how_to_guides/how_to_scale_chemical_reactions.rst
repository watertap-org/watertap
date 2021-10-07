.. _how_to_scale_chemical_reactions:

How to scale chemical reactions
===============================

.. warning::
    These scaling methods are meant to be used with the 'gradient-based' scaling
    functionality built into ipopt. It is not advised to use 'user-scaling' with
    the chemical modules, as this tends to have poor convergence. 

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

    The following guide and examples assume that your **GenericParameterBlock**
    is named **thermo_params**, your **GenericReactionParameterBlock** is named
    **rxn_params**, your unit model is named **unit**, and your **FlowsheetBlock**
    is named **fs**. For information on setup of the basic chemistry modules,
    see :ref:`How to setup simple chemistry<how_to_setup_simple_chemistry>`.


.. note::

    The scaling factor constants used in these demonstrations are the best found
    for these types of problems based on numerical testing of a variety of problems.
    These values can and may need to change for your particular problem.


.. note::

    IMPORTANT: Scaling of just the chemical reactions is insufficient for solving
    a chemistry module. User's MUST also scale the chemical species in the system
    due to the dilute nature of aqueous chemistry. See
    :ref:`How to scale chemical species<how_to_scale_chemical_species>`.


Types of Reactions
------------------

1. Inherent Reactions using log_power_law_equil form
2. Equilibrium Reactions using log_power_law_equil form
3. Rate Reactions using power_law_rate form
4. Stoichiometric Reactions (Documentation Pending)


Inherent Reactions using log_power_law_equil form
-------------------------------------------------

Inherent reactions (see :ref:`How to use inherent reactions<how_to_use_inherent_reactions>`)
are placed into the **GenericParameterBlock**
rather than the **GenericReactionParameterBlock**. As such, the scaling methods employed
must be applied to the **thermo_params** for this set of reactions. The code segment below
shows a demonstration of scaling factors that generally work well for these types of reactions.


Inherent Scaling Demonstration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block::

    # Iterate through the reactions to set appropriate eps values
    factor = 1e-4
    for rid in model.fs.thermo_params.inherent_reaction_idx:
        scale = value(model.fs.unit.control_volume.properties_out[0.0].k_eq[rid].expr)
        # Want to set eps in some fashion similar to this
        if scale < 1e-16:
            model.fs.thermo_params.component("reaction_"+rid).eps.value = scale*factor
        else:
            model.fs.thermo_params.component("reaction_"+rid).eps.value = 1e-16*factor

    for i in model.fs.unit.control_volume.inherent_reaction_extent_index:
        scale = value(model.fs.unit.control_volume.properties_out[0.0].k_eq[i[1]].expr)
        iscale.set_scaling_factor(model.fs.unit.control_volume.inherent_reaction_extent[0.0,i[1]], 10/scale)
        iscale.constraint_scaling_transform(model.fs.unit.control_volume.properties_out[0.0].
                inherent_equilibrium_constraint[i[1]], 0.1)


What this code segment is doing is to first assign an **eps** "smoothing" value for the **log_power_law_equil**
function. While this value is not directly impactful to scaling, it is necessary to provide a reasonable
estimate for the smoothing parameter, since a poor value will negatively impact convergence. The default value
is generally insufficient, so we change that value to always be at least 4 orders of magnitude less than
the **k_eq** value calculated for each reaction.

Next, we create scaling factors for the **inherent_reaction_extent** variable based on the inverse
of that calculated **k_eq** value.

Lastly, we perform a scaling transformation of the **inherent_equilibrium_constraint**. In our case,
we scale this constraint down by a factor of 0.1 because the **log_power_law_equil** constraint
form is already scaled up reasonably well.


Equilibrium Reactions using log_power_law_equil form
----------------------------------------------------

Equilibrium reactions (see :ref:`How to setup simple chemistry<how_to_setup_simple_chemistry>`)
are placed into the **GenericReactionParameterBlock**. As such, the scaling methods employed
must be applied to the **rxn_params** for this set of reactions. The code segment below
shows a demonstration of scaling factors that generally work well for these types of reactions.


Equilibrium Scaling Demonstration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block::

    # Equilibrium reactions have eps in the 'rxn_params'
    factor = 1e-4
    for rid in model.fs.rxn_params.equilibrium_reaction_idx:
        scale = value(model.fs.unit.control_volume.reactions[0.0].k_eq[rid].expr)
        # Want to set eps in some fashion similar to this
        if scale < 1e-16:
            model.fs.rxn_params.component("reaction_"+rid).eps.value = scale*factor
        else:
            model.fs.rxn_params.component("reaction_"+rid).eps.value = 1e-16*factor

    for i in model.fs.unit.control_volume.equilibrium_reaction_extent_index:
        scale = value(model.fs.unit.control_volume.reactions[0.0].k_eq[i[1]].expr)
        iscale.set_scaling_factor(model.fs.unit.control_volume.equilibrium_reaction_extent[0.0,i[1]], 10/scale)
        iscale.constraint_scaling_transform(model.fs.unit.control_volume.reactions[0.0].
                equilibrium_constraint[i[1]], 0.1)


.. note::

    These scaling arguments are identical to the **Inherent Reaction** scaling methods,
    however, because these reactions exist in a different location of the model, we
    showed this here for completeness. All reactions, regardless of location, need scaling.


Rate Reactions using power_law_rate form
----------------------------------------

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
        iscale.set_scaling_factor(model.fs.unit.control_volume.rate_reaction_extent[0.0,i], 10/scale)


Stoichiometric Reactions
------------------------

.. note::

    Documentation under development.
