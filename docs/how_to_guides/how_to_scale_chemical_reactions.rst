How to scale chemical reactions
===============================

No model in ProteusLib can solve without proper scaling of the constraints
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


Types of Reactions
------------------

1. Inherent Reactions using log_power_law_equil form
2. Equilibrium Reactions using log_power_law_equil form
3. Rate Reactions using power_law_rate form
4. Stoichiometric Reactions (Documentation Pending)


Inherent Reactions using log_power_law_equil form
-------------------------------------------------

Inherent reactions (see Documentation Pending) are placed into the **GenericParameterBlock**
rather than the **GenericReactionParameterBlock**. As such, the scaling methods employed
must be applied to the **thermo_params** for this set of reactions. The code segment below
shows a demonstration of scaling factors that generally work well for these types of reactions.


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
