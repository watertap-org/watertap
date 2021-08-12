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
