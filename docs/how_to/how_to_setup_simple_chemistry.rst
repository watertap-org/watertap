How to setup simple chemistry
-----------------------------

In ProteusLib, chemistry modules leverage the Generic Properties
(https://idaes-pse.readthedocs.io/en/stable/user_guide/components/property_package/general/index.html#generic-property-package-framework)
and Generic Reactions
(https://idaes-pse.readthedocs.io/en/stable/user_guide/components/property_package/general_reactions/index.html)
objects in IDAES. These objects can be used in conjunction with any unit process
where chemical reactions need to be considered. In this guide, we will cover how
to use these built-in objects within your own unit process.

What you will need
^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
