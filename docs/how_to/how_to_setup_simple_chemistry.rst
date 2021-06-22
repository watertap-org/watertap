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
