Pump
====

.. index::
   pair: watertap.unit_models.pump;pump

This pump unit model is a class of pressure changer which is based on the `IDAES pressure changer <https://idaes-pse.readthedocs.io/en/latest/reference_guides/model_libraries/generic/unit_models/pressure_changer.html>`_
and makes the following assumptions:

   * supports a single liquid phase only
   * supports steady-state only
   * supports isothermal pump only

Degrees of Freedom
------------------

Pump units generally have two or more degrees of freedom, depending on the pump efficiency model used.

Typical fixed variables are:

* outlet pressure, :math:`P_{out}` or :math:`\Delta P`
* pump efficiency

Model Structure
---------------

The pump unit model consists of a single ``ControlVolume0D`` (named ``control_volume``)
with one Inlet Port (named ``inlet``) and one Outlet Port (named ``outlet``).


Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq']"
   "Components", ":math:`j`", "['H2O', 'NaCl']*"

\*Solute depends on the imported property model; example shown here is for the NaCl property model.

Variables
---------

When the configuration ``variable_efficiency`` is set to default ``VariableEfficiency.none``, the variables are the same listed in Variable tab of
`IDAES pressure changer <https://idaes-pse.readthedocs.io/en/latest/reference_guides/model_libraries/generic/unit_models/pressure_changer.html>`_.

When the configuration ``variable_efficiency`` is not ``VariableEfficiency.none``, then there are three additional variables:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Best efficiency point flowrate of the centrifugal pump", ":math:`F_{bp}`", "``bep_flow``", "None", ":math:`\text{m}^3/\text{s}`"
   "Best efficiency of the centrifugal pump", ":math:`\eta_{bep}`", "``bep_eta``", "None", "dimensionless"
   "Ratio of pump flowrate to best efficiency point flowrate", ":math:`r_{bep, t}`", "``flow_ratio``", "[t]", "dimensionless"

Equations and Relationships
---------------------------
When the configuration ``variable_efficiency`` is set to default ``VariableEfficiency.none``, the constraints are the same listed in
`IDAES pressure changer <https://idaes-pse.readthedocs.io/en/latest/reference_guides/model_libraries/generic/unit_models/pressure_changer.html>`_.

When the configuration ``variable_efficiency`` is not ``VariableEfficiency.none``, there are two additional constraints:

.. csv-table::
   :header: "Description", "Equation"

   "Pump flow ratio", ":math:`F_{in, t} = r_{bp, t} * F_{bp}`"
   "Actual pump efficiency", ":math:`\eta_{pump, t} = \eta_{bp} * r_{bp, t}`"

When the configuration ``variable_efficiency`` is set to ``VariableEfficiency.flow``,
then the pump efficiency only depends on flow:

    .. math::

       \eta_{pump, t}=
       \begin{cases}
         0.4 & \text{for } r_{bp, t} < 0.6\\
         -0.995 * r_{bp, t}^{2} + 1.977 * r_{bp, t} + 0.018 & \text{for } 0.6 \le r_{bp, t} \le 1.4\\
         0.4 & \text{for } r_{bp, t} > 1.4
       \end{cases}

Class Documentation
-------------------
.. currentmodule:: watertap.unit_models.pressure_changer

.. autoclass:: Pump
    :members:
    :noindex: