.. _MVC:

Mechanical Vapor Compression (MVC)
==================================

.. code-block:: python

   from watertap.unit_models.mvc.components import Evaporator, Compressor, Condenser

Modeling a mechanical vapor compression (MVC) system is performed through flowsheet connectivity of individual model
block components. The model is simulated under the following criteria:

   * property package(s) must support liquid and vapor phases
   * all model components are 0D
   * all model components support steady-state only
   * the operating conditions of the evaporator are fixed
   * compressor performance is governed by isentropic efficiency
   * assumed complete condensation in the condenser

.. index::
   pair: watertap.unit_models.mvc;mvc

.. currentmodule:: watertap.unit_models.mvc

Introduction
------------
The MVC model is simulated as a modular combination of evaporator, compressor, and complete condenser unit models.
These three models may be connected in a flowsheet for the fundamental MVC desalination process. For additional
integration, pump and heat exchange unit operations may be added for pre- and post-processing of feed, brine, and
distillate streams. Due to this connectivity, specialized initialization routines and modified unit models have been
developed for this system and are located in a dedicated MVC directory.

Degrees of Freedom
------------------
Degrees of freedom for the MVC system are separated into the independent unit models. In this documentation, the
description of the MVC system will only include the evaporator, compressor, and complete condenser unit models. The
total system has 6 degrees of freedom in addition to the inlet feed state variables. Beginning from the feed, the
evaporator has 3 degrees of freedom:

   * 1 outlet condition (brine or vapor temperature or pressure)
   * overall heat transfer coefficient
   * area for heat exchange

The compressor which the vapor enters has 3 degrees of freedom:

   * pressure ratio
   * work (control volume variable)
   * isentropic efficiency

By utilizing the evaporator in the MVC directory a custom ``connect_to_condenser`` routine equates the heat exchange
between the evaporator and condenser unit models. Assuming the vapor will completely condense, there are no degrees of
freedom in the condenser. When combined in a flowsheet, water recovery may also be calculated. This may be used in place
of the specification of the evaporator outlet condition.

Model Structure
------------------
Construction of the ports, state blocks, and use of control volumes are defined in the variables section and separated
by unit. Property packages must be declared for the liquid and vapor phases of evaporator and condenser models. This
may be a two-phase property package or two different property packages respective to each phase.

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq']"
   "Components", ":math:`j`", "['H2O', other]"

.. _MVC_variables:

Variables
----------

Evaporator Variables
^^^^^^^^^^^^^^^^^^^^

The evaporator contains 3 state blocks corresponding to the feed inlet, brine outlet, and vapor outlet (inlet_feed,
outlet_brine, and outlet_vapor, respectively).

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Overall heat transfer coefficient", ":math:`U`", "U", "None", ":math:`W/\left(m^2K\right)`"
   "Heat transfer area", ":math:`A`", "area", "None", ":math:`m^2`"
   "Approach temperature in", ":math:`\Delta T_{in}`", "delta_temperature_in", "None", ":math:`K`"
   "Approach temperature out", ":math:`\Delta T_{out}`", "delta_temperature_in", "None", ":math:`K`"
   "Log-mean temperature difference", ":math:`LMTD`", "lmtd", "None", ":math:`K`"
   "Evaporator heat requirement", ":math:`Q_{evap}`", "heat_transfer", "None", ":math:`W`"

Compressor Variables
^^^^^^^^^^^^^^^^^^^^

The condenser consists of 1 ControlVolume0DBlock. Pressure differential and work variables are constructed on the
control volume.

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Pressure ratio", ":math:`PR`", "pressure_ratio", "None", ":math:`W/\left(m^2K\right)`"
   "Isentropic efficiency", ":math:`\eta`", "efficiency", "None", ":math:`\text{dimensionless}`"

Condenser Variables
^^^^^^^^^^^^^^^^^^^

The condenser consists of 1 ControlVolume0DBlock. No additional variables are constructed outside of those on the
control volume.

.. _MVC_equations:

Equations
-----------

Evaporator Equations
^^^^^^^^^^^^^^^^^^^^

.. csv-table::
   :header: "Description", "Equation"

   "Mass balance", ":math:`\dot{m}_{feed} = \dot{m}_{vapor}^{vap}+\dot{m}_{brine}^{liq}`"
   "Energy balance", ":math:`H_{vapor}^{vap}+H_{brine}^{vap}-H_{feed} = Q_{evap}`"
   "Vapor temperature", ":math:`T_{vapor} = T_{brine}`"
   "Log-mean temperature difference (Chen, 1987)", ":math:`LMTD = \left(\frac{1}{2}\left(\Delta T_{in}+\Delta T_{out}\right)\Delta T_{in}\Delta T_{out}\right)^\frac{1}{3}`"
   "Evaporator heat requirement", ":math:`Q_{evap} = UA(LMTD)`"
   "Approach temperature in*", ":math:`\Delta T_{in} = T_{condenser,in}-T_{brine}`"
   "Approach temperature out*", ":math:`\Delta T_{out} = T_{condenser,out}-T_{brine}`"
   "Heat transfer balance*", ":math:`Q_{evap} = -Q_{cond}`"

\*Equations are coupled with the condenser through the ``connect_to_condenser`` method.

Compressor Equations
^^^^^^^^^^^^^^^^^^^^

.. csv-table::
   :header: "Description", "Equation"

   "Pressure ratio", ":math:`P_{out} = (PR)P_{in}`"
   "Isentropic temperature", ":math:`T_{out}=T_{in}(PR)^{1-\frac{1}{\gamma}}`"
   "Efficiency", ":math:`\eta\left(H_{out}-H_{in}\right) = H_{isentropic,out}-H_{in}`"

Condenser Equations
^^^^^^^^^^^^^^^^^^^

The condenser performance is related through the equations denoted by the footnote in the evaporator section.

.. csv-table::
   :header: "Description", "Equation"

   "Complete condensation condition", ":math:`P_{out} >= P_{sat,out}\left(T_{out}\right)`"

Code Documentation
-------------------

* :mod:`watertap.unit_models.mvc.components.complete_condenser`
* :mod:`watertap.unit_models.mvc.components.compressor`
* :mod:`watertap.unit_models.mvc.components.evaporator`

References
-----------

Chen, J. J. J. (1987). Comments on improvements on a replacement for the logarithmic mean. Chemical Engineering Science,
42(10), 2488–2489. https://doi.org/10.1016/0009-2509(87)80128-8
