.. _OARO_flowsheet:

Osmotically Assisted Reverse Osmosis
====================================

Introduction
------------

Osmotically assisted reverse osmosis (OARO) is a non-evaporative membrane-based desalination technology that can treat
high-salinity brines. Compared to conventional reverse osmosis (RO), a saline sweep is added to reduce the osmotic pressure
difference across the membrane as well as enhance water transport. This OARO flowsheet includes numbers of
:ref:`OARO <OARO_0D>` units,
pumps, energy recovery devices (ERDs) and a
:ref:`RO <RO_0D>` unit
depending on the number of stages. It can be used to implement techno-economic analyses and optimize costing metrics
with specified number of stages and system water recovery.

Implementation
--------------

Figure 1 shows the process flow diagram for an n-stage OARO system.
After pressurized by the primary pump, the incoming feed water enters the feed-side of the first OARO module,
then the water permeates across the membrane and the resulting concentrated brine flows out of the system from the feed-side,
In addition, the permeate water will dilute a saline sweep solution in the permeate-side
and the diluted sweep solution will flow into the feed-side of the next stage of OARO module.
Similarly, the resulting concentrate flows out of the feed-side, gets energy recovered by the ERD and
pressurized by the recycle pump, and flows back into permeate-side of the last stage.
Meanwhile, the diluted sweep solution in the permeate-side flows into the next stage and finally into a conventional RO module.
Costing relationships for each of the unit models is described in the
:ref:`WaterTAP Costing Package Package <watertap_costing>`.
The flowsheet relies on the following key assumptions:

   * supports steady-state only
   * supports optimization and minimizes levelized cost of water (LCOW) with constraints
   * :ref:`NaCl Property Package <nacl>` is utilized
   * number of stages and system recovery should be specified for optimization


.. figure:: ../../_static/flowsheets/OARO.png
    :width: 800
    :align: center

    Figure 1. OARO flow diagram

Documentation for each of the unit models can be found here:
   * :ref:`OARO <OARO_0D>`
   * :ref:`RO <RO_0D>`

Documentation for the property model can be found here:
    * :ref:`NaCl Property Package <nacl>`

Documentation for the costing relationships can be found below.
    * :ref:`WaterTAP Costing Package <watertap_costing>`

This flowsheet aims to solve optimization problem that minimizes levelized cost of water (LCOW) with specified
number of stages and system mass recovery.
LCOW can be represented by the following equation where :math:`Q` represents product volumetric flow, :math:`f_{crf}` represents capital recovery factor
:math:`C_{cap,tot}` represents total capital cost, :math:`C_{op,tot}` represents total operating cost, and
:math:`f_{util}` represents the utilization factor:

    .. math::

        LCOW_{Q} = \frac{f_{crf}   C_{cap,tot} + C_{op,tot}}{f_{util} Q}

Degrees of Freedom
------------------
Firstly, the number of stages and water mass recovery of the system need to be specified. In addition, the following variables needs to be specified based on the default settings:
   * feed water conditions (flow, temperature, pressure, component concentrations)
   * pump efficiency of primary pumps and recycle pumps
   * ERD pump efficiency and outlet pressure
   * OARO solvent and solute permeability coefficients, membrane structural parameter, channel height and spacer porosity of both feed-side and permeate-side
   * RO solvent and solute permeability coefficients, feed-side channel height and spacer porosity, and permeate pressure

Flowsheet Specifications
------------------------

.. csv-table::
   :header: "Description", "Units", "Value"

   "**Primary pumps**"
   "Pump efficiency", ":math:`\text{dimensionless}`", "0.75"
   "**Recycle pumps**"
   "Pump efficiency", ":math:`\text{dimensionless}`", "0.75"
   "**ERDs**"
   "Pump efficiency", ":math:`\text{dimensionless}`", "0.75"
   "Outlet pressure", ":math:`\text{Pa}`", "101325"
   "**OAROs***"
   "Solvent permeability coefficient", ":math:`\text{m/Pa/s}`", "1E-12"
   "Solute permeability coefficient", ":math:`\text{m/s}`", "8E-8"
   "Membrane structural parameter", ":math:`\mu \text{m}`", "1200"
   "Feed-channel height", ":math:`\text{m}`", "2E-3"
   "Feed-side spacer porosity", ":math:`\text{dimensionless}`", "0.75"
   "Permeate-channel height", ":math:`\text{m}`", "2E-3"
   "Peremeate-side spacer porosity", ":math:`\text{dimensionless}`", "0.75"
   "**RO***"
   "Solvent permeability coefficient", ":math:`\text{m/Pa/s}`", "4.2E-12"
   "Solute permeability coefficient", ":math:`\text{m/s}`", "3.5E-8"
   "Feed-channel height", ":math:`\text{m}`", "2E-3"
   "Feed-side spacer porosity", ":math:`\text{dimensionless}`", "0.75"
   "Permeate pressure", ":math:`\text{Pa}`", "101325"

\*Settings for :ref:`OARO <OARO_0D>`
and :ref:`RO <RO_0D>`
can vary depending on the configurations.


Additional Variables
--------------------

.. csv-table::
   :header: "Description", "Symbol", "Value", "Units"

   "Maximum product concentration", ":math:`M_{out, max}`", "500", ":math:`\text{mg/L}`"

Additional Constraints
----------------------

There is an extra inequality constraint to ensure the product quality:

.. csv-table::
   :header: "Description", "Equation"

   "Product Quality", ":math:`M_{out, NaCl} \le M_{out, max}`"