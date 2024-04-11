Osmotically Assisted Reverse Osmosis
====================================

Introduction
------------

Osmotically assisted reverse osmosis (OARO) is a nonevaporative membrane-based desalination technology that can treat
high-salinity brines. Compared to conventional reverse osmosis (RO), a saline sweep is added to reduce the osmotic pressure
difference across the membrane as well as enhance water transport. This OARO flowsheet include numbers of
`OARO <https://watertap.readthedocs.io/en/latest/technical_reference/unit_models/osmotically_assisted_reverse_osmosis_0D.html>`_ units,
pumps, energy recovery devices (ERDs) and a
`RO <https://watertap.readthedocs.io/en/latest/technical_reference/unit_models/reverse_osmosis_0D.html>`_ unit
depending on the number of stages. It can be used to implement techno-economic analyses and optimize costing metrics
with specified number of stages and system water recovery.

Implementation
--------------

Figure 1 shows the process flow diagram for n-stage OARO system.
After pressurized by the primary pump, the incoming feed water enters the feed-side of the first OARO module,
then the water permeates across the membrane and the resulting concentrated brine flows out of the system from the feed-side,
In addition, the permeate water will dilute a saline sweep solution in the permeate-side
and the diluted sweep solution will flow into the feed-side of the next stage of OARO module.
Similarly, the resulting concentrate flows out of the feed-side, gets energy recovered by the ERD and
pressurized by the recycle pump, and flows back into permeate-side of the last stage.
Meantime, the diluted sweep solution in the permeate-side flows into the next stage until finally into a conventional RO module.
Costing relationships for each of the unit models is described in the
`WaterTAP Costing Package <https://watertap.readthedocs.io/en/latest/technical_reference/costing/watertap_costing.html>`_
The flowsheet relies on the following key assumptions:

   * supports steady-state only
   * supports optimization and minimizes levelized cost of water (LCOW) with constraints
   * `NaCl Property Package <https://watertap.readthedocs.io/en/latest/technical_reference/property_models/NaCl.html>`_ is utilized
   * number of stages and system recovery should be specified for optimization


.. figure:: ../../_static/flowsheets/OARO.png
    :width: 800
    :align: center

    Figure 1. OARO flow diagram

Documentation for each of the unit models can be found here:
   * `OARO <https://watertap.readthedocs.io/en/latest/technical_reference/unit_models/osmotically_assisted_reverse_osmosis_0D.html>`_
   * `RO <https://watertap.readthedocs.io/en/latest/technical_reference/unit_models/reverse_osmosis_0D.html>`_

Documentation for the property model can be found here:
    * `NaCl Property Package <https://watertap.readthedocs.io/en/latest/technical_reference/property_models/NaCl.html>`_

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Inlet/outlet", ":math:`x`", "['in', 'out', 'disposal']"
   "Phases", ":math:`p`", "['Liq']"
   "Components", ":math:`j`", "['H2O', 'NaCl']*"

\*Solute depends on the imported property model; example shown here is for the NaCl property model.

Flowsheet Specifications
------------------------
This flowsheet aims to solve optimization problem that minimizes levelized cost of water (LCOW) with specified
number of stages and system mass recovery. There is an extra inequality constraint to govern product quality:

.. csv-table::
   :header: "Description", "Equation"

   "Product Quality", ":math:`M_{out, NaCl} \le 500 \text{mg/L}`"

Degrees of Freedom
------------------
Firstly, the number of stages and water mass recovery of the system need to be specified. In addition,
aside from the inlet feed state variables (i.e., temperature, pressure, and component mass flowrates), each unit has
additional degrees of freedom (DOF) that needs to be specified.

For primary pumps and recycle pumps, there is 1 DOF that requires to be fixed:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units", "Default Value"

   "Pump efficiency", ":math:`\eta`", "efficiency_pump", "[t]", ":math:`\text{dimensionless}`", "0.75"

For ERDs, there are 2 DOF that require to be fixed:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units", "Default Value"

   "Pump efficiency", ":math:`\eta`", "efficiency_pump", "[t]", ":math:`\text{dimensionless}`", "0.9"
   "Outlet pressure", ":math:`P_{out}`", "pressure", "[t]", ":math:`\text{Pa}`", "101325"

For OARO units, typically there are 7 DOF that require to be fixed:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units", "Default Value"

   "Solvent permeability coefficient", ":math:`A`", "A_comp", "[t, j]", ":math:`\text{m/Pa/s}`", "1E-12"
   "Solute permeability coefficient", ":math:`B`", "B_comp", "[t, j]", ":math:`\text{m/s}`", "8E-8"
   "Membrane structural parameter", ":math:`S`", "structural_parameter", "[None]", ":math:`\text{\mu m}`", "1200"
   "Feed-channel height", ":math:`h_{ch,f}`", "feed_side.channel_height", "None", ":math:`\text{m}`", "2E-3"
   "Feed-side spacer porosity", ":math:`\epsilon_{sp,f}`", "feed_side.spacer_porosity", "None", ":math:`\text{dimensionless}`", "0.75"
   "Permeate-channel height", ":math:`h_{ch,p}`", "permeate_side.channel_height", "None", ":math:`\text{m}`", "2E-3"
   "Peremeate-side spacer porosity", ":math:`\epsilon_{sp,p}`", "permeate_side.spacer_porosity", "None", ":math:`\text{dimensionless}`", "0.75"

For RO unit, typically there are 5 DOF that require to be fixed:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units", "Default Value"

   "Solvent permeability coefficient", ":math:`A`", "A_comp", "[t, j]", ":math:`\text{m/Pa/s}`", "4.2E-12"
   "Solute permeability coefficient", ":math:`B`", "B_comp", "[t, j]", ":math:`\text{m/s}`", "3.5E-8"
   "Feed-channel height", ":math:`h_{ch,f}`", "feed_side.channel_height", "None", ":math:`\text{m}`", "2E-3"
   "Feed-side spacer porosity", ":math:`\epsilon_{sp,f}`", "feed_side.spacer_porosity", "None", ":math:`\text{dimensionless}`", "0.75"
   "Permeate pressure", ":math:`P_{p}`", "permeate.pressure", "[t]", ":math:`\text{Pa}`", "101325"
