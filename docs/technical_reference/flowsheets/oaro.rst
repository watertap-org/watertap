Osmotically Assisted Reverse Osmosis
====================================

Introduction
------------

Osmotically assisted reverse osmosis (OARO) is a nonevaporative membrane-based desalination technology that can treat
high-salinity brines. Compared to conventional reverse osmosis (RO), a saline sweep is added to reduce the pressure
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
   * supports optimization minimizes levelized cost of water (LCOW) with constraints
   * `NaCl Property Package <https://watertap.readthedocs.io/en/latest/technical_reference/property_models/NaCl.html>`_ is utilized
   * number of stages and system recovery should be specified for optimization


.. figure:: ../../_static/flowsheets/OARO.png
    :width: 800
    :align: center

    Figure 1. OARO flow diagram