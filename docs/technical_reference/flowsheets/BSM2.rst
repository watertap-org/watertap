Benchmark Simulation Model No.2
===============================

Introduction
------------

The Benchmark Simulation Model No.2 (BSM2) is an industry benchmark for modeling a full biological wastewater
treatment plant that includes a primary clarifier, the activated sludge process, and an anaerobic digester. These
unit processes are driven by biological reaction models that relate soluble and particulate wastewater
components to their respective process rate equations. While this BSM2 flowsheet can simply be used to
simulate and run techno-economic analyses on the operation of a conventional wastewater treatment plant,
an additional layer of utility can be derived from using BSM2 as a baseline for comparing alternative plant
configurations to a well-established standard and/or amongst the variations themselves by adding, removing,
or modifying unit processes using WaterTAP's flexible modeling capabilities.

Implementation
--------------

Figure 1 shows the process flow diagram for BSM2 where influent wastewater is fed
to a primary clarifier (primary treatment); the effluent is then passed to a series of activated sludge
reactors and a secondary clarifier (secondary treatment), and finally the sludge is passed through a thickener and
sent to the anaerobic digester. The anaerobic digester processes the sludge to produce
a biogas stream and residual sludge stream that passes through a dewatering unit which recycles liquid to
the headworks of the plant while sludge is released for disposal. The flowsheet relies on the following key assumptions:

   * supports steady-state only
   * property and reaction package are provided for the activated sludge model (ASM)
   * property and reaction package are provided for the anaerobic digester model (ADM)
   * interfaces are provided to convert between the properties of ASM and ADM

.. figure:: ../../_static/flowsheets/BSM2.png
    :width: 400
    :align: center

    Figure 1. BSM2 flowsheet

Documentation for each of the unit models can be found here:
    * `Thickener <https://watertap.readthedocs.io/en/latest/technical_reference/unit_models/thickener.html>`_
    * `Anaerobic digester <https://watertap.readthedocs.io/en/latest/technical_reference/unit_models/anaerobic_digester.html>`_
    * `Dewatering unit <https://watertap.readthedocs.io/en/latest/technical_reference/unit_models/dewatering_unit.html>`_

Documentation for each of the property models can be found here:
    * `ASM1 <https://watertap.readthedocs.io/en/latest/technical_reference/property_models/ASM1.html>`_
    * `ADM1 <https://watertap.readthedocs.io/en/latest/technical_reference/property_models/ADM1.html>`_