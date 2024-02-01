Extended Benchmark Simulation Model No.2
========================================

Introduction
------------

Like the conventional Benchmark Simulation Model No.2 (BSM2)`BSM2 <https://watertap.readthedocs.io/en/latest/technical_reference/flowsheets/BSM2.html>`_,
extended BSM2 is an industry benchmark for modeling a full biological wastewater
treatment plant that includes a primary clarifier, the activated sludge process, and an anaerobic digester.
These unit processes are driven by biological reaction models that relate soluble and particulate wastewater
components to their respective process rate equations. The main difference between conventional and extended BSM2
is that the latter uses modified ADM1 and ASM2d property packages as opposed to conventional ADM1 and ASM1. These modifications allow
key components like phosphorus, magnesium and calcium to be tracked throughout the system, which are essential in order to
accurately modeling certain novel technologies that can incorporated into BSM2. Thus, while this flowsheet can simply be used to
simulate and run techno-economic analyses on the operation of a conventional wastewater treatment plant,
an additional layer of utility can be derived from using BSM2 as a baseline for comparing alternative plant
configurations to a well-established standard and/or amongst the variations themselves by adding, removing,
or modifying unit processes using WaterTAP's flexible modeling capabilities.

Implementation
--------------

Figure 1 shows the process flow diagram for BSM2 where influent wastewater is fed
to a primary clarifier (primary treatment); the effluent is then passed to a series of activated sludge
reactors and a secondary clarifier (secondary treatment). CSTRs are used to model the anoxic reactors in the activated
sludge process and CSTRs with injection (which accounts for aeration tanks) are used to model the aerobic reactors.
Finally, the sludge is passed through a thickener and sent to the anaerobic digester. The anaerobic digester processes
the sludge to produce a biogas stream and residual sludge stream that passes through a dewatering unit which recycles
liquid to the headworks of the plant while sludge is released for disposal. Costing relationships for each of the unit
models is described in the `WaterTAP Costing Package <https://watertap.readthedocs.io/en/latest/technical_reference/costing/watertap_costing.html>`_
The flowsheet relies on the following key assumptions:

   * supports steady-state only
   * property and reaction package are provided for the activated sludge model (ASM)
   * property and reaction package are provided for the anaerobic digester model (ADM)
   * interfaces are provided to convert between the properties of ASM and ADM

.. figure:: ../../_static/flowsheets/BSM2.png
    :width: 800
    :align: center

    Figure 1. BSM2 flowsheet

Documentation for each of the unit models can be found here:
    * `Thickener <https://watertap.readthedocs.io/en/latest/technical_reference/unit_models/thickener.html>`_
    * `Anaerobic digester <https://watertap.readthedocs.io/en/latest/technical_reference/unit_models/anaerobic_digester.html>`_
    * `Dewatering unit <https://watertap.readthedocs.io/en/latest/technical_reference/unit_models/dewatering_unit.html>`_
    * `CSTR <https://idaes-pse.readthedocs.io/en/latest/reference_guides/model_libraries/generic/unit_models/cstr.html>`_
    * `ADM1 to ASM2d Translator <https://watertap.readthedocs.io/en/latest/technical_reference/unit_models/translators/translator_adm1_asm2d.html>`_
    * `ASM2d to ADM1 Translator <https://watertap.readthedocs.io/en/latest/technical_reference/unit_models/translators/translator_asm2d_adm1.html>`_
    * CSTR with injection
    * Primary clarifier
    * Secondary clarifier

Documentation for each of the property models can be found here:
    * `Modified ASM2d <https://watertap.readthedocs.io/en/latest/technical_reference/property_models/modified_ASM2D.html>`_
    * `Modified ADM1 <https://watertap.readthedocs.io/en/latest/technical_reference/property_models/modified_ADM1.html>`_

Future Refinements
------------------

The following modifications to extended BSM2 are planned for development:
    * Improving costing relationships in terms of detail, completeness, and reasonable validity
    * Accounting for temperature-dependence in the oxygen mass transfer coefficient (KLa) and oxygen concentration at saturation
    * Adding thermal energy requirements to the anaerobic digester and refining energy consumption estimates for units collectively
    * Accounting for mineral precipitation reactions
    * Accounting for ion speciation and activity
    * Accounting for sulfur components
    * Accounting for iron components
    * Replacing the ideal-separator formulation in the secondary clarifier with the widely used double-exponential settling model (i.e., the Takacs model)
