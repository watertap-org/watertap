Seawater RO Desalination
========================

Introduction
------------

This flowsheet represents a full-scale seawater reverse osmosis treatment facility.
The flowsheet includes pretreatment, desalination, post-treatment, and waste handling unit processes available in WaterTAP.
Unit models used on this flowsheet span the spectrum from simple to detailed.


Implementation
--------------

This flowsheet uses several different modeling features available in WaterTAP, including:


WaterTAP costing package
    * :doc:`/technical_reference/costing/watertap_costing`
Unit model costing packages
    * :doc:`/technical_reference/costing/detailed_unit_model_costing`
Zero order unit and costing models
    * :doc:`/technical_reference/unit_models/zero_order_unit_models/index`
Zero order costing package
    * :doc:`/technical_reference/costing/zero_order_costing`
Zero order property package
    * :doc:`/technical_reference/core/water_props`
Seawater property package
    * :doc:`/technical_reference/property_models/seawater`
Reverse osmosis model
    * :doc:`/technical_reference/unit_models/reverse_osmosis_0D`
Pressure exchanger model
    * :doc:`/technical_reference/unit_models/pressure_exchanger`
IDAES Translator blocks
    * :doc:`idaes:reference_guides/model_libraries/generic/unit_models/translator`
IDAES Product blocks
    * :doc:`idaes:reference_guides/model_libraries/generic/unit_models/product`
IDAES Separator blocks
    * :doc:`idaes:reference_guides/model_libraries/generic/unit_models/separator`
IDAES Mixer blocks
    * :doc:`idaes:reference_guides/model_libraries/generic/unit_models/mixer`
    
The demonstration file itself contains several core functions that are used to build, specify, initialize, and solve the model, as well as
some helper functions that group these core functions together for convenience. Building and solving the flowsheet proceeds in six steps:

1. Creating and instantiating the model using ``build()``:

    This function will create the core components and structure of the flowsheet. The keyword argument ``erd_type`` dictates the type
    of energy recovery device to be used, with options ``pressure_exchanger`` or ``pump_as_turbine``.
    First, the ``FlowsheetBlock``, ``Database``, property models, are created. The zero order property models require the user
    to provide a ``solute_list`` (e.g., TDS and TSS), while the seawater property model is pre-populated with TDS as the only solute.
    Separate ``Block`` are created to contain all the unit models required to model the pretreatment, desalination, and post-treatment
    parts of the treatment train:

        * Pre-treatment (``m.fs.pretreatment``): comprised entirely of zero order models, this block contains the intake, chemical addition, media, and cartridge filtration unit models.
        * Desalination (``m.fs.desalination``): contains all the unit models needed to represent the pumping, reverse osmosis (RO), and energy recovery device (ERD) processes.
        * Post-treatment (``m.fs.posttreatment``): includes post desalination disinfection, remineralization, and storage unit models.

    Outside of these, the there is a feed, distribution, landfill, and disposal block that are placed directly on the flowsheet.
    ``Translator`` blocks are added with appropriate constraints and ``Arc`` are used to connect the unit processes in the proper order.
    Finally, default scaling factors are set and scaling factors are calculated for all variables.

2. Specify the operating conditions with ``set_operating_conditions()``:

    This function begins by specifying the inlet conditions as outlined in Table 2. Then, starting with the ``pretreatment`` block, the operating 
    conditions for each unit model are set according to Table 3.

3. Initialize the unit and costing models with ``initialize_system()``:

    Starting with the ``Feed`` block and continuing sequentially through each part of the treatment train, this function sets the initial condition
    for all the unit models on the flowsheet and propagates ``Arcs`` that connect each block either to ``Translator`` blocks or treatment blocks. 
    Each block is initialized and solved by using the local ``solve()`` function to ensure each block solves optimally before trying to solve the next.
    Note that the order in which blocks/unit models are solved/initialized in WaterTAP is important because the initial conditions are *only* set 
    for the ``Feed`` block. For these conditions to cascade to downstream unit models, and for the downstream unit models to include the impacts of upstream 
    process (e.g., component removal), sequential initialization is necessary. Thus, initialization of this flowsheet proceeds as follows:

        #. ``Feed`` block where initial flow rates and solute concentrations are set.
        #. ``pretreatment`` block 
        #. translator block from ``pretreatment`` to ``desalination`` (i.e., ``m.fs.tb_prtrt_desal``)
        #. ``desalination`` block
        #. translator block from ``desalination`` to ``posttreatment`` (i.e., ``m.fs.tb_desal_psttrt``)
        #. ``posttreatment`` block

4. Add the system- and unit-level costing packages with ``add_costing()`` and initialize with ``initialize_costing()``:

    Because of the nature of the unit models used in this flowsheet (i.e., both zero order and detailed models), two separate system-level costing packages are required. 
    ``m.fs.zo_costing = ZeroOrderCosting()`` is used to aggregate costs for zero-order models, and ``m.fs.ro_costing = WaterTAPCosting`` is for the more detailed desalination models. 
    The costing block for each unit model is ``UnitModelCostingBlock`` that points to a system-level aggregation costing package via the configuration keyword ``flowsheet_costing_block``.
    Each system-level costing package has a ``.cost_process()`` method that is called to aggregate unit level costs and calculate overall process costs.
    To aggregate results from both costing packages, a separate ``Expression`` is created for ``total_capital_cost`` and ``total_operating_cost``, and each of these are used
    to calculate the ``LCOW``. Finally, like the unit models, the costing packages are initialized.

5. Solve the entire flowsheet and display final results with ``display_results()``:

    After building, specifying, initializing, and costing the models, the flowsheet is solved a final time with the local ``solve()`` function
    and the ``.report()`` method is called for each unit model using ``display_results()``.

The local function ``build_flowsheet()`` combines steps 1 and 2. The local function ``solve_flowsheet()`` combines 3, 4, and 5.

.. note::

    ``Translator`` blocks are used in flowsheets when more than one property package is used in different parts of the flowsheet.
    In this example, the zero order property package contains state variables that are only indexed by component (e.g., ``conc_mass_comp['TDS']``)
    while the seawater property package contains state variables indexed by both phase and component (e.g., ``conc_mass_phase_comp['Liq', 'TDS']``).
    The translator blocks in this flowsheet are used to communicate properties between unit models that use these two property packages and simply
    say, e.g., ``conc_mass_comp['TDS'] = conc_mass_phase_comp['Liq', 'TDS']``.

Pre-Treatment
^^^^^^^^^^^^^

Figure 1 presents the process flow diagram for ``m.fs.pretreatment``. The first unit model on this block, ``intake``, is connected to the
flowsheet level ``Feed`` block.

.. figure:: ../../_static/flowsheets/SW_RO_fs_pretreat.png
    :width: 800
    :align: center

    Figure 1: Process flow diagram for pre-treatment block.

Desalination
^^^^^^^^^^^^

Figure 2 presents the process flow diagram for ``m.fs.desalination`` if ``erd_type == "pressure_exchanger"``.
Figure 3 presents the process flow diagram for ``m.fs.desalination`` if ``erd_type == "pump_as_turbine"``.
In either case, the first unit model on this block is connected to the flowsheet level translator block ``tb_prtrt_desal``.

.. figure:: ../../_static/flowsheets/SW_RO_fs_desal_PXR.png
    :width: 800
    :align: center

    Figure 2: Process flow diagram for desalination block for pressure exchanger ERD.

.. figure:: ../../_static/flowsheets/SW_RO_fs_desal_turbine.png
    :width: 600
    :align: center

    Figure 3: Process flow diagram for desalination block for pump-as-turbine ERD.

Post-Treatment
^^^^^^^^^^^^^^

Figure 4 presents the process flow diagram for ``m.fs.posttreatment``. The first unit model on this block is connected to the
flowsheet level translator block ``tb_desal_psttrt``.

.. figure:: ../../_static/flowsheets/SW_RO_fs_posttreat.png
    :width: 800
    :align: center

    Figure 4: Process flow diagram for post-treatment block.

Full Flowsheet
^^^^^^^^^^^^^^

Figure 5 presents the process flow diagram for the entire flowsheet for both ERD options.

.. figure:: ../../_static/flowsheets/SW_RO_fs_full.png
    :width: 1500
    :align: center

    Figure 5: Process flow diagram for entire flowsheet.


Degrees of Freedom
------------------

The degrees of freedom (DOF) for the flowsheet can change depending on model configuration options.
For either ``erd_type``, after building the flowsheet with the provided ``build()`` function, there are 58 DOF.
Many of these DOF are related to the not-yet-propagated Arcs between each of the unit process blocks.
Others are related to influent conditions and specific operating conditions for each unit model that must be specified by the user.
The operating conditions for this demonstration are provided in the following section, and include:

* Influent conditions (component flows, temperature, pressure)
* Chemical doses
* RO membrane properties
* RO operating pressure
* Pump and ERD efficiencies
* Storage duration

Passing the default build to the provided function ``set_operating_conditions()`` will result in a model with zero DOF.

Flowsheet Specifications
------------------------

The influent conditions are defined from the case study used to develop this flowsheet. 
Additionally, some unit models have case-specific operating conditions.
The influent conditions and case-specific operating conditions for certain unit models are presented in the following table,
including the different build options for ``erd_type``:

.. csv-table::
   :header: "Description", "Value", "Units", "Flowsheet Model Name"

    **Influent Conditions**
   "Volumetric flow rate", "7.05", ":math:`\text{MGD}`"
   "TDS :sup:`1` concentration", "35", ":math:`\text{g/L}`"
   "TSS :sup:`2` concentration", "0.03", ":math:`\text{g/L}`"
   "Temperature", "298", ":math:`\text{K}`"
   "Pressure", "100000", ":math:`\text{Pa}`"

   **Pre-Treatment**
   "Ferric chloride dose", "20", ":math:`\text{mg/L}`", "``m.fs.pretreatment.ferric_chloride_addition``"
   "Storage tank 1 storage time", "2", ":math:`\text{hr}`", "``m.fs.pretreatment.storage_tank_1``"
   
   **Desalination**
   "RO water permeability coefficient", "4.2e-12", ":math:`\text{m/Pa/s}`", "``m.fs.desalination.RO``"
   "RO salt permeability coefficient", "3.5e-8", ":math:`\text{m/s}`", "``m.fs.desalination.RO``"
   "RO spacer porosity", "0.97", ":math:`\text{dimensionless}`", "``m.fs.desalination.RO``"
   "RO channel height", "1e-3", ":math:`\text{m}`", "``m.fs.desalination.RO``"
   "RO membrane width per stage", "1000", ":math:`\text{m}`", "``m.fs.desalination.RO``"
   "RO total membrane area per stage", "13914", ":math:`\text{m}^2`", "``m.fs.desalination.RO``"
   "RO permeate side pressure", "101325", ":math:`\text{Pa}`", "``m.fs.desalination.RO``"
   "Pump 1 efficiency", "0.8", ":math:`\text{dimensionless}`", "``m.fs.desalination.P1``"
   "Pump 1 operating pressure", "70e5", ":math:`\text{Pa}`", "``m.fs.desalination.P1``"
   
   *if* ``erd_type == "pressure_exchanger"``
   "Pressure exchanger efficiency", "0.95", ":math:`\text{dimensionless}`", "``m.fs.desalination.PXR``"
   "Pump 2 efficiency", "0.8", ":math:`\text{dimensionless}`", "``m.fs.desalination.P2``"
   
   *if* ``erd_type == "pump_as_turbine"``
   "Energy recovery device pump efficiency", "0.95", ":math:`\text{dimensionless}`", "``m.fs.desalination.ERD``"
   "Energy recovery device permeate side pressure", "101325", ":math:`\text{Pa}`", "``m.fs.desalination.ERD``"
   
   **Post-Treatment**
   "Anti-scalant dose", "5", ":math:`\text{mg/L}`", "``m.fs.posttreatment.anti_scalant_addition``"
   "Lime dose", "2.3", ":math:`\text{mg/L}`", "``m.fs.posttreatment.lime_addition``"
   "Storage tank 2 storage time", "1", ":math:`\text{hr}`", "``m.fs.posttreatment.storage_tank_2``"
   "Storage tank 3 storage time", "1", ":math:`\text{hr}`", "``m.fs.posttreatment.storage_tank_3``"
   "UV/AOP :sup:`3` reduction equivalent dose", "350", ":math:`\text{mJ/}\text{cm}^2`", "``m.fs.posttreatment.uv_aop``"
   "UV/AOP :sup:`3` UV transmittance", "0.95", ":math:`\text{dimensionless}`", "``m.fs.posttreatment.uv_aop``"

.. note::
   
   :sup:`1`  TDS = total dissolved solids
   |
   :sup:`2`  TSS = total suspended solids
   |
   :sup:`3`  UV = Ultraviolet; AOP = Advanced oxidation process
