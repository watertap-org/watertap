Reverse Osmosis with Energy Recovery
====================================

Introduction
------------

This flowsheet represents a reverse osmosis (RO) with energy recovery device (ERD). 
The flowsheet is similar to the desalination block of the full-scale seawater RO treatment facility flowsheet,
but provides different options for optimization of the process.

Implementation
--------------

This flowsheet uses several different modeling features available in WaterTAP, including:


WaterTAP costing package
    * :doc:`/technical_reference/costing/watertap_costing`
Unit model costing packages
    * :doc:`/technical_reference/costing/detailed_unit_model_costing`
Sodium chloride (NaCl) property package
    * :doc:`/technical_reference/property_models/NaCl`
Reverse osmosis model
    * :doc:`/technical_reference/unit_models/reverse_osmosis_0D`
Pressure exchanger model
    * :doc:`/technical_reference/unit_models/pressure_exchanger`
IDAES Product blocks
    * :doc:`idaes:reference_guides/model_libraries/generic/unit_models/product`
IDAES Separator blocks
    * :doc:`idaes:reference_guides/model_libraries/generic/unit_models/separator`
IDAES Mixer blocks
    * :doc:`idaes:reference_guides/model_libraries/generic/unit_models/mixer`

The demonstration file itself contains several core functions that are used to build, specify, initialize, solve, and (optionally) optimize the model, 
as well as helper functions that group these core functions together for convenience.

* Creating and instantiating the model using ``build()``:
    
    This function will create the core components and structure of the flowsheet including unit, property, and costing models.
    Different components will be added depending on the type of ERD used (dictated by the ``erd_type`` argument).
    An ``Arc`` is used to connect the different component flows to one another and default scaling factors are set.

* Specify the operating conditions with ``set_operating_conditions()``:

    Influent component flows and RO pump operating pressures are set with this function.
    The user can set via keyword arguments the volumetric flow rate (``flow_vol``), mass fraction of NaCl (``salt_mass_conc``),
    a target fractional mass-based water recovery ``water_recovery``, and a fractional over pressure (``over_pressure``). 
    Over pressure is the fractional increase in pressure over of the brine osmotic pressure to use as the
    operating pressure for the RO process.




Degrees of Freedom 
------------------

The degrees of freedom (DOF) for the flowsheet can change depending on model configuration options.
For either ``pump_as_turbine`` or ``pressure_exchanger`` as ``erd_type``, there are 15 DOF. Running
the model with ``no_ERD`` results in 13 DOF.

* Influent conditions (component flows, temperature, pressure)
* RO membrane properties
* RO operating pressure
* Pump and ERD efficiencies

Passing any model build to the provided function ``set_operating_conditions()`` will result in a model with zero DOF.

Flowsheet Specifications
------------------------


The influent conditions are defined from the case study used to develop this flowsheet. 
Additionally, some unit models have case-specific operating conditions.
The influent conditions and case-specific operating conditions for certain unit models are presented in the following table,
including the different build options for ``erd_type``:

.. csv-table::
   :header: "Description", "Default Value", "Units"

    **Influent Conditions**
   "Volumetric flow rate", "1e-3", ":math:`\text{m}^3/\text{s}`"
   "TDS mass fraction", "0.035", ":math:`\text{dimensioneless}`"
   "Temperature", "298", ":math:`\text{K}`"
   "Pressure", "101325", ":math:`\text{Pa}`"
   
   **Desalination**
   "RO water permeability coefficient", "4.2e-12", ":math:`\text{m/Pa/s}`"
   "RO salt permeability coefficient", "3.5e-8", ":math:`\text{m/s}`"
   "RO spacer porosity", "0.85", ":math:`\text{dimensionless}`"
   "RO channel height", "1e-3", ":math:`\text{m}`"
   "RO membrane width per stage", "5", ":math:`\text{m}`"
   "RO total membrane area per stage", "50", ":math:`\text{m}^2`"
   "RO permeate side pressure", "101325", ":math:`\text{Pa}`"
   "Pump 1 efficiency", "0.8", ":math:`\text{dimensionless}`"
   "Pump 1 operating pressure", "70e5", ":math:`\text{Pa}`"
   
   *if* ``erd_type == "pressure_exchanger"``
   "Pressure exchanger efficiency", "0.95", ":math:`\text{dimensionless}`"
   "Pump 2 efficiency", "0.8", ":math:`\text{dimensionless}`"
   
   *if* ``erd_type == "pump_as_turbine"``
   "Energy recovery device pump efficiency", "0.95", ":math:`\text{dimensionless}`"
   "Energy recovery device permeate side pressure", "101325", ":math:`\text{Pa}`"
   


Code Documentation
------------------



References
----------