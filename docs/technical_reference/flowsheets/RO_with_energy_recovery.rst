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

Degrees of Freedom 
------------------

DOF details

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
   "TDS concentration", "35", ":math:`\text{g/L}`"
   "TSS concentration", "0.03", ":math:`\text{g/L}`"
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
   


Code Documentation
------------------



References
----------