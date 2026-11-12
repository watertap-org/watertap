Zero Order (ZO) Properties
===========================

.. index::
   pair: watertap.property_models.zero_order_prop_pack;ZOParameterBlock

.. currentmodule:: watertap.property_models.zero_order_prop_pack

The ZO (zero-order) properties module contains a simple property package for saline waters which is intended for use with the WaterTAP zero-order unit model library. 
The ZO property package can be used in a flowsheet as shown below:

.. doctest::

  import pyomo.environ as pyo # Pyomo environment

  from idaes.core import FlowsheetBlock

  # Import ZO property package
  from watertap.property_models import ZOParameterBlock

  # Create a flowsheet
  m = pyo.ConcreteModel()
  m.fs = FlowsheetBlock(dynamic=False)

  # Add an instance of the ZO property package with three solutes: A, B and C
  m.fs.zo_props = ZOParameterBlock(solute_list=["A", "B", "C"])

Package Details
---------------

The ZO property package assumes that the solution is at approximately ambient conditions and that contributions to most properties from the dissolved solutes are minimal. 
Thus, most properties of the solution are assumed to be constant and equal to those of water at ambient conditions (these are defined via parameters in the ZOParameterBlock and can be adjusted if required to model different conditions).

The state variables used in the ZO property package are:

* Volumetric flowrate (`flow_vol`, :math:`Q`, units :math:`m^3/s`),
* Mass concentration of solutes (`conc_mass_comp`, :math:`C`, units :math:`kg/m^3`),
* Pressure (`pressure`, :math:`P`, units :math:`Pa`), and
* Temperature (`temperature`, :math:`T`, units :math:`K`).

The ZO property package supports a single liquid phase (named "Liq") and automatically includes water (named "H2O") as a solvent in the component list.

Dissolved Solutes
-----------------

The ZO property package requires users to define the list of dissolved solutes present in the solution, which is done using the `solute_list` configuration argument (as shown above). The solutes defined in the `solute_list` configuration argument are automatically added to the property package component list and concentration terms will be created for each of these.

Class Documentation
-------------------

* :class:`ZOParameterBlock`
* :class:`ZOParameterData`
* :class:`ZOStateBlock`
* :class:`ZOStateBlockData`
