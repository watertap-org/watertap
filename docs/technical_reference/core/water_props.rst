Ideal Water Properties
======================

.. index::
   pair: watertap.core.zero_order_properties;WaterParameterBlock

.. currentmodule:: watertap.core.zero_order_properties

The ideal water properties module contains a simple property package for saline waters which is intended for use with the WaterTap zero-order unit model library. The ideal water property package can be used in a flowsheet as shown below:

.. doctest::

  import pyomo.environ as pyo # Pyomo environment

  from idaes.core import FlowsheetBlock

  # Import ideal water property package from watertap.core
  from watertap.core.zero_order_properties import WaterParameterBlock

  # Create a flowsheet
  m = pyo.ConcreteModel()
  m.fs = FlowsheetBlock(dynamic=False)

  # Add an instance of the ideal water property package with three solutes A, B and C
  m.fs.water_props = ZOParameterBlock(solute_list=["A", "B", "C"])

Package Details
---------------

The ideal water property package assumes that the solution is at approximately ambient conditions and that contributions to most properties from the dissolved solutes are minimal. Thus, most properties of the solution are assumed to be constant and equal to those of water at ambient conditions (these are defined via parameters in the WaterParameterBlock and can be adjusted if required to model different conditions).

The state variables used in the ideal water property package are:

* Volumetric flowrate (`flow_vol`, :math:`Q`, units :math:`m^3/s`),
* Mass concentration of solutes (`conc_mass_comp`, :math:`C`, units :math:`kg/m^3`),
* Pressure (`pressure`, :math:`P`, units :math:`Pa`), and
* Temperature (`temperature`, :math:`T`, units :math:`K`).

The ideal water property package supports a single liquid phase (named "Liq") and automatically includes water (named "H2O") as a solvent in the component list.

Dissolved Solutes
-----------------

The ideal water property package requires users to define the list of dissolved solutes present in the solution, which is done using the `solute_list` configuration argument (as shown above). The solutes defined in the `solute_list` configuration argument are automatically added to the property package component list and concentration terms will be created for each of these.

Class Documentation
-------------------

* :class:`WaterParameterBlock`
* :class:`WaterParameterBlockData`
* :class:`WaterStateBlock`
* :class:`WaterStateBlockData`
