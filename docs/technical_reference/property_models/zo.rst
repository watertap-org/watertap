Zero Order Property Package
============================

.. code-block:: python

   from watertap.property_models import ZOParameterBlock

.. index::
   pair: watertap.property_models.zero_order_prop_pack;ZOParameterBlock

.. currentmodule:: watertap.property_models.zero_order_prop_pack


The zero order (ZO) property package is a simple property model that is intended for use with the WaterTAP :ref:`zo_unit_models_home`.
The ZO property package can be used in a flowsheet as shown below:

.. doctest::

  from pyomo.environ import ConcreteModel

  from idaes.core import FlowsheetBlock

  # Import ZO property package
  from watertap.property_models import ZOParameterBlock

  # Create a flowsheet
  m = ConcreteModel()
  m.fs = FlowsheetBlock(dynamic=False)

  # Add an instance of the ZO property package with three solutes A, B and C
  m.fs.zo_props = ZOParameterBlock(solute_list=["A", "B", "C"])

Package Details
---------------

The ZO property package assumes that the solution is at approximately ambient conditions and that contributions to most properties from the dissolved solutes are minimal. 
Thus, most properties of the solution are assumed to be constant and equal to those of water at ambient conditions (these are defined via parameters in the ZOParameterBlock and can be adjusted if required to model different conditions).

The only state variable in the ZO property package is the mass flow rate of the solutes, `flow_mass_comp`.
This includes water as `H2O` and any additional solutes defined by the user via the `solute_list` configuration argument when creating an instance of the ZOParameterBlock.

The ZO property package supports a single liquid phase (named "Liq") and automatically includes water 
(named "H2O") as a solvent in the component list.

Additional properties available in the ZO property package include:

.. csv-table::
   :header: "Description", "Variable", "Equation"

   "Mass flow rate", "``flow_mass_comp``", ":math:`M_j`"
   "Mass density", "``dens_mass``", ":math:`\rho`"
   "Component mass concentration", "``conc_mass_comp``", ":math:`C_j = \frac{M_j}{\sum_{j} M_j} \times \rho`"
   "Volumetric flowrate", "``flow_vol``", ":math:`Q = \frac{\sum_{j} M_j}{\rho}`"
   "Dynamic viscosity", "``visc_d``", ":math:`\mu = \mu_{H2O}`"


Dissolved Solutes
-----------------

The zero order water property package requires users to define the list of dissolved solutes present in the solution, which is done using the `solute_list` configuration argument (as shown above). The solutes defined in the `solute_list` configuration argument are automatically added to the property package component list and concentration terms will be created for each of these.

Class Documentation
-------------------

* :class:`ZOParameterBlock`
* :class:`ZOParameterBlockData`
* :class:`ZOStateBlock`
* :class:`ZOStateBlockData`