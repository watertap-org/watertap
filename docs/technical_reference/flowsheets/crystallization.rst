NaCl Crystallization
===============================

Introduction
------------

Crystallization is the precipitation and extraction of the crystals from a mother liquor based on solute supersaturation. Crystallization is a potential solution to desalination's environmental brine management problem; brine crystallization systems shown to achieve zero-liquid discharge and salt recovery.

The objective of the mathematical model is to simulate crystallization from a NaCl mother liquor. The example is based on (but does not replicate exactly) Example 11.4 in `Tavare, N.S. (1995). <https://link.springer.com/chapter/10.1007/978-1-4899-0233-7_11>`_ 

Implementation
--------------

The modeled crystallization proecess is illustrated by Figure 1. Thermal energy for the evaporative crystallizer is added in the suspension heater, with the water vapor product sent off to a condenser (not modeled here). 

The flowsheet relies on the following key assumptions:

   * supports steady-state only.
   * a property package supporting all three phases (e.g., the crystallizer property package) is provided for the crystallizer model.

.. figure:: ../../_static/flowsheets/crystallizer.jpg
    :width: 500
    :align: center

    Figure 1. Crystallizer flowsheet

Documentation for the crystallizer unit model can be found below. 
    * `Crystallizer <https://watertap.readthedocs.io/en/latest/technical_reference/unit_models/crystallizer_0D.html>`_

Degrees of Freedom
------------------
For the crystallizer unit model, if the inlet feed condition is fully specified, the user is left with two degrees of freedom. In this example, the following variables are initially specified for simulating the crystallizer (i.e., degrees of freedom = 0):

    * NaCl mother liqour conditions (i.e., flow, temperature, pressure, compositions)
    * crystallizer operating temperature
    * Solid outlet mass flow

However, it should be noted that any of the following variables could have been fixed as alternatives to the solid outlet mass flow:
    * Crystallization yield
    * Solids volumetric fraction of product slurry
    * Magma density
    * Vapour outlet mass flow
    * Crystallizer thermal energy input

We demonstrate how this may be done in our flowsheet by updating and re-solving the model. Thus, the example shows five cases:
#. **Case 1**: fixing the crystallizer operating temperature and solid outlet mass flow
#. **Case 2**: fixing the crystallizer operating temperature and crystallization yield
#. **Case 3**: fixing the crystallizer operating temperature and the product slurry solids volumetric fraction
#. **Case 4**: fixing the crystallizer operating temperature and magma density
#. **Case 5**: fixing the crystallizer operating temperature and crystallizer thermal energy input

Flowsheet Specifications
------------------------

.. csv-table::
   :header: "Description", "Value", "Units"

   "**NaCl mother liquor**:math:`^1`"
   "NaCl feed mass flow (liquid phase)","10.5119", ":math:`\text{kg}\text{/s}`"
   "Water feed mass flow (liquid phase)","38.9326", ":math:`\text{kg}\text{/s}`"
   "NaCl feed mass flow (solid phase)","0", ":math:`\text{kg}\text{/s}`"
   "Water feed mass flow (vapor phase)","0", ":math:`\text{kg}\text{/s}`"
   "Feed temperature", "293.15", ":math:`\text{K}`"
   "Pressure", "101325", ":math:`\text{Pa}`"

   "**Crystallizer**:math:`^1`"
   "Crystallizer operating temperature", "348.15", ":math:`\text{K}`"

   "**Other degree of freedom**"
   "Solid outlet mass flow :math:`^1`", "5.556", ":math:`\text{kg/}\text{s}`"
   "Crystallization yield", "70", ":math:`\text{%}`"
   "Solids volumetric fraction of product slurry", "0.1182", ":math:`\text{dimensionless}`"
   "Magma density :math:`^1`", "250", ":math:`\text{kg/}\text{m}^3`"
   "Crystallizer thermal energy input :math:`^1`", "55,000", ":math:`\text{kW}`"



References
----------
[1] Tavare, N.S. (1995). Crystallizer Design and Operation. In: Industrial Crystallization. The Springer Chemical Engineering Series. Springer, Boston, MA. https://doi.org/10.1007/978-1-4899-0233-7_11
