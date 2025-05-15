.. _MD_1D:

Membrane Distillation (1D)
=========================================

.. code-block:: python

   from watertap.unit_models.MD.membrane_distillation_1D import MembraneDistillation1D

This Membrane Distillation (MD) unit model:
   * supports the following configurations: 

     - **DCMD** (Direct Contact Membrane Distillation)
     - **VMD** (Vacuum Membrane Distillation)
     - **GMD** (Permeate Gap/Conductive Gap Membrane Distillation)

   * is 1-dimensional
   * supports steady-state only
   * assumes heat loss in equipment is negligible
   * assumes permeate exits the membrane pores with zero salinity
   * assumes no concentration polarization for the cold channel
   * assumes complete vapor condensation for the cold channel (in DCMD and GMD)
   * accounts for vapor expansion in VMD
   * assumes linear temperature change across gap channel (in GMD)
   * assumes no pressure change and temperature polarization in VMD vaccuum channel


Degrees of Freedom
------------------
In addition to the hot channel and cold channel inlet state variables (i.e, temperature, pressure, and component flowrates) for the **DCMD** and **GMD** configurations, the MD model has at least **4 degrees of freedom** for all configurations that should be fixed for the unit to be fully specified. Typically, the following variables are fixed:

- Membrane permeability coefficient
- Membrane thickness
- Membrane thermal conductivity
- Recovery *or* membrane area

**Additional degress of freedom**:

- **VMD** introduces vacuum pressure at the cold side.
- **GMD** introduces gap thermal conductivity and gap thickness.

Configuring the MD unit to calculate temperature polarization, concentration polarization, mass transfer
coefficient, and pressure drop would result in five additional degrees of freedom. In this case, in addition to the
previously fixed variables, we typically fix the following variables to fully specify the unit:

    * Hot channel spacer porosity
    * Hot channel height
    * Cold channel spacer porosity (in DCMD and GMD)
    * Cold channel height (in DCMD and GMD)
    * Membrane length *or* membrane width

Model Structure
---------------
The MD model consists of a separate `MDchannel1Dblock` for each channel depending on the configuration:

- **DCMD**: Includes **hot channel** and **cold channel**.
- **VMD**: Includes **hot channel** and **vacuum (cold) channel**.
- **GMD**: Includes **hot channel**, **gap channel**, and **cold channel**.

- **hot and cold channels in all configurations** includes bulk properties StateBlocks indexed by time and space which are used for mass, energy, and momentum balances
- **hot channel in all configurations, cold channel in DCMD and GMD, and gap channel in GMD** includes StateBlocks indexed by time and space for the conditions at the membrane interface and gap interface
- **hot channel in all configurations, cold channel in DCMD, and gap channel in GMD** includes StateBlocks indexed by time and space for Vapor properties at the membrane interface (for **DCMD** and **VMD** configurations).

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Space", ":math:`x`", "None"
   "Phases", ":math:`p`", "['Liq', 'Vap']"
   "Components", ":math:`j`", "['H2O', solute]*"

\*Solute depends on the imported property model.

Variables
----------

Refer to the :any:`0MD_variables` section in the  0D MD model.

.. _1MD_equations:

Equations
-----------

Refer to the :any:`0MD_equations` section in the  0D MD model.

Class Documentation
-------------------

* :mod:`watertap.unit_models.MD.membrane_distillation_1D`
* :mod:`watertap.unit_models.MD.membrane_distillation_base`
* :mod:`watertap.unit_models.MD.MD_channel_1D`
* :mod:`watertap.unit_models.MD.MD_channel_base`
