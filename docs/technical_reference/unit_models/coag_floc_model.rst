.. _coagulation_flocculation:

Coagulation-Flocculation
========================

.. code-block:: python

   from watertap.unit_models.coag_floc_model import CoagulationFlocculation

Introduction
------------

Coagulation-Flocculation is a water treatment process  designed to remove total suspended
solids (TSS) from solution by converting small suspended and/or colloidal matter, and some
natural organic matter, into larger 'floccs' (or 'Sludge') that can then be separated by
sedimentation and/or filtration (or other separation processes) later in the treatment train.
This is accomplished in a 2-stage process:

Stage 1) Coagulation (where chemical coagulants and other aids are rapidly mixed into solution)
Stage 2) Flocculation (where the 'floccs' are formed through slow, gentle mixing/agitation)

In this implementation of the model, the user MUST provide a measured final Turbidity (in NTU) made
during a Jar Test for the given water source. This measurement is then used to estimate how much
TSS would be removed during the Coagulation-Flocculation process. Users may also
provide the specific chemical composition of additives used to achieve this final Turbidity, and
can provide information on the level of salts in those additives. This information can be used
to estimate an increase in total dissolved salts (TDS) that may occur due to the addition of
those chemicals.

This model also includes relationships for power usage in the rapid mixing and flocculation
basins. Users will need to provide information for power usage such as retention times of
each basin, mixing paddle sizes, number of mixers, etc.

The main assumptions of the implemented model are as follows:

1) Coagulation-Flocculation can be modeled together in a single combined unit (Figure 1)
2) Model dimensionality is limited to a 0D control volume
3) Predicted levels of suspended solid removal can be determined solely by Jar Test measurements
4) Single liquid phase only
5) Isothermal operation

.. figure:: ../../_static/unit_models/coagulation_flocculation.png
    :width: 600
    :align: center

    Figure 1. Schematic representation of a coagulation-flocculation unit modeled in WaterTAP

Ports
-----

The model provides two ports (Pyomo notation in parenthesis):

* Inlet port (inlet)
* Outlet port (outlet)

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq']"
   "Components", ":math:`j`", "['H2O', 'TDS', 'TSS', 'Sludge', ...]"
   "Chemical Additives", ":math:`i`", "['chem_A', 'chem_B', ...]"

**Users are responsible for naming any chemical additives and defining all parameters associated with them**

Degrees of Freedom and Variables
--------------------------------
Aside from the inlet feed state variables (i.e., temperature, pressure, component mass flowrates),

the Coagulation-Flocculation model has at least an additional 13 degrees of freedom that
the user must specify. The table below gives an outline of these.

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Fluid temperature", ":math:`T`", "temperature", "[t]", ":math:`\text{K}`"
   "Fluid pressure", ":math:`P`", "pressure", "[t]", ":math:`\text{Pa}`"
   "Mass flowrate of components", ":math:`M_j`", "flow_mass_phase_comp", "[t, 'Liq', j]", ":math:`\text{kg/s}`"
   "Slope relationship between measured Turbidity and TSS", ":math:`a`", "slope", "[t]", ":math:`\text{mg/L/NTU}`"
   "Intercept relationship between measured Turbidity and TSS", ":math:`b`", "intercept", "[t]", ":math:`\text{mg/L}`"
   "Turbidity measured before Jar Test", ":math:`Turb_o`", "initial_turbidity_ntu", "[t]", ":math:`\text{NTU}`"
   "Turbidity measured after Jar Test", ":math:`Turb_f`", "final_turbidity_ntu", "[t]", ":math:`\text{NTU}`"
   "Chemical Doses added during Jar Test", ":math:`D_i`", "chemical_doses", "[t, i]", ":math:`\text{mg/L}`"
   "Retention time for each rapid mixer", ":math:`\tau_r`", "rapid_mixing_retention_time", "[t]", ":math:`\text{s}`"
   "Number of rapid mixers in series", ":math:`n_r`", "num_rapid_mixing_basins", "None", "None"
   "Rapid mixer velocity gradient", ":math:`G_r`", "rapid_mixing_vel_grad", "[t]", ":math:`\text{s}^{-1}`"
   "Retention time of flocculation basin", ":math:`\tau_f`", "floc_retention_time", "[t]", ":math:`\text{s}`"
   "Flocculation single paddle length (from center of rotation to blade edge)", ":math:`L`", "single_paddle_length", "None", ":math:`\text{m}`"
   "Flocculation single paddle width", ":math:`w`", "single_paddle_width", "None", ":math:`\text{m}`"
   "Flocculation paddle rotational speed", ":math:`\omega`", "paddle_rotational_speed", "[t]", ":math:`\text{revolutions/s}`"
   "Flocculation paddle drag coefficient", ":math:`C_D`", "paddle_drag_coef", "[t]", "None"
   "Flocculation paddle velocity fraction", ":math:`f`", "vel_fraction", "None", "None"
   "Number of rotating paddle wheels", ":math:`n_w`", "num_paddle_wheels", "None", "None"
   "Number of paddles per wheel", ":math:`n_p`", "num_paddles_per_wheel", "None", "None"

**Users must provide values for and 'fix' these variables to solve the model with DOF=0. However, users may also leave variables unfixed for optimization purposes.**


**NOTE: Default values are provided for the slope and intercept relationships between Turbidity and TSS. These come from Rugner et al. (2013) but can be substituted as needed to match any data available relating turbidity to TSS.**

**NOTE: Variables for 'temperature', 'pressure', and 'flow_mass_phase_comp' come from the associated property package as state variables and are accessed via {port_name}.{state_var_name}**


Chemical Dosing Parameters
--------------------------
In addition to providing and fixing values for chemical additives, the users will
need to provide parameter information for each additive including molecular weight (:math:`MW_{a,i}`),
moles of salt that would be added per mole of additive (:math:`N_i`), and a representative molecular
weight of the salt species that would be formed from addition of the additive (:math:`MW_{s,i}`). If a user
does not have this information off hand, then the user can simply give a value of '0' for
the moles of salt added per mole of additive (and dummy values for the molecular weights).
This information is only used to estimate the rise in TDS when salts are added, so it
is not critical for the determination of the main objective of Coagulation-Flocculation,
which is the removal of TSS.

To provide this information to the unit model, users must add a 'chemical_additives'
dictionary to the initialization of the unit model. That dictionary must have the
following format.

.. code-block::

   chem_dict = {'chem_A':
                  {'parameter_data':
                    {'mw_additive': (value, units),
                     'moles_per_mole_additive': value,
                     'mw_salt': (value, units)
                    }
                  },
                'chem_B':
                  {'parameter_data':
                    {'mw_additive': (value, units),
                     'moles_per_mole_additive': value,
                     'mw_salt': (value, units)
                    }
                  }
              }

For example, this 'chem_dict' would be passed into the model on construction as
one of the configuration options as shown below.

.. code-block::

    model.fs.unit = CoagulationFlocculation(
        property_package=model.fs.properties,
        chemical_additives=chem_dict,
    )

**NOTE: The above example assumes you have already constructed a pyomo model named 'model' and attached an IDAES flowsheet named 'fs' to it, as well as a properties block named 'properties'**

Equations and Relationships
---------------------------

.. csv-table::
   :header: "Description", "Equation"

   "TSS relationship with initial Turbidity", ":math:`TSS_o = b + a(Turb_o)`"
   "TSS relationship with final Turbidity", ":math:`TSS_f = b + a(Turb_f)`"
   "TSS loss rate", ":math:`S_{TSS} = M_{TSS,in} - Q \cdotp TSS_f`"
   "TSS mass balance", ":math:`0 = M_{TSS,in} - M_{TSS,out} - S_{TSS}`"
   "Sludge mass balance", ":math:`0 = M_{Sludge,in} - M_{Sludge,out} + S_{TSS}`"
   "TDS gain rate", ":math:`S_{TDS} = Q \cdotp {\sum_{i} \frac{D_i}{MW_{a,i}} \cdotp N_i \cdotp MW_{s,i} }`"
   "TDS mass balance", ":math:`0 = M_{TDS,in} - M_{TDS,out} + S_{TDS}`"
   "Rapid Mixer Total Volume", ":math:`V_r = Q \cdotp \tau_r \cdotp n_r`"
   "Rapid Mixer Total Power Usage", ":math:`P_r = {G_r}^2 \cdotp \mu \cdotp V_r`"
   "Flocculation Basin Total Volume", ":math:`V_f = Q \cdotp \tau_f`"
   "Paddle Wheel Speed", ":math:`v_p = \pi \cdotp L \cdotp \omega`"
   "Flocculation Power Usage", ":math:`P_p = 0.5 \cdotp C_D \cdotp L \cdotp w \cdotp n_w \cdotp n_p \cdotp \rho {(f \cdotp v_p)}^3`"
   "Total Power Usage", ":math:`P_T = P_p + P_r`"

**Relationships for power usage all come from Mines (2014)**

**NOTE:** :math:`Q` **is defined as the total volumetric flow rate and** :math:`S_{j}` **is the source/sink term for component** :math:`j`

References
----------
H. Rugner, M. Schwientek, B. Beckingham, B. Kuch, P. Grathwohl,
Environ. Earth Sci. 69 (2013) 373-380. DOI:`10.1007/s12665-013-2307-1 <https://doi.org/10.1007/s12665-013-2307-1>`_

R.O. Mines `Environmental Engineering: Principles and Practice <https://www.biblio.com/9781118801451>`_,
1st Ed, John Wiley & Sons, 2014. Ch. 6.
