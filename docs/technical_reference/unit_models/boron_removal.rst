Boron Removal
=============

Introduction
------------

Boron removal is generally applied as a desalination post-treatment after the first
Reverse Osmosis (RO) stage to shift the equilibrium of Boron (B[OH]\ :sub:`3`\) to Borate (B[OH]\ :sub:`4`\ :sup:`-`)
so that a second RO stage can remove Borate ions from solution such that the total boron
concentration can go below ~1 mg/L per the World Health Organization standards. This
is accomplished by adding caustic chemicals (such as sodium hydroxide (NaOH) or lime
(Ca[OH]\ :sub:`2`\)) to the stream just before the 2nd RO stage. The goal is to shift
the pH of the stream above around 9.5 so that the majority of total boron is in the
ionic state of Borate. Those Borate ions would then be rejected in the 2nd RO stage.

In this implementation of the model, the user MUST provide a state variable for the
molar flow of Boron and Borate (at a minimum). Optionally, users may also have molar
flows for protons, hydroxide anions, and the cation species associated with the caustic
additive used in the unit. This will help track added salts and changes in pH when
the caustic chemical is added. The unit will always provide an approximation to the
mixture pH whether or not the user has state variables for protons in their property
package. Users must also specify a dosage of the caustic chemical added to the reactor
as well as the molecular weight of the chemical additive (see Chemical Dosing Parameters
section below for more detailed information).

The main assumptions of the implemented model are as follows:

1) Model dimensionality is limited to a 0D control volume
2) Single liquid phase only
3) Isothermal operation
4) Only major reactions occurring within this process are Boron/Borate dissociation and water dissociation
5) No other residual ion species change significantly when the caustic chemical is added

.. figure:: ../../_static/unit_models/BoronRemovalDiagram.png
    :width: 600
    :align: center

    Figure 1. Schematic representation of a Boron removal unit modeled in WaterTAP

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
   "Components", ":math:`j`", "['H2O', 'Boron', 'Borate', 'Protons', 'Hydroxide', 'Caustic_Cation', ...]"

**Users can name the 'Components' however they want, but they must also provide the unit model with a mapping of component names to an associated state variable name used in the property package**

**State variables must have 'H2O', 'Boron', and 'Borate'. All other components are optional**

Degrees of Freedom and Variables
--------------------------------
Aside from the inlet feed state variables (i.e., temperature, pressure, component molar flowrates),

the boron removal model has 1 additional degree of freedom that
the user must specify.

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Fluid temperature", ":math:`T`", "temperature", "[t]", ":math:`\text{K}`"
   "Fluid pressure", ":math:`P`", "pressure", "[t]", ":math:`\text{Pa}`"
   "Molar flowrate of components", ":math:`M_j`", "flow_mass_phase_comp", "[t, 'Liq', j]", ":math:`\text{mol/s}`"
   "Caustic Chemical Dose", ":math:`D`", "caustic_dose", "[t]", ":math:`\text{mg/L}`"

**Users must provide values for and 'fix' these variables to solve the model with DOF=0. However, users may also leave variables unfixed for optimization purposes.**

**NOTE: Variables for 'temperature', 'pressure', and 'flow_mol_phase_comp' come from the associated property package as state variables and are accessed via {port_name}.{state_var_name}**


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
            default={
                "property_package": model.fs.properties,
                "chemical_additives": chem_dict,
            }
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
M.M. Benjamin, `Water Chemistry <https://www.biblio.com/9781577666677>`_, Waveland Press,
Inc.: Illinois, 2010, Ch. 1, 18-51.

Lenntech, `Desalination Post-treatment: Boron Removal Process <https://www.lenntech.com/processes/desalination/post-treatment/post-treatments/boron-removal.htm>`_,
Accessed May 16, 2022.
