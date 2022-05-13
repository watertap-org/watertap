Electrodialysis (0D)
====================

Introduction
------------

Electrodialysis, an electrochemical seperation technology for water desalination, has been 
utilized for decades.  Compared with other technologies such as reverse osmosis (RO),
electridialysis showes advtanges in desalinating brackish waters with 
moderate salinity by its less intense energy consumption, higher recovery, and robust 
tolerance for adverse non-ionic components (e.g., silica and biological substances) in water.
Between the electrodes of a eletrodialysis stack (a reactor module), multiple anion- and 
cation-exchange membranes are alternatively positioned and seprated by spacers to form individual 
cells. When operated, electrolidalysis converts electrical current to ion flux and, assisted by
the opposite ion-exchage membranes, moves ion from one cell (diluate channel) to its adjencent cell
(concentrate channel) in a cell-pair treatment unit (Figure 1). Recovered (desalianted) water is 
collected from diluate channles of all cell pairs while the concentrate product can be disposed as brine 
or retreated. 

.. figure:: ../../_static/unit_models/xxx.png
    :width: 600
    :align: center

    Figure 1. Schematic representation of a electrodialysis cell pair as a model unit


One cell pair in a electrodialysis stack can thus be treated as a modeling unit that can multiply to
larger scale systems.  This model establishes mathetmaical descriptions of ion and water transport in
a cell pair and expands to simulate a stack with a specified cell-pair number.  Modeled mass transfer 
mechanisms inlcude electrical migration and diffusion of ions and osmosis and electroosmosis of water. 
The following key assumptions are based on. 

* The concentrate and diluate channels ahve identical geometry and fluidic conditions. 
* For each channel, component fluxes are uniform in the bulk solutions (0-dimensional assumption) 
    and are set as the average of inlete and outlect of each channel. 
* Co-current flow operation. 
* Electrical current is operated below the limiting current. 
* Ideality assumptions: activity, osmotic, and van't Hoff coefficients are set at one. 
* All ion-exchange membrane properties (ion and water transport number, resistence, permeabiliity) are
    constant. 
* Detailed concentration gradient effect at membrane-water interfaces is neglected. 
* Constant pressure and temperature through each channel. 

Control Volumes
-----

This model has two control volumes for the concentrate and diluate channels.

* diluate_channel
* concentrate_channel 

Ports
-----

This model provides four ports (Pyomo notation in parenthesis) on the two :

* inlet_diluate (inlet)
* outlet_diluate (outlet)
* inlet_concentrate (inlet)
* outlet_concentrate (outlet)

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq']"
   "Components", ":math:`j`", "['H2O', 'Na_+', 'Cl_-', ...]"
   "Membrane", :math: `na`, ['cem', ;'aem']

**Users are responsible for naming any chemical additives and defining all parameters associated with them**

Degrees of Freedom
------------------
Aside from the inlet feed state variables (i.e. temperature, pressure, component flowrates), this model has
31 additional degerees of freedom that need to be fixed to fully solve the model. 15 of them are repeated over 
the two channels (making up 30) and 1 (the electrical input) is general for both channels.  The total 31 are 
tabulated below. 


.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units", "Number"

   "Water transport number", ":math:`t_w`", "water_trans_number_membrane", "["cem", "aem"], ":math:`\text{dimensionless}`", 2
   "Water permeability", ":math:`L`", "water_permeability_membrane", "["cem", "aem"]", ":math:`\ `", 

   "Voltage", ":math:`U`", "voltage", "[t]", ":math:``"
   "Resistence of the electrode compartments", ":math:`R_{el}`", "electrodes_resistence", "[t]", ":math:`\text{mg/L}`"
   "Cell pair number", ":math:`N`", "cell_pair_numb", "None", ":math:``"
   "Current utilization for all mass transfer", ":math:`\xi`", "final_turbidity_ntu", "None", ":math:``"
   "Spacer thickness", ":math:`s`", "spacer_thickness", "none", ":math:``"
   "Area resistence of ion exchnge membranes", ":math:`r`", "membrane_surface_resistence", "["cem", "aem"]", ":math:`\text{s}`"
   "Cell width", ":math:`b`", "cell_width", "None", ":math: `\text{m}`"
   "Cell length", ":math:`l`", "cell_length", "None", ":math: `\text{m}`"
   "Thickness of ion exchange membranes", ":math:`\delta`", "memebrane_thickness", "["cem", "aem"]", ":math:`\text{m}`"
   "diffusivity of solute in the membrane phase", ":math:`\D`", "solute_diffusivity_membrane", "[["cem", "aem"], j]", ":math:`\text{m^2/s}`"
   "transport number of ions in the membrane phase", ":math:`t`", "ion_trans_number_membrane", "[["cem", "aem"], j]", ":math:`\text{dimensionless}`"
   
at least 4 degrees of freedom that should be fixed for the unit to be fully specified.

Typically, the following variables are fixed, in addition to state variables at the inlet:
    * membrane water permeability, A
    * membrane salt permeability, B
    * permeate pressure
    * membrane area

On the other hand, configuring the RO unit to calculate concentration polarization effects, mass transfer
coefficient, and pressure drop would result in 3 additional degrees of freedom. In this case, in addition to the
previously fixed variables, we typically fix the following variables to fully specify the unit:

    * feed-spacer porosity
    * feed-channel height
    * membrane length *or* membrane width *or* inlet Reynolds number

Model Structure
------------------
This RO model consists of 2 ControlVolume0DBlocks: one for the feed-side and one for the permeate-side.
 
* The feed-side includes 2 StateBlocks (properties_in and properties_out) which are used for mass, energy, and momentum balances, and 2 additional StateBlocks for the conditions at the membrane interface (properties_interface_in and properties_interface_out).
* The permeate-side includes 3 StateBlocks (properties_in, properties_out, and properties_mixed). The inlet and outlet StateBlocks are used to only determine the permeate solute concentration for solvent and solute flux at the feed-side inlet and outlet, while the mixed StateBlock is used for mass balance based on the average flux.

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Inlet/outlet", ":math:`x`", "['in', 'out']"
   "Phases", ":math:`p`", "['Liq']"
   "Components", ":math:`j`", "['H2O', 'NaCl']*"

\*Solute depends on the imported property model; example shown here is for the NaCl property model.

.. _0dro_variables:

Variables
----------

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Solvent permeability coefficient", ":math:`A`", "A_comp", "[t, j]", ":math:`\text{m/Pa/s}`"
   "Solute permeability coefficient", ":math:`B`", "B_comp", "[t, j]", ":math:`\text{m/s}`"
   "Mass density of solvent", ":math:`\rho_{solvent}`", "dens_solvent", "[p]", ":math:`\text{kg/}\text{m}^3`"
   "Mass flux across membrane", ":math:`J`", "flux_mass_phase_comp", "[t, x, p, j]", ":math:`\text{kg/s}\text{/m}^2`"
   "Membrane area", ":math:`A_m`", "area", "None", ":math:`\text{m}^2`"
   "Component recovery rate", ":math:`R_j`", "recovery_mass_phase_comp", "[t, p, j]", ":math:`\text{dimensionless}`"
   "Volumetric recovery rate", ":math:`R_{vol}`", "recovery_vol_phase", "[t, p]", ":math:`\text{dimensionless}`"
   "Observed solute rejection", ":math:`r_j`", "rejection_phase_comp", "[t, p, j]", ":math:`\text{dimensionless}`"
   "Over-pressure ratio", ":math:`P_{f,out}/Δ\pi_{out}`", "over_pressure_ratio", "[t]", ":math:`\text{dimensionless}`"
   "Mass transfer to permeate", ":math:`M_p`", "mass_transfer_phase_comp", "[t, p, j]", ":math:`\text{kg/s}`"

The following variables are only built when specific configuration key-value pairs are selected.

if ``has_pressure_change`` is set to ``True``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Pressure drop", ":math:`ΔP`", "deltaP", "[t]", ":math:`\text{Pa}`"

if ``concentration_polarization_type`` is set to ``ConcentrationPolarizationType.fixed``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Concentration polarization modulus", ":math:`CP_{mod}`", "cp_modulus", "[t, j]", ":math:`\text{dimensionless}`"

if ``concentration_polarization_type`` is set to ``ConcentrationPolarizationType.calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Mass transfer coefficient in feed channel", ":math:`k_f`", "Kf", "[t, x, j]", ":math:`\text{m/s}`"

if ``mass_transfer_coefficient`` is set to ``MassTransferCoefficient.calculated``
or ``pressure_change_type`` is set to ``PressureChangeType.calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Feed-channel height", ":math:`h_{ch}`", "channel_height", "None", ":math:`\text{m}`"
   "Hydraulic diameter", ":math:`d_h`", "dh", "None", ":math:`\text{m}`"
   "Spacer porosity", ":math:`\epsilon_{sp}`", "spacer_porosity", "None", ":math:`\text{dimensionless}`"
   "Reynolds number", ":math:`Re`", "N_Re", "[t, x]", ":math:`\text{dimensionless}`"


if ``mass_transfer_coefficient`` is set to ``MassTransferCoefficient.calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Schmidt number", ":math:`Sc`", "N_Sc", "[t, x]", ":math:`\text{dimensionless}`"
   "Sherwood number", ":math:`Sh`", "N_Sh", "[t, x]", ":math:`\text{dimensionless}`"

if ``mass_transfer_coefficient`` is set to ``MassTransferCoefficient.calculated``
or ``pressure_change_type`` is **NOT** set to ``PressureChangeType.fixed_per_stage``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Membrane length", ":math:`L`", "length", "None", ":math:`\text{m}`"
   "Membrane width", ":math:`W`", "width", "None", ":math:`\text{m}`"

if ``pressure_change_type`` is set to ``PressureChangeType.fixed_per_unit_length``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Average pressure drop per unit length of feed channel", ":math:`(\frac{ΔP}{Δx})_{avg}`", "dP_dx", "[t]", ":math:`\text{Pa/m}`"

if ``pressure_change_type`` is set to ``PressureChangeType.calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Feed-channel velocity", ":math:`v_f`", "velocity", "[t, x]", ":math:`\text{m/s}`"
   "Friction factor", ":math:`f`", "friction_factor_darcy", "[t, x]", ":math:`\text{dimensionless}`"
   "Pressure drop per unit length of feed channel at inlet/outlet", ":math:`ΔP/Δx`", "dP_dx", "[t, x]", ":math:`\text{Pa/m}`"

.. _0dro_equations:

Equations
-----------

.. csv-table::
   :header: "Description", "Equation"

   "Solvent flux across membrane", ":math:`J_{solvent} = \rho_{solvent} A(P_{f} - P_p - (\pi_{f}-\pi_{p}))`"
   "Solute flux across membrane", ":math:`J_{solute} = B(C_{f} - C_{p})`"
   "Average flux across membrane", ":math:`J_{avg, j} = \frac{1}{2}\sum_{x} J_{x, j}`"
   "Permeate mass flow by component j", ":math:`M_{p, j} = A_m J_{avg,j}`"
   "Permeate-side solute mass fraction", ":math:`X_{x, j} = \frac{J_{x, j}}{\sum_{x} J_{x, j}}`"
   "Feed-side membrane-interface solute concentration", ":math:`C_{interface} = CP_{mod}C_{bulk}=C_{bulk}\exp(\frac{J_{solvent}}{k_f})-\frac{J_{solute}}{J_{solvent}}(\exp(\frac{J_{solvent}}{k_f})-1)`"
   "Concentration polarization modulus",":math:`CP_{mod} = C_{interface}/C_{bulk}`"
   "Mass transfer coefficient",":math:`k_f = \frac{D Sh}{d_h}`"
   "Sherwood number",":math:`Sh = 0.46 (Re Sc)^{0.36}`"
   "Schmidt number",":math:`Sc = \frac{\mu}{\rho D}`"
   "Reynolds number",":math:`Re = \frac{\rho v_f d_h}{\mu}`"
   "Hydraulic diameter",":math:`d_h = \frac{4\epsilon_{sp}}{2/h_{ch} + (1-\epsilon_{sp})8/h_{ch}}`"
   "Cross-sectional area",":math:`A_c = h_{ch}W\epsilon_{sp}`"
   "Membrane area",":math:`A_m = LW`"
   "Pressure drop",":math:`ΔP = (\frac{ΔP}{Δx})_{avg}L`"
   "Feed-channel velocity",":math:`v_f = Q_f/A_c`"
   "Friction factor",":math:`f = 0.42+\frac{189.3}{Re}`"
   "Pressure drop per unit length",":math:`\frac{ΔP}{Δx} = \frac{1}{2d_h}f\rho v_f^{2}`"
   "Component recovery rate",":math:`R_j = \frac{M_{p,j}}{M_{f,in,j}}`"
   "Volumetric recovery rate",":math:`R_{vol} = \frac{Q_{p}}{Q_{f,in}}`"
   "Observed solute rejection", ":math:`r_j = 1 - \frac{C_{p,mix}}{C_{f,in}}`" 

Class Documentation
-------------------

.. automodule:: watertap.unit_models.reverse_osmosis_0D
    :members:
    :noindex:

