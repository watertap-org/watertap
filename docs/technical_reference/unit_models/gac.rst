Granular Activated Carbon (GAC)
====================
.. note::

    The GAC model is an ongoing project and subject to change, as is this documentation.

This is a granular activated carbon (GAC) model that works under the following criteria and assumptions
   * supports a single liquid phase only
   * supports a single solute and single solvent (water) only
   * supports steady-state only
   * assumes isothermal conditions
   * model performance is independent of a gravity-fed or pressurized GAC unit, therefore assumes isobaric

.. index::
   pair: watertap.unit_models.gac;gac

.. currentmodule:: watertap.unit_models.gac

Degrees of Freedom
------------------
In the default configuration of the GAC unit model there are 18 degrees of freedom in addition to the inlet state variables
(i.e. temperature, pressure, component flowrates) that should be fixed for the model to be fully specified.
In association with using the Freundlich adsorption isotherm and 5-parameter model, the following 9 variables are almost always fixed
and may be derived from experimental data:

   * Freundlich isotherm k parameter
   * Freundlich isotherm 1/n parameter
   * Stanton number equation parameters :math:`\text{a}_0` and :math:`\text{a}_1`
   * throughput ratio equation parameters :math:`\text{b}_0`, :math:`\text{b}_1`, :math:`\text{b}_2`, :math:`\text{b}_3` and :math:`\text{b}_4`

Additionally, the following 9 variables are traditionally fixed:

   * the target dimensionless concentration or the dimensionless concentration at the time of replacement
   * empty bed contact time
   * bed voidage
   * superficial velocity
   * GAC particle porosity
   * GAC particle apparent density
   * GAC particle diameter
   * liquid phase film transfer coefficient
   * surface diffusion coefficient

When setting the configuration options to calculate the liquid phase film transfer coefficient and surface diffusion coefficient,
these respective variables are no longer specified and 4 newly introduced variables must be fixed. This is a net result of 20 degrees of freedom.
Newly utilized variables that must be fixed include:

   * molal volume of the solute
   * GAC particle sphericity
   * tortuosity of the path that the adsorbate must take as compared to the radius
   * surface-to-pore diffusion flux ratio

.. note::

    Pending explanation of figure

.. figure:: ../../_static/unit_models/gac.png
    :width: 1200
    :align: center
    
    Figure 1.

Model Structure
------------------
The GAC model consists of 1 ControlVolume0DBlock (process_flow) for the process flow of the water treatment train.
The process flow includes 2 StateBlocks (inlet and outlet) which are used for mass and momentum balances.
It also includes 1 StateBlock (adsorbed) for the solute that is adsorbed into the GAC parties.

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq']"
   "Components", ":math:`j`", "['H2O', 'Solute']*"

\*Solute is provided in the imported property model.

.. _GAC_variables:

Variables
----------

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Freundlich isotherm k parameter", ":math:`k`", "freund_k", "None", ":math:`\left(\text{m}^3\text{/kg}\right)^\left( \frac{1}{n} \right)`"
   "Freundlich isotherm 1/n parameter", ":math:`\frac{1}{n}`", "freund_ninv", "None", ":math:`\text{dimensionless}`"
   "Equilibrium concentration of adsorbed phase with liquid phase", ":math:`q_e`", "equil_conc", "None", ":math:`\left( \text{kg}_\text{adsorbate}\text{/kg}_\text{adsorbent} \right)`"
   "Mass of contaminant absorbed at the time of replacement", ":math:`M_{solute}`", "mass_adsorbed", "None", ":math:`\text{kg}`"
   "Mass of contaminant adsorbed if fully saturated", ":math:`M_\text{solute\text{,}e}`", "mass_adsorbed_saturated", "None", ":math:`\text{kg}`"
   "Adsorber bed void fraction", ":math:`\epsilon`", "bed_voidage", "None", ":math:`\text{dimensionless}`"
   "Adsorber bed volume", ":math:`V`", "bed_volume", "None", ":math:`\text{m}^3`"
   "Adsorber bed area", ":math:`A`", "bed_area", "None", ":math:`\text{m}^2`"
   "Adsorber bed length", ":math:`L`", "bed_length", "None", ":math:`\text{m}`"
   "Mass of fresh GAC in the bed", ":math:`M`", "bed_mass_gac", "None", ":math:`\text{kg}`"
   "Superficial velocity", ":math:`v_s`", "velocity_sup", "None", ":math:`\text{m/s}`"
   "Interstitial velocity", ":math:`v_i`", "velocity_int", "None", ":math:`\text{m/s}`"
   "GAC particle porosity", ":math:`\epsilon_p`", "particle_porosity", "None", ":math:`\text{dimensionless}`"
   "GAC particle apparent density", ":math:`\rho_a`", "particle_dens_app", "None", ":math:`\text{kg/}\text{m}^3`"
   "GAC particle bulk density", ":math:`\rho_b`", "particle_dens_bulk", "None", ":math:`\text{kg/}\text{m}^3`"
   "GAC particle solid density", ":math:`\rho_s`", "particle_dens_sol", "None", ":math:`\text{kg/}\text{m}^3`"
   "GAC particle diameter", ":math:`d_p`", "particle_dia", "None", ":math:`\text{m}`"
   "Average dimensionless concentration in the duration of operation", ":math:`\bar{\left( \frac{C}{C_{0}} \right)}`", "conc_ratio_avg", "None", ":math:`\text{dimensionless}`"
   "Dimensionless concentration of the effluent at the time of replacement", ":math:`\left( \frac{C}{C_{0}} \right)|_{z=L}`", "conc_ratio_avg", "None", ":math:`\text{dimensionless}`"
   "Approximate saturation of the GAC in the bed at the time of replacement", ":math:`\left( \frac{q}{q_{e}} \right)`", "gac_saturation_replace", "None", ":math:`\text{dimensionless}`"
   "Empty bed contact time", ":math:`EBCT`", "ebct", "None", ":math:`\text{s}`"
   "Mass throughput ratio", ":math:`T`", "mass_throughput", "None", ":math:`\text{dimensionless}`"
   "Residence time", ":math:`\tau`", "res_time", "None", ":math:`\text{s}`"
   "Elapsed operation time", ":math:`t`", "elap_time", "None", ":math:`\text{s}`"
   "Steady state GAC replacement rate", ":math:`\dot{m}_{GAC}`", "gac_mass_replace_rate", "None", ":math:`\text{s}`"
   "Liquid phase film transfer coefficient", ":math:`k_f`", "kf", "None", ":math:`\text{m/s}`"
   "Surface diffusion coefficient", ":math:`D_s`", "ds", "None", ":math:`\text{m}^2\text{/s}`"
   "Solute distribution parameter", ":math:`D_g`", "dg", "None", ":math:`\text{dimensionless}`"
   "Biot number", ":math:`Bi`", "N_Bi", "None", ":math:`\text{dimensionless}`"
   "Minimum Stanton number for CPS", ":math:`St_{min}`", "min_N_St", "None", ":math:`\text{dimensionless}`"
   "Minimum empty bed contact time for CPS", ":math:`EBCT_{min}`", "min_ebct", "None", ":math:`\text{s}`"
   "Minimum residence time for CPS", ":math:`\tau_{min}`", "min_res_time", "None", ":math:`\text{s}`"
   "Minimum elapsed operation time for CPS", ":math:`t_{min}`", "min_elap_time", "None", ":math:`\text{s}`"
   "Mass throughput ratio for the upstream edge of the MTZ", ":math:`T|\left( \frac{C}{C_{0}} = 0.95 \right)`", "mass_throughput_mtz_upstream", "None", ":math:`\text{dimensionless}`"
   "Empty bed contact time of the partial MTZ at the time of replacement", ":math:`EBCT_{MTZ}`", "ebct_mtz_replace", "None", ":math:`\text{s}`"
   "Length of the partial MTZ at the time of replacement", ":math:`L_{MTZ}`", "length_mtz_replace", "None", ":math:`\text{m}`"
   "Stanton equation parameter 0", ":math:`a_0`", "a0", "None", ":math:`\text{dimensionless}`"
   "Stanton equation parameter 1", ":math:`a_1`", "a1", "None", ":math:`\text{dimensionless}`"
   "Throughput equation parameter 0", ":math:`b_0`", "b0", "None", ":math:`\text{dimensionless}`"
   "Throughput equation parameter 1", ":math:`b_1`", "b1", "None", ":math:`\text{dimensionless}`"
   "Throughput equation parameter 2", ":math:`b_2`", "b2", "None", ":math:`\text{dimensionless}`"
   "Throughput equation parameter 3", ":math:`b_3`", "b3", "None", ":math:`\text{dimensionless}`"
   "Throughput equation parameter 4", ":math:`b_4`", "b4", "None", ":math:`\text{dimensionless}`"
   "Molecular diffusion coefficient", ":math:`D_l`", "diffus_liq", "None", ":math:`\text{m}^2\text{\s}`"
   "Molal Volume", ":math:`V_b`", "molal_volume", "None", ":math:`\text{m}^3\text{/mol}`"
   "Reynolds number", ":math:`Re`", "N_Re", "None", ":math:`\text{dimensionless}`"
   "Schmidt number", ":math:`Sc`", "N_Sc", "None", ":math:`\text{dimensionless}`"
   "GAC particle sphericity", ":math:`\phi`", "sphericity", "None", ":math:`\text{dimensionless}`"

.. _GAC_equations:

Equations
-----------


Class Documentation
-------------------

* :mod:`watertap.unit_models.gac`
