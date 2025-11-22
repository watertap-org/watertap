.. _GAC:

Granular Activated Carbon (GAC)
===============================

.. code-block:: python

   from watertap.unit_models.gac import GAC

This is an empirical, performance-based granular activated carbon (GAC) model that works under the following criteria and assumptions:
   * simulation of this unit model is only supported with the Multi-Component Aqueous Solution (MCAS) property package
   * supports a single liquid phase only
   * supports adsorption of a single solute species only while other species are considered inert
   * supports steady-state only
   * assumes isothermal conditions
   * model performance is independent of a gravity-fed or pressurized GAC unit, therefore assumes isobaric conditions

.. index::
   pair: watertap.unit_models.gac;gac

.. currentmodule:: watertap.unit_models.gac

Introduction
------------
The implemented model for estimating GAC performance is adapted from a simplified model originally presented in Hand, 1984 and
further elaborated in Crittenden, 2012. This formulation is denoted as the constant-pattern homogeneous surface diffusion model (CPHSDM).
As a GAC system is operated as a batch process, a mass transfer zone (MTZ) is formed in the bed where a concentration profile, or
breakthrough curve, develops as a function of the adsorption properties. This MTZ is bounded by saturated GAC upstream and
fresh GAC downstream. The CPHSDM is valid under the assumption that the shape of the MTZ is constant as it travels through the
bed and a constant pattern solution (CPS) may be determined. The CPS is calculated through a multistep procedure utilizing common
dimensionless groups applied in polynomial fits to determine performance. Therefore, coefficients used in the polynomial must be
derived from experimental data of the intended system to produce valid results. Coefficients for common compounds treated by GAC may be
found in both Hand, 1984 and Crittenden, 2012, replicated in :ref:`Table 1 <min_St_params>` and :ref:`Table 2 <throughout_params>` below.
The model is estimated to have within 10% error and therefore may be applied to bed lengths shorter than the minimum
length determined by the CPHSDM within the error threshold (in addition to being applicable to bed lengths greater than the minimum
length determined by the CPHSDM). As an alternative to the polynomial fits for these dimensionless groups, the WaterTAP GAC model 
includes surrogate models that can be used to predict performance without the need for the polynomial coefficients.

Model Structure
------------------
The GAC model consists of one ControlVolume0DBlock (``process_flow``) for the process flow of the water treatment train.
The process flow includes two StateBlocks (``inlet`` and ``outlet``) which are used for mass and momentum balances.
It also includes one StateBlock (``adsorbed``) for the solute that is adsorbed into the GAC particles.
The material removed in the ``adsorbed`` state block is simulated as liquid phase solute but should be interpreted as solute that has adsorbed
into the solid phase GAC particles. The steady state mass removal and replacement rate of the GAC itself is provided as a unit model
variable and excluded from flowsheet material balances. Therefore, GAC is never included as a component in property package specifications.

Effluent Concentration
----------------------
The model includes two approaches for determining the effluent concentration of the model as passed via the
``add_trapezoidal_effluent_approximation`` configuration argument. If ``add_trapezoidal_effluent_approximation`` is set to ``True`` (the default setting),
the batch operation results of the CPS are converted to approximate steady-state results for intuitive use of the model
for flowsheet purposes. A visualization of the transformation is provided in Figure 1. For a traditional breakthrough
curve, the CPHSDM method calculates a single value for ``conc_ratio_replace`` and ``operational_time``, highlighted on the
breakthrough curve. This operational time is the amount of elapsed time after startup that the bed is refreshed with new
GAC adsorbent. Steady state concentration can be analogous to all of the effluent in this operational time being stored
as holdup, therefore the average concentration ratio is significantly less than concentration ratio at the time of
bed replacement, as many days pass before the start of the breakthrough. To approximate the average effluent
concentration in this time frame, the breakthrough curve is numerically integrated with the trapezoid rule. The curve is
discretized with respect to the concentration ratio instead of the (traditionally done 'x' variable) operational time
due to simplicity of solving the model equations.

.. figure:: ../../_static/unit_models/gac.png
    :width: 800
    :align: center

    Figure 1. Discretization of the breakthrough curve for the steady state approximation. The ratio of the shaded area
    to the area left of the vertical dotted line corresponding to ``operational_time`` is the calculated
    ``conc_ratio_avg``. Expected values of ``conc_ratio_avg`` are often less than 0.25 depending on the
    ``conc_ratio_replace`` setpoint.

If ``add_trapezoidal_effluent_approximation`` is set to ``False``, the model sets the effluent concentration to that calculated via ``conc_ratio_replace``.
This configuration may not be suitable for a flowsheet applications, but the lower model complexity can be beneficial for certain use cases,
such as preliminary design or parameter estimation.

CPHSDM Calculations
-------------------
The GAC model relies on calculation of a minimum Stanton number and throughput parameter to determine the CPS.
The original presentation of the CPHSDM in Hand, 1984 includes a polynomial fit of these two parameters, with
regressed coefficients presented in :ref:`Table 1 <min_St_params>` and :ref:`Table 2 <throughout_params>`, adopted from Hand, 1984.
The WaterTAP GAC model includes two options for calculating the minimum Stanton number and throughput parameter, passed via the 
``cphsdm_calculation_method`` argument in the model configuration:

    * ``input``: The user fixes the necessary parameters to the proper values.
    * ``surrogate``: The model uses a surrogate model to calculate the necessary parameters.

Because the parameters are only provided for discrete values of the input variables, the surrogate models can be used for values between these discrete points.
The surrogate models were generated by first calculating :math:`St_{min}` and :math:`T` for all the discrete points, then using that 
data to create surrogate models using `PySMO <https://idaesplus.readthedocs.io/projects/idaes/en/latest/explanations/modeling_extensions/surrogate/api/pysmo/pysmo_radialbasisfunctions.html>`_. 
These two surrogate models take the general form:

.. math::

   St_{min} = f\left(\frac{1}{n}, Bi\right)

   T = g\left(\frac{1}{n}, Bi, \frac{C}{C_0}\right)

The throughput surrogate is used for both the throughput at the relative effluent concentration of interest (``conc_ratio_replace``) and for the
discretized relative effluent concentration points used to estimate steady-state performance (``ele_conc_ratio_replace[i]``).
Surrogate models are loaded from files in the `watertap/data/surrogate_defaults/gac` directory. They are accessible from the
unit model as e.g., ``m.fs.unit.min_N_St_surrogate`` and ``m.fs.unit.throughput_surrogate``. The surrogate constraints themselves are on
``SurrogateBlock`` objects accessible from the unit model as e.g., ``m.fs.unit.min_N_St_surrogate_blk`` and ``m.fs.unit.throughput_surrogate_blk``.


Degrees of Freedom
------------------
In the default configuration of the GAC unit model there are 17 degrees of freedom in addition to the inlet state variables
(i.e., temperature, pressure, component flowrates) that should be fixed for the model to be fully specified.
In association with using the Freundlich adsorption isotherm and empirical model, the following 2 equilibrium variables are almost always fixed
and may be derived from experimental data:

   * Freundlich isotherm :math:`k` parameter
   * Freundlich isotherm :math:`\frac{1}{n}` parameter

If ``cphsdm_calculation_method`` is set to ``input``, the following variables should also be fixed:

   * Stanton number equation parameters :math:`a_0` and :math:`a_1` (Hand, 1984), (Crittenden, 1987)
   * throughput ratio equation parameters :math:`b_0`, :math:`b_1`, :math:`b_2`, :math:`b_3` and :math:`b_4` (Hand, 1984)

Additionally, the following 9 variables are traditionally fixed:

   * particle apparent density
   * particle diameter
   * empty bed contact time
   * bed voidage *or* particle bulk density
   * superficial velocity *or* bed length
   * effluent to inlet concentration ratio at operational time *or* steady state approximation of average effluent to inlet concentration ratio in operational time by trapezoid rule *or* bed volumes treated
   * surface diffusion coefficient
   * liquid phase film transfer coefficient

When setting the configuration options to calculate the surface diffusion coefficient and liquid phase film transfer coefficient,
these respective variables are no longer specified and 4 newly introduced variables must be fixed. This excludes new variables
or parameters that may be required to be specified within the property package when called by the GAC model. This is a net result of 19 degrees of freedom.
Newly utilized variables that must be fixed include:

   * shape correction factor
   * particle porosity
   * tortuosity of the path that the adsorbate must take as compared to the radius
   * surface-to-pore diffusion flux ratio

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Time", ":math:`t`", "[0]"
   "Phases", ":math:`p`", "['Liq']"
   "Components", ":math:`j`", "['H2O', target_species, background solutes]*"
   "Species adsorbed", ":math:`\text{target_species}`", "[target_species]"
   "Inert species", ":math:`\text{inert_species}`", "['H2O', target_species, background solutes] - [target_species]"


If ``add_trapezoidal_effluent_approximation`` is set to ``True``:

.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "Number of discretized operational time elements used for steady state approximation", ":math:`\text{ele_disc}`", "[0:elements_ss_approx]"
   "Number of discretized trapezoidal area terms for steady state approximation", ":math:`\text{ele_index}`", "[1:elements_ss_approx]"

| \* target_species is provided in the ``target_species`` argument of the unit model and corresponds to the single solute which is adsorbed.
| \* ``inert_species`` are the difference in ``component_list - target_species``.

.. _GAC_variables:

Variables
----------
Supporting only a single solute, variables concerning the adsorption of a species are respective to the target species
although not explicitly indexed. All other species are inert from a mass balance perspective, but effects of background
species to the adsorption of the target species may be modeled by adjusting the Freundlich isotherm parameters and other
variables in the model.

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Freundlich isotherm k parameter", ":math:`k`", "freund_k", "None", ":math:`\left(\text{m}^3\text{/kg}\right)^\left( \frac{1}{n} \right)`"
   "Freundlich isotherm 1/n parameter", ":math:`\frac{1}{n}`", "freund_ninv", "None", ":math:`\text{dimensionless}`"
   "Surface diffusion coefficient", ":math:`D_s`", "ds", "None", ":math:`\text{m}^2\text{/s}`"
   "Liquid phase film transfer coefficient", ":math:`k_f`", "kf", "None", ":math:`\text{m/s}`"
   "Equilibrium concentration of adsorbed phase with liquid phase", ":math:`q_e`", "equil_conc", "None", ":math:`\left( \text{kg}_\text{adsorbate}\text{/kg}_\text{adsorbent} \right)`"
   "Solute distribution parameter", ":math:`D_g`", "dg", "None", ":math:`\text{dimensionless}`"
   "Biot number", ":math:`Bi`", "N_Bi", "None", ":math:`\text{dimensionless}`"
   "Superficial velocity", ":math:`u_s`", "velocity_sup", "None", ":math:`\text{m/s}`"
   "Interstitial velocity", ":math:`u_i`", "velocity_int", "None", ":math:`\text{m/s}`"
   "Bed void fraction", ":math:`\epsilon`", "bed_voidage", "None", ":math:`\text{dimensionless}`"
   "Bed length", ":math:`L`", "bed_length", "None", ":math:`\text{m}`"
   "Bed diameter", ":math:`D`", "bed_diameter", "None", ":math:`\text{m}`"
   "Bed area", ":math:`A`", "bed_area", "None", ":math:`\text{m}^2`"
   "Bed volume", ":math:`V`", "bed_volume", "None", ":math:`\text{m}^3`"
   "Empty bed contact time", ":math:`EBCT`", "ebct", "None", ":math:`\text{s}`"
   "Fluid residence time in the bed", ":math:`\tau`", "residence_time", "None", ":math:`\text{s}`"
   "Mass of fresh GAC in the bed", ":math:`M_{GAC}`", "bed_mass_gac", "None", ":math:`\text{kg}`"
   "GAC apparent density", ":math:`\rho_a`", "particle_dens_app", "None", ":math:`\text{kg/}\text{m}^3`"
   "GAC bulk density", ":math:`\rho_b`", "particle_dens_bulk", "None", ":math:`\text{kg/}\text{m}^3`"
   "GAC particle diameter", ":math:`d_p`", "particle_dia", "None", ":math:`\text{m}`"
   "Minimum Stanton number to achieve a constant pattern solution", ":math:`St_{min}`", "min_N_St", "None", ":math:`\text{dimensionless}`"
   "Minimum empty bed contact time to achieve a constant pattern solution", ":math:`EBCT_{min}`", "min_ebct", "None", ":math:`\text{s}`"
   "Specific throughput", ":math:`T`", "throughput", "None", ":math:`\text{dimensionless}`"
   "Minimum fluid residence time in the bed to achieve a constant pattern solution", ":math:`\tau_{min}`", "min_residence_time", "None", ":math:`\text{s}`"
   "Minimum operational time of the bed from fresh to achieve a constant pattern solution", ":math:`t_{min}`", "min_operational_time", "None", ":math:`\text{s}`"
   "Effluent to inlet concentration ratio at operational time", ":math:`\frac{C}{C_{0}}\bigg{|}_{z=L,/,t=t_{op}}`", "conc_ratio_replace", "None", ":math:`\text{dimensionless}`"
   "Operational time of the bed from fresh", ":math:`t_{op}`", "operational_time", "None", ":math:`\text{s}`"
   "Bed volumes treated at operational time", ":math:`BVT`", "bed_volumes_treated", "None", ":math:`\text{dimensionless}`"

The following variables are only built when specific configuration options are selected.

If ``add_trapezoidal_effluent_approximation`` is set to ``True``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Specific throughput from empirical equation by discrete element", ":math:`T_{ele}`", "ele_throughput", "None", ":math:`x`"
   "Minimum operational time of the bed from fresh to achieve a constant pattern solution by discrete element", ":math:`t_{min, ele}`", "ele_min_operational_time", "ele_index", ":math:`\text{s}`"
   "Effluent to inlet concentration ratio at operational time by discrete element", ":math:`\left(\frac{C}{C_{0}}\right)_{ele}\bigg{|}_{z=L,/,t=t_{op_ e}}`", "ele_conc_ratio_replace", "ele_index", ":math:`\text{dimensionless}`"
   "Operational time of the bed from fresh by discrete element", ":math:`t_{op, ele}`", "ele_operational_time", "ele_disc", ":math:`\text{s}`"
   "Trapezoid rule of elements for numerical integration of average concentration ratio", ":math:`term_{ele}`", "ele_conc_ratio_avg", "ele_disc", ":math:`\text{dimensionless}`"
   "Steady state approximation of average effluent to inlet concentration ratio in operational time by trapezoid rule", ":math:`\left(\frac{C}{C_{0}}\right)_{avg}`", "conc_ratio_avg", "None", ":math:`\text{dimensionless}`"
   "Total mass of adsorbed species at operational time", ":math:`M`", "mass_adsorbed", "None", ":math:`\text{kg}`"
   "GAC usage/replacement/regeneration rate", ":math:`\dot{m}_{GAC}`", "gac_usage_rate", "None", ":math:`\text{m/s}`"

If ``cphsdm_calculation_method`` is set to ``input``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Stanton equation parameter 0", ":math:`a_0`", "a0", "None", ":math:`\text{dimensionless}`"
   "Stanton equation parameter 1", ":math:`a_1`", "a1", "None", ":math:`\text{dimensionless}`"
   "Throughput equation parameter 0", ":math:`b_0`", "b0", "None", ":math:`\text{dimensionless}`"
   "Throughput equation parameter 1", ":math:`b_1`", "b1", "None", ":math:`\text{dimensionless}`"
   "Throughput equation parameter 2", ":math:`b_2`", "b2", "None", ":math:`\text{dimensionless}`"
   "Throughput equation parameter 3", ":math:`b_3`", "b3", "None", ":math:`\text{dimensionless}`"
   "Throughput equation parameter 4", ":math:`b_4`", "b4", "None", ":math:`\text{dimensionless}`"

If ``film_transfer_coefficient_type`` is set to ``calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Reynolds number", ":math:`Re`", "N_Re", "None", ":math:`\text{dimensionless}`"
   "Schmidt number", ":math:`Sc`", "N_Sc", "None", ":math:`\text{dimensionless}`"
   "Shape correction factor", ":math:`SCF`", "shape_correction_factor", "None", ":math:`\text{dimensionless}`"

If ``surface_diffusion_coefficient_type`` is set to ``calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "GAC particle porosity", ":math:`\epsilon_p`", "particle_porosity", "None", ":math:`\text{dimensionless}`"
   "Tortuosity of the path that the adsorbate must take as compared to the radius", ":math:`\tau_p`", "tort", "None", ":math:`\text{dimensionless}`"
   "Surface-to-pore diffusion flux ratio", ":math:`S\!P\!D\!F\!R`", "spdfr", "None", ":math:`\text{dimensionless}`"

.. _GAC_equations:

Equations
-----------

.. csv-table::
   :header: "Description", "Equation"

   "Equilibrium concentration", ":math:`q_e = kC_0^{1/n}`"
   "Solute distribution parameter", ":math:`D_g=\frac{\rho_aq_e\left( 1-\epsilon \right)}{\epsilon C_0}`"
   "Biot number", ":math:`Bi=\frac{k_fd_p\left( 1-\epsilon \right)}{2D_sD_g\epsilon}`"
   "Bed void fraction based on GAC particle densities", ":math:`\epsilon=1-\frac{\rho_b}{\rho_a}`"
   "Relating velocities based on bed voidage", ":math:`u_i=\frac{u_s}{\epsilon}`"
   "Bed length based on velocity and EBCT", ":math:`L=(EBCT)u_s`"
   "Bed diameter and area relation", ":math:`A=\pi\left(\frac{D}{2}\right)^2`"
   "Bed area based on velocity and volumetric flow", ":math:`A=\frac{Q}{u_s}`"
   "Bed volume based on cylindrical dimensions", ":math:`V=AL`"
   "Fluid residence time in the bed", ":math:`\tau=(EBCT)\epsilon`"
   "Total mass of GAC in the bed", ":math:`M_{GAC}=V\rho_b`"
   "Minimum empty bed contact time to achieve constant pattern solution", ":math:`EBCT_{min}=\frac{St_{min}d_p}{2k_f\left( 1-\epsilon \right)}`"
   "Minimum fluid residence time in the bed to achieve a constant pattern solution", ":math:`\tau_{min}=EBCT_{min}\epsilon`"
   "Minimum operational time of the bed from fresh to achieve a constant pattern solution", ":math:`t_{min}=\tau_{min}\left( D_g+1 \right)T`"
   "Elapsed operational time between a fresh bed and the theoretical bed replacement", ":math:`t_{op}=t_{min}+\left( \tau-\tau_{min} \right)\left( D_g+1 \right)`"
   "Bed volumes treated", ":math:`BVT=\frac{t_{op}\epsilon}{\tau}`"

If ``add_trapezoidal_effluent_approximation`` is set to ``True``:

.. csv-table::
   :header: "Description", "Equation"

   "Minimum operational time of the bed from fresh to achieve a constant pattern solution by discretized element", ":math:`t_{min, ele}=\tau_{min}\left( D_g+1 \right)T`"
   "Creating evenly spaced discretized elements", ":math:`\frac{C}{C_{0}}\bigg{|}_{t=t_{op, ele}}=0.01+(ele-1)*\frac{\left(\frac{C}{C_{0}}\bigg{|}_{t=t_{op}}-0.01\right)}{num\text{_}ele}`"
   "Finite element discretization of concentration ratios over time", ":math:`term_{ele}=\left(\frac{t_{op, ele}-t_{op, (ele-1)}}{t_{op}}\right)\frac{\left(\frac{C}{C_{0}}\bigg{|}_{t=t_{op, ele}}+{C_{0}}\bigg{|}_{t=t_{op, (ele-1)}}\right)}{2}`"
   "Summation of finite elements for average concentration during operating time", ":math:`\left(\frac{C}{C_{0}}\right)_{avg}=\sum_{ele\text{_}index}term_{ele}`"
   "Mass adsorbed in the operational time", ":math:`M=\frac{\dot{m}_{j}}{t_{op}}`"
   "Steady state rate of new GAC mass required", ":math:`\dot{m}_{GAC}=\frac{M_{GAC}}{t_{op}}`"


If ``cphsdm_calculation_method`` is set to ``input``:

.. csv-table::
   :header: "Description", "Equation"

   "Minimum Stanton number to achieve constant pattern solution", ":math:`St_{min}=a_0Bi+a_1`"
   "Throughput based on empirical 5-parameter regression", ":math:`T=b_0+b_1\left( \frac{C}{C_0} \right)^{b_2}+\frac{b_3}{1.01-\left( \frac{C}{C_0} \right)^{b_4}}`"
   "Throughput based on empirical 5-parameter regression by discretized element", ":math:`T_{ele}=b_0+b_1\left(\frac{C}{C_{0}}\right)_{ele}^{b_2}+\frac{b_3}{1.01-\left(\frac{C}{C_{0}}\right)_{ele}^{b_4}}`"


If ``film_transfer_coefficient_type`` is set to ``calculated``:

.. csv-table::
   :header: "Description", "Equation"

   "Reynolds number calculation*", ":math:`Re=\frac{\rho_ld_pu_i}{\mu_l}`"
   "Schmidt number calculation*", ":math:`Sc=\frac{\mu_l}{\rho_sD_l}`"
   "Liquid phase film transfer rate from the Gnielinski correlation*", ":math:`k_f=(SCF)\frac{\left[ 1+1.5\left( 1-\epsilon \right) \right]D_l}{d_p}\left( 2+0.644Re^{\frac{1}{2}}Sc^{\frac{1}{3}} \right)`"

\*Subscript :math:`l` denotes bulk liquid phase properties, here those are supplied by the property package.


If ``surface_diffusion_coefficient_type`` is set to ``calculated``:

.. csv-table::
   :header: "Description", "Equation"

   "surface diffusion parameter (Crittenden, 1987)", ":math:`D_s=\left( S\!P\!D\!F\!R \right)\left( \frac{\epsilon_pC_0D_l}{\rho_aq_e\tau_p} \right)`"


.. _GAC_CPHSDM_params:

CPHSDM Empirical Parameters
--------------------------- 

If ``cphsdm_calculation_method`` is set to ``input``, the following parameters can be used for ``a0`` and ``a1`` to calculate the minimum Stanton number 
and ``b0``, ``b1``, ``b2``, ``b3``, and ``b4`` for the throughput equation. These tables were adapted directly from Hand, 1984.


.. _min_St_params:

.. table:: Table 1: Empirical Parameters for Minimum Stanton Number Equation
   :align: left

   +------------------+---------------+---------+-------------+--------+
   |  ``freund_ninv`` | 0.5 ≤ Bi ≤ 10           | 10 ≤ Bi ≤ ∞          |
   +                  +---------------+---------+-------------+--------+
   |                  | ``a0``        | ``a1``  | ``a0``      | ``a1`` |
   +==================+===============+=========+=============+========+
   | 0.05             | 0.02105       | 1.98947 | 0.22        | 0      |
   +------------------+---------------+---------+-------------+--------+
   | 0.1              | 0.02105       | 2.18947 | 0.24        | 0      |
   +------------------+---------------+---------+-------------+--------+
   | 0.2              | 0.04211       | 2.37895 | 0.28        | 0      |
   +------------------+---------------+---------+-------------+--------+
   | 0.3              | 0.10526       | 2.54737 | 0.36        | 0      |
   +------------------+---------------+---------+-------------+--------+
   | 0.4              | 0.23158       | 2.68421 | 0.50        | 0      |
   +------------------+---------------+---------+-------------+--------+
   | 0.5              | 0.52632       | 2.73684 | 0.80        | 0      |
   +------------------+---------------+---------+-------------+--------+
   | 0.6              | 1.15789       | 3.42105 | 1.50        | 0      |
   +------------------+---------------+---------+-------------+--------+
   | 0.7              | 1.78947       | 7.10526 | 2.50        | 0      |
   +------------------+---------------+---------+-------------+--------+
   | 0.8              | 3.68421       | 13.1579 | 5.00        | 0      |
   +------------------+---------------+---------+-------------+--------+
   | 0.9              | 6.31579       | 56.8421 | 12.00       | 0      |
   +------------------+---------------+---------+-------------+--------+



.. _throughout_params:

.. table:: Table 2: Empirical Parameters for Throughput Equation
   :align: left

   +------------------+------+-----------+----------+----------+----------+-----------+
   |  ``freund_ninv`` | Bi   | ``b0``    | ``b2``   | ``b2``   | ``b3``   | ``b4``    |
   +==================+======+===========+==========+==========+==========+===========+
   | 0.05             | 0.5  | -5.447214 | 6.598598 | 0.026569 | 0.019384 | 20.450470 |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 2    | -5.465811 | 6.592484 | 0.004989 | 0.004988 | 0.503250  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 4    | -5.531155 | 6.584935 | 0.023580 | 0.009019 | 0.273076  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 6    | -5.606508 | 6.582188 | 0.022088 | 0.013126 | 0.214246  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 8    | -5.606500 | 6.504701 | 0.020872 | 0.017083 | 0.189537  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 10   | -5.664173 | 6.456597 | 0.018157 | 0.019935 | 0.149314  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 14   | -0.662780 | 1.411252 | 0.060709 | 0.020229 | 0.143293  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 25   | -0.662783 | 1.350940 | 0.031070 | 0.020350 | 0.129998  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | ≥100 | 0.665879  | 0.711310 | 2.987309 | 0.016783 | 0.361023  |
   +------------------+------+-----------+----------+----------+----------+-----------+
   | 0.1              | 0.5  | -1.919873 | 3.055368 | 0.055488 | 0.024284 | 15.311766 |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 2    | -2.278950 | 3.393925 | 0.046838 | 0.004751 | 0.384675  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 4    | -2.337178 | 3.379926 | 0.043994 | 0.008650 | 0.243412  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 6    | -2.407407 | 3.374131 | 0.041322 | 0.012552 | 0.196565  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 8    | -2.477819 | 3.370954 | 0.038993 | 0.016275 | 0.176437  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 10   | -2.566414 | 3.370950 | 0.035003 | 0.019386 | 0.150788  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 16   | -2.567201 | 3.306341 | 0.020940 | 0.019483 | 0.136813  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 30   | -2.568618 | 3.241783 | 0.009595 | 0.019610 | 0.121746  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | ≥100 | -2.568360 | 3.191482 | 0.001555 | 0.019682 | 0.110113  |
   +------------------+------+-----------+----------+----------+----------+-----------+
   | 0.2              | 0.5  | -1.441000 | 2.569000 | 0.060920 | 0.002333 | 0.371100  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 2    | -1.474313 | 2.558300 | 0.058480 | 0.005026 | 0.241265  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 4    | -1.506696 | 2.519259 | 0.055525 | 0.008797 | 0.187510  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 6    | -1.035395 | 1.983018 | 0.069283 | 0.012302 | 0.167924  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 8    | -0.169192 | 1.077521 | 0.144879 | 0.015500 | 0.168083  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 10   | -1.402932 | 2.188339 | 0.052191 | 0.018422 | 0.133574  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 13   | -1.369220 | 2.118545 | 0.039492 | 0.018453 | 0.127565  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 25   | -1.514159 | 2.209450 | 0.017937 | 0.018510 | 0.118517  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | ≥100 | 0.680346  | 0.649006 | 2.570086 | 0.014947 | 0.369818  |
   +------------------+------+-----------+----------+----------+----------+-----------+
   | 0.3              | 0.5  | -1.758696 | 2.846576 | 0.049530 | 0.003022 | 0.156816  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 2    | -1.657862 | 2.688895 | 0.048409 | 0.005612 | 0.140937  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 4    | -0.565664 | 1.537833 | 0.084451 | 0.008808 | 0.199086  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 6    | -0.197077 | 1.118564 | 0.117894 | 0.011527 | 0.135874  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 8    | -0.197070 | 1.069216 | 0.119760 | 0.013925 | 0.132691  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 10   | -0.173358 | 1.000000 | 0.120311 | 0.015940 | 0.133973  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 15   | -0.173350 | 0.919411 | 0.071768 | 0.014156 | 0.086270  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 35   | 0.666471  | 0.484570 | 7.719440 | 0.013444 | 0.259545  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | ≥100 | 0.696161  | 0.516951 | 2.054587 | 0.012961 | 0.303218  |
   +------------------+------+-----------+----------+----------+----------+-----------+
   | 0.4              | 0.5  | -0.534251 | 1.603834 | 0.094055 | 0.004141 | 0.137797  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 2    | -0.166270 | 1.190897 | 0.122280 | 0.006261 | 0.134278  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 4    | -0.166270 | 1.131946 | 0.115513 | 0.008634 | 0.126813  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 6    | -0.166270 | 1.089789 | 0.112284 | 0.010463 | 0.124307  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 9    | 0.491912  | 0.491833 | 0.487414 | 0.011371 | 0.147747  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 12   | 0.564119  | 0.419196 | 0.639819 | 0.011543 | 0.149005  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 15   | 0.640669  | 0.432466 | 1.048056 | 0.011616 | 0.212726  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 25   | 0.672353  | 0.397007 | 1.153169 | 0.011280 | 0.216883  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | ≥100 | 0.741435  | 0.448054 | 1.929879 | 0.010152 | 0.306448  |
   +------------------+------+-----------+----------+----------+----------+-----------+
   | 0.5              | 0.5  | -0.040800 | 1.099652 | 0.158995 | 0.005467 | 0.139116  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 4    | -0.040800 | 0.982757 | 0.111618 | 0.008072 | 0.111404  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 10   | 0.094602  | 0.754878 | 0.092069 | 0.009877 | 0.090763  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 14   | 0.023000  | 0.802068 | 0.057545 | 0.009662 | 0.084532  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 25   | 0.023000  | 0.793673 | 0.039324 | 0.009326 | 0.082751  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | ≥100 | 0.529213  | 0.291801 | 0.082428 | 0.008317 | 0.075461  |
   +------------------+------+-----------+----------+----------+----------+-----------+
   | 0.6              | 0.5  | 0.352536  | 0.692114 | 0.263134 | 0.005482 | 0.121775  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 2    | 0.521979  | 0.504220 | 0.327290 | 0.005612 | 0.128679  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 6    | 0.676253  | 0.334583 | 0.482297 | 0.005898 | 0.138946  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 14   | 0.769531  | 0.259497 | 0.774068 | 0.005600 | 0.165513  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 50   | 0.849057  | 0.215799 | 1.343183 | 0.004725 | 0.223759  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | ≥100 | 0.831231  | 0.227304 | 1.174756 | 0.004961 | 0.212109  |
   +------------------+------+-----------+----------+----------+----------+-----------+
   | 0.7              | 0.5  | 0.575024  | 0.449062 | 0.278452 | 0.004122 | 0.121682  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 4    | 0.715269  | 0.307172 | 0.442104 | 0.004371 | 0.138351  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 12   | 0.787940  | 0.243548 | 0.661599 | 0.004403 | 0.162595  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 25   | 0.829492  | 0.204078 | 0.784529 | 0.004050 | 0.179003  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | ≥100 | 0.847012  | 0.190678 | 0.931686 | 0.003849 | 0.183239  |
   +------------------+------+-----------+----------+----------+----------+-----------+
   | 0.8              | 0.5  | 0.708905  | 0.314101 | 0.357499 | 0.003276 | 0.119300  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 4    | 0.784576  | 0.239663 | 0.484422 | 0.003206 | 0.134987  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 14   | 0.839439  | 0.188966 | 0.648124 | 0.003006 | 0.157697  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | ≥100 | 0.882747  | 0.146229 | 0.807987 | 0.002537 | 0.174543  |
   +------------------+------+-----------+----------+----------+----------+-----------+
   | 0.9              | 0.5  | 0.865453  | 0.157618 | 0.444973 | 0.001650 | 0.148084  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 4    | 0.854768  | 0.171434 | 0.495042 | 0.001910 | 0.142251  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | 16   | 0.866180  | 0.163992 | 0.573946 | 0.001987 | 0.157594  |
   +                  +------+-----------+----------+----------+----------+-----------+
   |                  | ≥100 | 0.893192  | 0.133039 | 0.624100 | 0.001740 | 0.164248  |
   +------------------+------+-----------+----------+----------+----------+-----------+

Code Documentation
-------------------

* :mod:`watertap.unit_models.gac`
* :mod:`watertap.costing.unit_models.gac`

References
-----------
| Hand, D. W., Crittenden, J. C., & Thacker, W. E. (1984). 
| Simplified models for design of fixed-bed adsorption systems.
| *Journal of Environmental Engineering*, 110(2), 440-456.
|
| Crittenden, J., Rhodes, R., Hand, D., Howe, K., & Tchobanoglous, G. (2012). 
| MWHs Water Treatment. Principles and Design.
| John Wiley & Sons.

| Crittenden, J. C., Berrigan, J. K., Hand, D. W., & Lykins, B. (1987). 
| Design of Rapid Fixed‐Bed Adsorption Tests for Nonconstant Diffusivities. 
| *Journal of Environmental Engineering*, 113(2), 243–259.
