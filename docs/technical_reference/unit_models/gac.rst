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
found in both Hand, 1984 and Crittenden, 2012. The model is estimated to have within 10% error and therefore may be applied to bed lengths
shorter than the minimum length determined by the CPHSDM within the error threshold (in addition to being applicable to bed lengths greater than the minimum length determined by the CPHSDM).

The batch operation results of the CPS are converted to approximate steady-state results for intuitive use of the model
for flowsheet purposes. A visualization of the transformation is provided in Figure 1. For a traditional breakthrough
curve the CPHSDM method calculates the single point, single ``conc_ratio_replace`` and ``operational_time``, highlighted on the
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

Degrees of Freedom
------------------
In the default configuration of the GAC unit model there are 17 degrees of freedom in addition to the inlet state variables
(i.e., temperature, pressure, component flowrates) that should be fixed for the model to be fully specified.
In association with using the Freundlich adsorption isotherm and empirical model, the following 9 variables are almost always fixed
and may be derived from experimental data:

   * Freundlich isotherm :math:`k` parameter
   * Freundlich isotherm :math:`\frac{1}{n}` parameter
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

Model Structure
------------------
The GAC model consists of 1 ControlVolume0DBlock (``process_flow``) for the process flow of the water treatment train.
The process flow includes 2 StateBlocks (``inlet`` and ``outlet``) which are used for mass and momentum balances.
It also includes 1 StateBlock (``adsorbed``) for the solute that is adsorbed into the GAC particles.
The material removed in the ``adsorbed`` state block is simulated as liquid phase solute but should be interpreted as solute that has adsorbed
into the solid phase GAC particles. The steady state mass removal and replacement rate of the GAC itself is provided as a unit model
variable and excluded from flowsheet material balances. Therefore, GAC is never included as a component in property package specifications.


CPHSDM Calculations
-------------------
The GAC model relies on calculation of a minimum Stanton number and throughput to determine the CPS.
The original presentation of the CPHSDM in Hand, 1984 and Crittenden, 2012 is based on a polynomial fit of these two parameters, with 
regressed parameters presented in :ref:`Table 1 <min_St_params>` and :ref:`Table 2 <throughout_params>`, adopted from Hand, 1984.
The WaterTAP GAC model includes two options for calculating the minimum Stanton number and throughput parameter, passed via the 
``cphsdm_calculation_method`` argument in the model configuration:

    * ``input``: The user fixes the necessary parameters to the proper values.
    * ``surrogate``: The model uses a surrogate model to calculate the necessary parameters.

Because the parameters are only provided for discrete values of the input variables, the surrogate models can be used for values between these discrete points.
The surrogate models were generated by first calculating :math:`St_{min}` and :math:`T` for all the discrete points, then using that 
data to create surrogate models using a radial basis function `PySMO <https://idaesplus.readthedocs.io/projects/idaes/en/latest/explanations/modeling_extensions/surrogate/api/pysmo/pysmo_radialbasisfunctions.html>`_ trainer. 
These two surrogate models take the general form:

.. math::

   St_{min} = f\left(\frac{1}{n}, Bi\right)

   T = g\left(\frac{1}{n}, Bi, \frac{C}{C_0}\right)

The throughput surrogate is used for both the throughput at the effluent concentration of interest (``conc_ratio_replace``) and for all the
discretized points used to estimate steady-state performance (``ele_conc_ratio_replace[i]``)
Surrogate models are loaded from files in the `watertap/data/surrogate_defaults/gac` directory. They are accessible from the
unit model as e.g., ``m.fs.unit.min_N_St_surrogate`` and ``m.fs.unit.throughput_surrogate``. The surrogate constraints themselves are on
``SurrogateBlock`` objects accessible from the unit model as e.g., ``m.fs.unit.min_N_St_surrogate_blk`` and ``m.fs.unit.throughput_surrogate_blk``.

Sets
----
.. csv-table::
   :header: "Description", "Symbol", "Indices"

   "time", ":math:`t`", "[0]"
   "phases", ":math:`p`", "['Liq']"
   "components", ":math:`j`", "['H2O', target_species, background solutes]*"
   "species adsorbed", ":math:`\text{target_species}`", "[target_species]"
   "inert species", ":math:`\text{inert_species}`", "['H2O', target_species, background solutes] - [target_species]"
   "number of discretized operational time elements used for steady state approximation", ":math:`\text{ele_disc}`", "[0:elements_ss_approx]"
   "number of discretized trapezoidal area terms for steady state approximation", ":math:`\text{ele_index}`", "[1:elements_ss_approx]"

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
   "surface diffusion coefficient", ":math:`D_s`", "ds", "None", ":math:`\text{m}^2\text{/s}`"
   "liquid phase film transfer coefficient", ":math:`k_f`", "kf", "None", ":math:`\text{m/s}`"
   "equilibrium concentration of adsorbed phase with liquid phase", ":math:`q_e`", "equil_conc", "None", ":math:`\left( \text{kg}_\text{adsorbate}\text{/kg}_\text{adsorbent} \right)`"
   "solute distribution parameter", ":math:`D_g`", "dg", "None", ":math:`\text{dimensionless}`"
   "Biot number", ":math:`Bi`", "N_Bi", "None", ":math:`\text{dimensionless}`"
   "superficial velocity", ":math:`u_s`", "velocity_sup", "None", ":math:`\text{m/s}`"
   "interstitial velocity", ":math:`u_i`", "velocity_int", "None", ":math:`\text{m/s}`"
   "bed void fraction", ":math:`\epsilon`", "bed_voidage", "None", ":math:`\text{dimensionless}`"
   "bed length", ":math:`L`", "bed_length", "None", ":math:`\text{m}`"
   "bed diameter", ":math:`D`", "bed_diameter", "None", ":math:`\text{m}`"
   "bed area", ":math:`A`", "bed_area", "None", ":math:`\text{m}^2`"
   "bed volume", ":math:`V`", "bed_volume", "None", ":math:`\text{m}^3`"
   "empty bed contact time", ":math:`EBCT`", "ebct", "None", ":math:`\text{s}`"
   "fluid residence time in the bed", ":math:`\tau`", "residence_time", "None", ":math:`\text{s}`"
   "mass of fresh gac in the bed", ":math:`M_{GAC}`", "bed_mass_gac", "None", ":math:`\text{kg}`"
   "gac apparent density", ":math:`\rho_a`", "particle_dens_app", "None", ":math:`\text{kg/}\text{m}^3`"
   "gac bulk density", ":math:`\rho_b`", "particle_dens_bulk", "None", ":math:`\text{kg/}\text{m}^3`"
   "gac particle diameter", ":math:`d_p`", "particle_dia", "None", ":math:`\text{m}`"
   "minimum Stanton number to achieve a constant pattern solution", ":math:`St_{min}`", "min_N_St", "None", ":math:`\text{dimensionless}`"
   "minimum empty bed contact time to achieve a constant pattern solution", ":math:`EBCT_{min}`", "min_ebct", "None", ":math:`\text{s}`"
   "specific throughput from empirical equation", ":math:`T`", "throughput", "None", ":math:`\text{dimensionless}`"
   "minimum fluid residence time in the bed to achieve a constant pattern solution", ":math:`\tau_{min}`", "min_residence_time", "None", ":math:`\text{s}`"
   "minimum operational time of the bed from fresh to achieve a constant pattern solution", ":math:`t_{min}`", "min_operational_time", "None", ":math:`\text{s}`"
   "effluent to inlet concentration ratio at operational time", ":math:`\frac{C}{C_{0}}\bigg{|}_{z=L,/,t=t_{op}}`", "conc_ratio_replace", "None", ":math:`\text{dimensionless}`"
   "operational time of the bed from fresh", ":math:`t_{op}`", "operational_time", "None", ":math:`\text{s}`"
   "bed volumes treated at operational time", ":math:`BVT`", "bed_volumes_treated", "None", ":math:`\text{dimensionless}`"
   "specific throughput from empirical equation by discrete element", ":math:`T_{ele}`", "ele_throughput", "None", ":math:`x`"
   "minimum operational time of the bed from fresh to achieve a constant pattern solution by discrete element", ":math:`t_{min, ele}`", "ele_min_operational_time", "ele_index", ":math:`\text{s}`"
   "effluent to inlet concentration ratio at operational time by discrete element", ":math:`\left(\frac{C}{C_{0}}\right)_{ele}\bigg{|}_{z=L,/,t=t_{op_ e}}`", "ele_conc_ratio_replace", "ele_index", ":math:`\text{dimensionless}`"
   "operational time of the bed from fresh by discrete element", ":math:`t_{op, ele}`", "ele_operational_time", "ele_disc", ":math:`\text{s}`"
   "trapezoid rule of elements for numerical integration of average concentration ratio", ":math:`term_{ele}`", "ele_conc_ratio_avg", "ele_disc", ":math:`\text{dimensionless}`"
   "steady state approximation of average effluent to inlet concentration ratio in operational time by trapezoid rule", ":math:`\left(\frac{C}{C_{0}}\right)_{avg}`", "conc_ratio_avg", "None", ":math:`\text{dimensionless}`"
   "total mass of adsorbed species at operational time", ":math:`M`", "mass_adsorbed", "None", ":math:`\text{kg}`"
   "gac usage/replacement/regeneration rate", ":math:`\dot{m}_{GAC}`", "gac_usage_rate", "None", ":math:`\text{m/s}`"


The following variables are only built when specific configuration options are selected.

if ``cphsdm_calculation_method`` is set to ``input``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Stanton equation parameter 0", ":math:`a_0`", "a0", "None", ":math:`\text{dimensionless}`"
   "Stanton equation parameter 1", ":math:`a_1`", "a1", "None", ":math:`\text{dimensionless}`"
   "throughput equation parameter 0", ":math:`b_0`", "b0", "None", ":math:`\text{dimensionless}`"
   "throughput equation parameter 1", ":math:`b_1`", "b1", "None", ":math:`\text{dimensionless}`"
   "throughput equation parameter 2", ":math:`b_2`", "b2", "None", ":math:`\text{dimensionless}`"
   "throughput equation parameter 3", ":math:`b_3`", "b3", "None", ":math:`\text{dimensionless}`"
   "throughput equation parameter 4", ":math:`b_4`", "b4", "None", ":math:`\text{dimensionless}`"
   
if ``film_transfer_coefficient_type`` is set to ``calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Reynolds number", ":math:`Re`", "N_Re", "None", ":math:`\text{dimensionless}`"
   "Schmidt number", ":math:`Sc`", "N_Sc", "None", ":math:`\text{dimensionless}`"
   "shape correction factor", ":math:`SCF`", "shape_correction_factor", "None", ":math:`\text{dimensionless}`"

if ``surface_diffusion_coefficient_type`` is set to ``calculated``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "gac particle porosity", ":math:`\epsilon_p`", "particle_porosity", "None", ":math:`\text{dimensionless}`"
   "tortuosity of the path that the adsorbate must take as compared to the radius", ":math:`\tau_p`", "tort", "None", ":math:`\text{dimensionless}`"
   "surface-to-pore diffusion flux ratio", ":math:`S\!P\!D\!F\!R`", "spdfr", "None", ":math:`\text{dimensionless}`"

.. _GAC_equations:

Equations
-----------

.. csv-table::
   :header: "Description", "Equation"

   "equilibrium concentration", ":math:`q_e = kC_0^{1/n}`"
   "solute distribution parameter", ":math:`D_g=\frac{\rho_aq_e\left( 1-\epsilon \right)}{\epsilon C_0}`"
   "Biot number", ":math:`Bi=\frac{k_fd_p\left( 1-\epsilon \right)}{2D_sD_g\epsilon}`"
   "bed void fraction based on gac particle densities", ":math:`\epsilon=1-\frac{\rho_b}{\rho_a}`"
   "relating velocities based on bed voidage", ":math:`u_i=\frac{u_s}{\epsilon}`"
   "bed length based on velocity and ebct", ":math:`L=(EBCT)u_s`"
   "bed diameter and area relation", ":math:`A=\pi\left(\frac{D}{2}\right)^2`"
   "bed area based on velocity and volumetric flow", ":math:`A=\frac{Q}{u_s}`"
   "bed volume based on cylindrical dimensions", ":math:`V=AL`"
   "fluid residence time in the bed", ":math:`\tau=(EBCT)\epsilon`"
   "total mass of gac in the bed", ":math:`M_{GAC}=V\rho_b`"
   "minimum empty bed contact time to achieve constant pattern solution", ":math:`EBCT_{min}=\frac{St_{min}d_p}{2k_f\left( 1-\epsilon \right)}`"
   "minimum fluid residence time in the bed to achieve a constant pattern solution", ":math:`\tau_{min}=EBCT_{min}\epsilon`"
   "minimum operational time of the bed from fresh to achieve a constant pattern solution", ":math:`t_{min}=\tau_{min}\left( D_g+1 \right)T`"
   "elapsed operational time between a fresh bed and the theoretical bed replacement", ":math:`t_{op}=t_{min}+\left( \tau-\tau_{min} \right)\left( D_g+1 \right)`"
   "bed volumes treated", ":math:`BVT=\frac{t_{op}\epsilon}{\tau}`"
   "minimum operational time of the bed from fresh to achieve a constant pattern solution by discretized element", ":math:`t_{min, ele}=\tau_{min}\left( D_g+1 \right)T`"
   "creating evenly spaced discretized elements", ":math:`\frac{C}{C_{0}}\bigg{|}_{t=t_{op, ele}}=0.01+(ele-1)*\frac{\left(\frac{C}{C_{0}}\bigg{|}_{t=t_{op}}-0.01\right)}{num\text{_}ele}`"
   "finite element discretization of concentration ratios over time", ":math:`term_{ele}=\left(\frac{t_{op, ele}-t_{op, (ele-1)}}{t_{op}}\right)\frac{\left(\frac{C}{C_{0}}\bigg{|}_{t=t_{op, ele}}+{C_{0}}\bigg{|}_{t=t_{op, (ele-1)}}\right)}{2}`"
   "summation of finite elements for average concentration during operating time", ":math:`\left(\frac{C}{C_{0}}\right)_{avg}=\sum_{ele\text{_}index}term_{ele}`"
   "mass adsorbed in the operational time", ":math:`M=\frac{\dot{m}_{j}}{t_{op}}`"
   "steady state rate of new gac mass required", ":math:`\dot{m}_{GAC}=\frac{M_{GAC}}{t_{op}}`"


if ``cphsdm_calculation_method`` is set to ``input``:

.. csv-table::
   :header: "Description", "Equation"

   "minimum Stanton number to achieve constant pattern solution", ":math:`St_{min}=a_0Bi+a_1`"
   "throughput based on empirical 5-parameter regression", ":math:`T=b_0+b_1\left( \frac{C}{C_0} \right)^{b_2}+\frac{b_3}{1.01-\left( \frac{C}{C_0} \right)^{b_4}}`"
   "throughput based on empirical 5-parameter regression by discretized element", ":math:`T_{ele}=b_0+b_1\left(\frac{C}{C_{0}}\right)_{ele}^{b_2}+\frac{b_3}{1.01-\left(\frac{C}{C_{0}}\right)_{ele}^{b_4}}`"


if ``film_transfer_coefficient_type`` is set to ``calculated``:

.. csv-table::
   :header: "Description", "Equation"

   "Reynolds number calculation*", ":math:`Re=\frac{\rho_ld_pu_i}{\mu_l}`"
   "Schmidt number calculation*", ":math:`Sc=\frac{\mu_l}{\rho_sD_l}`"
   "liquid phase film transfer rate from the Gnielinshi correlation*", ":math:`k_f=(SCF)\frac{\left[ 1+1.5\left( 1-\epsilon \right) \right]D_l}{d_p}\left( 2+0.644Re^{\frac{1}{2}}Sc^{\frac{1}{3}} \right)`"

\*Subscript :math:`l` denotes bulk liquid phase properties, here those are supplied by the property package.


if ``surface_diffusion_coefficient_type`` is set to ``calculated``:

.. csv-table::
   :header: "Description", "Equation"

   "surface diffusion parameter (Crittenden, 1987)", ":math:`D_s=\left( S\!P\!D\!F\!R \right)\left( \frac{\epsilon_pC_0D_l}{\rho_aq_e\tau_p} \right)`"


.. _GAC_CPHSDM_params:

CPHSDM Empirical Parameters
--------------------------- 

If ``cphsdm_calculation_method`` is set to ``input``, the following parameters can be used for ``a0`` and ``a1`` to calculate the minimum Stanton number 
and ``b0``, ``b1``, ``b2``, ``b3``, and ``b4`` for the throughput equation. These tables were adapted directly from Hand, 1984.


.. _min_St_params:

.. table:: Empirical Parameters for Minimum Stanton Number Equation
   :align: left

   +---------------+-------------------------------+-----------------------+
   | 1/n           | 0.5 ≤ Bi ≤ 10                 | 10 ≤ Bi ≤ ∞           |
   +---------------+---------------+---------------+---------------+-------+
   |``freund_ninv``| ``a0``        | ``a1``        | ``a0``        | ``a1``|
   +===============+===============+===============+===============+=======+
   | 0.05          | 2.10526 × 10⁻²| 1.98947       | 0.22          | 0     |
   +---------------+---------------+---------------+---------------+-------+
   | 0.10          | 2.10526 × 10⁻²| 2.18947       | 0.24          | 0     |
   +---------------+---------------+---------------+---------------+-------+
   | 0.20          | 4.21053 × 10⁻²| 2.37985       | 0.28          | 0     |
   +---------------+---------------+---------------+---------------+-------+
   | 0.30          | 1.05263 × 10⁻¹| 2.54737       | 0.36          | 0     |
   +---------------+---------------+---------------+---------------+-------+
   | 0.40          | 2.31579 × 10⁻¹| 2.68421       | 0.50          | 0     |
   +---------------+---------------+---------------+---------------+-------+
   | 0.50          | 5.26316 × 10⁻¹| 2.73684       | 0.80          | 0     |
   +---------------+---------------+---------------+---------------+-------+
   | 0.60          | 1.15789 × 10⁰ | 3.42105       | 1.50          | 0     |
   +---------------+---------------+---------------+---------------+-------+
   | 0.70          | 1.78947 × 10⁰ | 7.10526       | 2.50          | 0     |
   +---------------+---------------+---------------+---------------+-------+
   | 0.80          | 3.68421 × 10⁰ | 13.1579       | 5.00          | 0     |
   +---------------+---------------+---------------+---------------+-------+
   | 0.90          | 6.31579 × 10⁰ | 56.8421       | 12.00         | 0     |
   +---------------+---------------+---------------+---------------+-------+

.. _throughout_params:

.. table:: Empirical Parameters for Throughput Equation
   :align: left

   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   |``freund_ninv``                | Bi        | ``b0``     | ``b1``     | ``b2``     | ``b3``     | ``b4``     |
   +===============================+===========+============+============+============+============+============+
   | 0.05                          | 0.5       | -5.447214  | 6.598598   | 0.026569   | 0.019384   | 20.450470  |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.05                          | 2.0       | -5.465811  | 6.592484   | 0.004989   | 0.004988   | 0.503520   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.05                          | 4.0       | -5.531155  | 6.584935   | 0.022380   | 0.009109   | 0.273705   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.05                          | 6.0       | -5.606508  | 6.582188   | 0.022088   | 0.013126   | 0.214246   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.05                          | 8.0       | -5.605003  | 6.504701   | 0.028072   | 0.017803   | 0.189537   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.05                          | 10.0      | -5.664173  | 6.456597   | 0.018157   | 0.019935   | 0.149314   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.05                          | 14.0      | -0.662780  | 1.411252   | 0.060079   | 0.020709   | 0.142393   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.05                          | 25.0      | -0.662783  | 1.350940   | 0.031007   | 0.020350   | 0.129998   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.05                          | ≥100.0    | -0.667873  | 0.713110   | 0.287309   | 0.016823   | 0.117271   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.10                          | 0.5       | -1.919873  | 3.055368   | 0.055488   | 0.024284   | 15.311766  |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.10                          | 2.0       | -2.278950  | 3.399925   | 0.046838   | 0.004751   | 0.384675   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.10                          | 4.0       | -2.338078  | 3.379926   | 0.143994   | 0.005308   | 0.248492   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.10                          | 6.0       | -2.407407  | 3.374131   | 0.041392   | 0.012552   | 0.229345   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.10                          | 8.0       | -2.477819  | 3.370954   | 0.038993   | 0.012257   | 0.194358   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.10                          | 10.0      | -2.566414  | 3.370950   | 0.035003   | 0.019386   | 0.150788   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.10                          | 16.0      | -2.567201  | 3.063341   | 0.020490   | 0.019483   | 0.136813   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.10                          | 25.0      | -2.568618  | 3.241783   | 0.009595   | 0.019962   | 0.121746   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.10                          | ≥100.0    | -2.568630  | 3.191482   | 0.015555   | 0.016781   | 0.101075   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.20                          | 0.5       | -1.441000  | 2.569000   | 0.069020   | 0.020333   | 0.211706   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.20                          | 2.0       | -1.474313  | 2.558000   | 0.058480   | 0.020536   | 0.217104   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.20                          | 4.0       | -1.506696  | 2.519259   | 0.055355   | 0.008797   | 0.182742   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.20                          | 8.0       | -1.533595  | 1.983018   | 0.069238   | 0.021302   | 0.188033   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.20                          | 20.0      | -0.161992  | 1.077521   | 0.144879   | 0.015500   | 0.168083   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.20                          | 25.0      | -1.409232  | 2.188339   | 0.152191   | 0.018142   | 0.156048   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.20                          | 13.0      | -1.369220  | 2.118545   | 0.039492   | 0.018453   | 0.127565   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.20                          | 25.0      | -1.514159  | 2.209450   | 0.017937   | 0.018150   | 0.107466   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.20                          | ≥100.0    | -0.680346  | 0.649006   | 2.570086   | 0.014947   | 0.369818   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.30                          | 0.5       | -1.758696  | 2.846576   | 0.049530   | 0.003022   | 0.131927   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.30                          | 2.0       | -1.657826  | 2.688895   | 0.048490   | 0.005216   | 0.139368   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.30                          | 4.0       | -0.565664  | 1.537383   | 0.084351   | 0.008088   | 0.139906   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.30                          | 8.0       | -0.197077  | 1.118564   | 0.117894   | 0.011527   | 0.137546   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.30                          | 10.0      | -0.197070  | 1.069216   | 0.117690   | 0.013925   | 0.133691   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.30                          | 15.0      | -0.173358  | 1.000001   | 0.120311   | 0.015490   | 0.124052   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.30                          | 20.0      | -0.173350  | 0.919141   | 0.071760   | 0.014546   | 0.085279   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.30                          | ≥100.0    | 0.696617   | 0.516957   | 2.054587   | 0.012961   | 0.333578   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.40                          | 0.5       | -0.534251  | 1.603834   | 0.095492   | 0.014624   | 0.212861   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.40                          | 2.0       | -0.166270  | 1.190897   | 0.122280   | 0.006261   | 0.134278   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.40                          | 4.0       | -0.166273  | 1.131946   | 0.115513   | 0.008634   | 0.136997   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.40                          | 6.0       | -0.166270  | 1.089783   | 0.112284   | 0.010645   | 0.141626   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.40                          | 8.0       | -0.491219  | 0.491833   | 0.487414   | 0.013717   | 0.144115   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.40                          | 12.0      | -0.564119  | 0.419196   | 0.639819   | 0.011543   | 0.149005   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.40                          | 20.0      | -0.640669  | 0.432466   | 1.048506   | 0.009892   | 0.166825   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.40                          | 25.0      | -0.672533  | 0.672533   | 1.153169   | 0.011280   | 0.212683   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.40                          | ≥100.0    | 0.741435   | 0.448054   | 1.929879   | 0.010152   | 0.306448   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.50                          | 0.5       | -0.048000  | 1.099652   | 0.158995   | 0.005467   | 0.139116   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.50                          | 4.0       | -0.048000  | 0.982757   | 0.111618   | 0.008072   | 0.111404   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.50                          | 10.0      | 0.094602   | 0.754878   | 0.092069   | 0.009877   | 0.090763   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.50                          | 14.0      | 0.023000   | 0.802068   | 0.057545   | 0.009662   | 0.084532   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.50                          | 25.0      | 0.023000   | 0.793673   | 0.039324   | 0.009326   | 0.082751   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.50                          | ≥100.0    | 0.529213   | 0.291801   | 0.082428   | 0.008317   | 0.075461   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.60                          | 0.5       | 0.352536   | 0.692114   | 0.263134   | 0.005482   | 0.121775   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.60                          | 2.0       | 0.521979   | 0.504220   | 0.327290   | 0.005612   | 0.128679   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.60                          | 6.0       | 0.676253   | 0.334583   | 0.482297   | 0.005898   | 0.138946   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.60                          | 14.0      | 0.769531   | 0.259497   | 0.774068   | 0.005600   | 0.165513   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.60                          | 50.0      | 0.849057   | 0.215799   | 1.343183   | 0.004725   | 0.223759   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.60                          | ≥100.0    | 0.831231   | 0.227304   | 1.174756   | 0.004961   | 0.212109   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.60                          | 0.5       | 0.575024   | 0.449062   | 0.278452   | 0.004122   | 0.121682   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.60                          | 4.0       | 0.715269   | 0.307172   | 0.442104   | 0.004371   | 0.138351   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.70                          | 12.0      | 0.787940   | 0.243548   | 0.661599   | 0.004403   | 0.162595   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.70                          | 25.0      | 0.829492   | 0.204078   | 0.784529   | 0.004050   | 0.179005   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.70                          | ≥100.0    | 0.847012   | 0.190678   | 0.931686   | 0.003849   | 0.183239   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.80                          | 0.5       | 0.708905   | 0.314101   | 0.357499   | 0.003276   | 0.119300   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.60                          | 4.0       | 0.715269   | 0.307172   | 0.442104   | 0.004371   | 0.138351   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.70                          | 12.0      | 0.787940   | 0.243548   | 0.661599   | 0.004403   | 0.162595   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.70                          | 25.0      | 0.829492   | 0.204078   | 0.784529   | 0.004050   | 0.179005   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.70                          | ≥100.0    | 0.847012   | 0.190678   | 0.931686   | 0.003849   | 0.183239   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.80                          | 0.5       | 0.708905   | 0.314101   | 0.357499   | 0.003276   | 0.119300   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.80                          | 4.0       | 0.784576   | 0.239663   | 0.484422   | 0.003206   | 0.134987   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.80                          | 14.0      | 0.839439   | 0.188966   | 0.648124   | 0.003306   | 0.157697   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.80                          | ≥100.0    | 0.882747   | 0.146229   | 0.807987   | 0.002537   | 0.174543   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.90                          | 0.5       | 0.865453   | 0.157618   | 0.444973   | 0.001650   | 0.148084   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.90                          | 4.0       | 0.854768   | 0.171434   | 0.495042   | 0.001910   | 0.142251   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.90                          | 16.0      | 0.866180   | 0.163992   | 0.573946   | 0.001987   | 0.157594   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+
   | 0.90                          | ≥100.0    | 0.893192   | 0.133039   | 0.624100   | 0.001740   | 0.164248   |
   +-------------------------------+-----------+------------+------------+------------+------------+------------+


Code Documentation
-------------------

* :mod:`watertap.unit_models.gac`
* :mod:`watertap.costing.unit_models.gac`

References
-----------
Hand, D. W., Crittenden, J. C., & Thacker, W. E. (1984). Simplified models for design of fixed-bed adsorption systems.
Journal of Environmental Engineering, 110(2), 440-456.

Crittenden, J., Rhodes, R., Hand, D., Howe, K., & Tchobanoglous, G. (2012). MWHs Water Treatment. Principles and Design.
John Wiley & Sons.

Crittenden, J. C., Berrigan, J. K., Hand, D. W., & Lykins, B. (1987). Design of Rapid Fixed‐Bed Adsorption Tests for
Nonconstant Diffusivities. Journal of Environmental Engineering, 113(2), 243–259.
