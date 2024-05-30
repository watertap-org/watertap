Ion Exchange
============

Introduction
------------

The simple ion exchange (IX) flowsheet can be simulated to predict the performance of an IX system to remove targeted ions and components. This flowsheet can
be useful to expedite the set-up, usage, and costing of an IX system for conventional water treatment applications using the Langmuir isotherm and
the constant pattern assumption.

Implementation
--------------

The ion exchange flowsheet demonstration proceeds through four steps:

    1. ``ix_build``: This function builds the flowsheet with a list of ions as input. If the keyword argument ``target_ion`` is not 
    used, the first ion in the list of ions provided is used as the ``target_ion`` configuration argument for the ``IonExchange0D`` model.
    The ion used in the demonstration is calcium (``Ca_2+``), but the local function ``get_ion_config`` can be used to get diffusivity, molecular weight, 
    and charge data for sodium (``Na_+``), chloride (``Cl_-``), magnesium (``Mg_2+``), and sulfate (``SO4_2-``). 

    2. ``



Only consisting of a single unit operation, the assumptions for the flowsheet are aligned with those detailed in the :doc:`ion exchange unit model documentation </technical_reference/unit_models/ion_exchange>`.
The code-based naming of modeling objects for the inlets, outlets, units, and streams are shown in Figure 1.

.. .. figure:: ../../_static/flowsheets/gac.png
..     :width: 500
..     :align: center

..     Figure 1. GAC flowsheet

The following modeling components are used within the flowsheet:

Documentation for property models:
    * :doc:`/technical_reference/property_models/mc_aq_sol`
Documentation for unit models:
    * :doc:`/technical_reference/unit_models/ion_exchange_0D`
Documentation for unit models from IDAES:
    * :doc:`idaes:reference_guides/model_libraries/generic/unit_models/feed`
    * :doc:`idaes:reference_guides/model_libraries/generic/unit_models/product`
Documentation for costing models:
    * :doc:`/technical_reference/costing/watertap_costing`
    * :doc:`/technical_reference/costing/ion_exchange`

Degrees of Freedom
------------------

The degrees of freedom for the flowsheet can change depending on the configuration options specified during the build. Excluding those variables which are
only necessary for specific configuration options, the following variables are initially fixed for simulating the GAC flowsheet (i.e., degrees of freedom = 0):

    * feed conditions (component flows, temperature, pressure)
    * 

Flowsheet Specifications
------------------------

.. csv-table::
   :header: "Description", "Value", "Units"

   "Feed molar flowrate of water", "2777.5", ":math:`\text{mol}/\text{s}`"
   "Feed molar flowrate of target ion", "0.125", ":math:`\text{mol}/\text{s}`"
   "feed temperature", "298.15", ":math:`\text{K}`"
   "feed pressure", "101325", ":math:`\text{Pa}`"
   "Freundlich isotherm k parameter", "10", ":math:`\left(\text{m}^3\text{/kg}\right)^\left( \frac{1}{n} \right)`"
   "Freundlich isotherm 1/n parameter", "0.9", ":math:`\text{dimensionless}`"
   "liquid phase film transfer coefficient", "5e-5", ":math:`\text{m/s}`"
   "surface diffusion coefficient", "2e-13", ":math:`\text{m}^2\text{/s}`"
   "gac apparent density", "750", ":math:`\text{kg/}\text{m}^3`"
   "gac particle diameter", "0.001", ":math:`\text{m}`"
   "empty bed contact time", "600", ":math:`\text{s}`"
   "bed void fraction", "0.4", ":math:`\text{dimensionless}`"
   "bed length", "6", ":math:`\text{m}`"
   "effluent to inlet concentration ratio at operational time", "0.50", ":math:`\text{dimensionless}`"
   "Stanton equation parameter 0", "3.68421", ":math:`\text{dimensionless}`"
   "Stanton equation parameter 1", "13.1579", ":math:`\text{dimensionless}`"
   "throughput equation parameter 0", "0.784576", ":math:`\text{dimensionless}`"
   "throughput equation parameter 1", "0.239663", ":math:`\text{dimensionless}`"
   "throughput equation parameter 2", "0.484422", ":math:`\text{dimensionless}`"
   "throughput equation parameter 3", "0.003206", ":math:`\text{dimensionless}`"
   "throughput equation parameter 4", "0.134987", ":math:`\text{dimensionless}`"

Future Refinements
------------------

The following modifications to the GAC flowsheet are planned for development:

    * Add surrogate models to lessen the need for numerous empirical parameters
    * Improve auto-scaling of model for ease of use

Code Documentation
------------------

* :mod:`watertap.examples.flowsheets.gac`

References
----------
Hand, D. W., Crittenden, J. C., & Thacker, W. E. (1984). Simplified models for design of fixed-bed adsorption systems.
Journal of Environmental Engineering, 110(2), 440-456.

Crittenden, J., Rhodes, R., Hand, D., Howe, K., & Tchobanoglous, G. (2012). MWHs Water Treatment. Principles and Design.
John Wiley & Sons.

United States Environmental Protection Agency. (2021). Work Breakdown Structure-Based Cost Model for Granular Activated
Carbon Drinking Water Treatment.