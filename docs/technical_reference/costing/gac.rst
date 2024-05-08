Granular Activated Carbon Costing Method
=========================================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., `m.fs.costing.gac_pressure` or `m.fs.costing.gac_gravity`) when applying the `cost_gac` costing
method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "Number of contactors operating in parallel", ":math:`N_{op}`", "num_contactors_op", "1", ":math:`\text{dimensionless}`"
   "Number of redundant contactors in parallel", ":math:`N_{red}`", "num_contactors_redundant", "1", ":math:`\text{dimensionless}`"
   "Fraction of spent GAC adsorbent to be regenerated and reused", ":math:`f_{regen}`", "regen_frac", "0.70", ":math:`\text{dimensionless}`"
   "Reference maximum value of GAC initial charge mass where economy of scale no longer discounts the unit price (U.S. EPA, 2021)", ":math:`M_{GAC}^{ref}`", "bed_mass_max_ref", "18143.7", ":math:`kg`"
   "Contactor polynomial cost coefficients", ":math:`x`", "contactor_cost_coeff", "tabulated", ":math:`\text{dimensionless}`"
   "Adsorbent exponential cost coefficients", ":math:`y`", "adsorbent_unit_cost_coeff", "tabulated", ":math:`\text{dimensionless}`"
   "Other process costs power law coefficients", ":math:`z`", "other_cost_param", "tabulated", ":math:`\text{dimensionless}`"
   "Unit cost to regenerate spent GAC adsorbent", ":math:`C_{unit,regen}`", "regen_unit_cost", "4.28352", ":math:`USD\_2020/kg`"
   "Unit cost to makeup spent GAC adsorbent with fresh adsorbent", ":math:`C_{unit,makeup}`", "makeup_unit_cost", "4.58223", ":math:`USD\_2020/kg`"
   "Energy consumption polynomial coefficients", ":math:`alpha`", "energy_consumption_coeff", "tabulated", ":math:`\text{dimensionless}`"

Parameters which are tabulated have costing methods available for either steel pressure vessel contactors (default) or concrete gravity basin contactors. Given that the form of the costing
component equations are different (polynomial, exponential, and power law), the units associated with the parameters are embedded in the constraints and not directly applied to the variable.
Additionally, the index is generalized to its position ``([0:len(parameter_data)])`` in the list, although some parameters are coefficients while others are exponents (see equations below for details).
Variables with the (U.S. EPA, 2021) citation are directly taken from previously determined expressions. Other variables are regressed from higher detailed costing methods in (U.S. EPA, 2021).

.. csv-table::
   :header: "Variable Name", "Contactor Type", "Index 0", "Index 1", "Index 2", "Index 3"

   "adsorbent_unit_cost_coeff (U.S. EPA, 2021)", "n/a", "4.58342", "-1.25311e-5", "", ""
   "contactor_cost_coeff (U.S. EPA, 2021)", "pressure", "10010.9", "2204.95", "-15.9378", "0.110592"
   "contactor_cost_coeff", "gravity", "75131.3", "735.550", "-1.01827", "0"
   "other_cost_param", "pressure", "16660.7", "0.552207", "", ""
   "other_cost_param", "gravity", "38846.9", "0.490571", "", ""
   "energy_consumption_coeff_data", "pressure", "8.09926e-4", "8.70577e-4", "0", ""
   "energy_consumption_coeff_data", "gravity", "0.123782", "0.132403", "-1.41512e-5", ""

\*Energy consumption is the sum of energy required to operate booster, backwash, and residual pumps.

Costing Method Variables
++++++++++++++++++++++++

The following variables are constructed on the unit block (e.g., m.fs.unit.costing) when applying the `cost_gac` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Default Value", "Units"

   "Total unit capital cost", ":math:`C_{cap}`", "capital_cost", "10,000", ":math:`USD\_2020`"
   "Capital cost of contactor(s)", ":math:`C_{cap,bed}`", "contactor_cost", "10,000", ":math:`USD\_2020`"
   "Mass of GAC initial charge used to determine the unit cost", ":math:`M_{GAC}^{min}`", "bed_mass_gac_ref", "4", ":math:`kg`"
   "Unit cost of GAC adsorbent for initial charge", ":math:`C_{unit,carbon}`", "adsorbent_unit_cost", "2", ":math:`USD\_2020/kg`"
   "Capital cost of GAC adsorbent", ":math:`C_{cap,carbon}`", "adsorbent_cost", "10,000", ":math:`USD\_2020`"
   "Capital costs of other process supplements", ":math:`C_{cap,other}`", "other_process_cost", "10,000", ":math:`USD\_2020`"
   "Approximate GAC system energy consumption*", ":math:`P`", "energy_consumption", "100", ":math:`kW`"
   "Fixed operating costs", ":math:`C_{op}`", "fixed_operating_cost", "10,000", ":math:`USD\_2020/yr`"
   "Operating costs to regenerate spent GAC adsorbent", ":math:`C_{op,regen}`", "gac_regen_cost", "10,000", ":math:`USD\_2020/yr`"
   "Operating costs to makeup spent GAC adsorbent with fresh adsorbent", ":math:`C_{op,makeup}`", "gac_makeup_cost", "10,000", ":math:`USD\_2020/yr`"

Capital Cost Calculations
+++++++++++++++++++++++++

Costing GAC contactors is defaulted to purchasing 1 operational and 1 redundant contactor for alternating operation with minimal downtown. For large systems, this may be a poor
assumption considering vessel sizing and achieving pseudo-steady state. The number of contactors input by the user should justify reasonable (commercially available) dimensions
of identical modular contactors in parallel. When costing several operational vessels, the area reported in the unit model should be interpreted as the sum of the areas across
all operating GAC contactors. The costing parameters may be selected from either steel pressure-fed vessels or concrete gravity-fed basins by the ``contactor_type`` argument.
Note this only affects costing calculations. Volume dimensions calculations within the model remain assuming a cylindrical bed. Capital costs are determined by the summation of
three costing terms. Each term is is calculated by a one parameter (different for each term) function considering economy of scale.

    .. math::

        C_{cap,tot} = C_{cap,bed}+C_{cap,carbon}+C_{cap,other}

Contactor and GAC adsorbent capital costs are estimated using functions and parameters reported in US EPA, 2021. Contactors are assumed to be carbon steel pressure vessels with
plastic internals and are determined as a polynomial function of individual contactor volume. The unit cost per kilogram of GAC adsorbent needed is calculated using an exponential
function. A maximum reference mass is imposed in the costing method to define a best available price where above this required charge, the price would no longer be discounted.
Other process costs (vessels, pipes, instrumentation, and controls) included in the US EPA, 2021 model are aggregated into a separate term. The parameters for the power law function
with respect to the total system contactor volume were regressed using results from the US EPA, 2021 model.

    .. math::

        & C_{cap,bed} = \left( N_{op}+N_{red} \right)\left( x_0+x_1\left( \frac{V}{N_{op}} \right)+x_2\left( \frac{V}{N_{op}} \right)^2+x_3\left( \frac{V}{N_{op}} \right)^3 \right) \\\\
        & M_{GAC}^{min} = \text{min}\left(M_{GAC}^{model}, M_{GAC}^{ref}\right) \\\\
        & C_{carbon} = y_0e^{y_1M_{GAC}^{min}} \\\\
        & C_{cap,carbon} = C_{carbon}M_{GAC} \\\\
        & C_{cap,other} = z_0\left( \left( N_{op}+N_{red} \right)\frac{V}{N_{op}} \right)^{z_1}


Note that given the the ability to alter the parameters in these correlations, GAC adsorbent unit costs (:math:`C_{carbon}`) may be fixed to a value (:math:`y_0`) by setting :math:`y_1=0`.

Operating Cost Calculations
+++++++++++++++++++++++++++

Operating costs are calculated as the cost to replace spent GAC adsorbent in the contactor beds. Energy is costed as a flow term by the WaterTAP costing method.

    .. math::

        C_{op,tot} = C_{op,regen}+C_{op,makeup}

Since the replacement adsorbent purchases are expected to be purchased in bulk at smaller quantities than the initial charge, the cost of fresh GAC adsorbent for makeup has an different
cost per unit mass, expected to be higher than the initial charge unit cost.

    .. math::

        & C_{op,regen} = f_{regen}C_{unit,regen}\dot{m}_{GAC}^{model} \\\\
        & C_{op,makeup} = \left( 1-f_{regen} \right)C_{unit,makeup}\dot{m}_{GAC}^{model} \\\\
        & P = \alpha_0+\alpha_1V+\alpha_2V^2
 
Code Documentation
------------------

* :mod:`watertap.costing.unit_models.gac`

References
----------
United States Environmental Protection Agency. (2021). Work Breakdown Structure-Based Cost Model for Granular Activated
Carbon Drinking Water Treatment.