
Flexible Reverse Osmosis Flowsheet
============================

Introduction
------------

The flowsheets represents the operational window and energy consumption of two different RO plants.
The flowsheet is used with the Pricetaker model to determine the cost-optimal operation based on the treatment 
system energy surrogates and constraints as well as the variable grid electricity costs. 


File Structure
---------------

The flowsheet is 'wrd_ro_flowsheet' in watertap.flowsheets.flex_desal.wrd_ro_flowsheet.
The parameters of each unit model are defined in params. They are passed during the model build.
The specific unit models for the UF and RO are in watertap.flowsheets.flex_desal.wrd_unit_models. 
The general unit model are in watertap.flowsheets.flex_desal.unit_models. 
A few additional function are included in the watertap.flowsheets.flex_desal.utils

Parameters
----------

The top-level ``FlexDesalParams`` values include:

.. csv-table::
  :header: "Parameter", "Default value", "Units", "Description"

  "``start_date``", "``2022-07-05 00:00:00``", "timestamp", "Start of the simulation horizon."
  "``end_date``", "``2022-07-06 00:00:00``", "timestamp", "End of the simulation horizon."
  "``timestep_hours``", "``0.25``", "h", "Length of each simulation time step."
  "``product_water_price``", "``0``", "$ / m3", "Unit revenue for produced water."
  "``fixed_monthly_cost``", "``766000``", "$ / month", "Fixed monthly customer/facility charge."
  "``customer_rate``", "``100``", "-", "Customer-rate multiplier used in tariff calculations."
  "``constrain_to_baseline_production``", "``False``", "-", "Whether to enforce baseline production tracking."
  "``curtailment_fraction``", "``0.0``", "dimensionless", "Allowed fractional curtailment relative to baseline production."
  "``annual_production_AF``", "``3125``", "acre-ft/year", "Annual production target used for absolute production constraints."
  "``production_constraint_to_objective``", "``False``", "-", "Whether production compliance is enforced through objective penalization."
  "``production_constraint_penalty``", "``0.6``", "$ / m3 (effective penalty scale)", "Penalty weight applied when production target is incorporated in the objective."
  "``emissions_cost``", "``0``", "$ / kg", "Cost assigned to emissions associated with grid electricity use."
  "``include_demand_response``", "``False``", "-", "Enable demand-response price/revenue terms."
  "``include_battery``", "``False``", "-", "Enable battery operation model."
  "``include_onsite_solar``", "``False``", "-", "Enable onsite solar generation model."
  "``onsite_capacity``", "``0``", "kW", "Installed onsite generation capacity."
  "``nonworking_hours``", "``[]``", "hour index list", "Hours where startup/shutdown actions may be restricted."
  "``rainy_days``", "``None``", "day count", "Optional count/indicator used by scenario-specific logic."
  "``CAPEX_yr``", "``None``", "$ / year", "Optional annualized CAPEX value for economic reporting."
  "``max_daily_shutdowns``", "``None``", "count/day", "Optional limit on shutdown events over a daily rolling window."



Unit Models
-----------
The existing flowsheets are built with the following unit models:

Intake
~~~~~~~~

The intake unit uses the ``IntakeParams`` dataclass in
``watertap.flowsheets.flex_desal.params``.

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default value
     - Units
     - Description
   * - ``energy_intensity``
     - ``0.157121734``
     - kWh/m3
     - Specific intake energy intensity.
   * - ``minimum_flowrate``
     - ``1063.5``
     - m3/h
     - Minimum intake flowrate.
   * - ``nominal_flowrate``
     - ``1063.5``
     - m3/h
     - Nominal intake flowrate.
   * - ``maximum_flowrate``
     - ``1063.5``
     - m3/h
     - Maximum intake flowrate.
   * - ``feed_cost``
     - ``None``
     - $/m3
     - Optional feed-water cost.
   * - ``chemical_cost``
     - ``None``
     - $/m3
     - Optional chemical cost.

Pretreatment (UF)
~~~~~~~~~~~~~~~~~

The UF unit uses the ``WRD_UFParams`` dataclass in
``watertap.flowsheets.flex_desal.params``.

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default value
     - Units
     - Description
   * - ``num_uf_pumps``
     - ``4``
     - None
     - Number of UF pumps represented in the model.
   * - ``minimum_operating_pumps``
     - ``1``
     - None
     - Minimum number of UF pumps that must operate.
   * - ``allow_shutdown``
     - ``True``
     - None
     - Enables UF unit on/off commitment logic.
   * - ``minimum_flowrate``
     - ``344``
     - m3/h
     - Minimum UF pump flowrate when operating (m3/h).
   * - ``nominal_flowrate``
     - ``900``
     - m3/h
     - Nominal UF pump flowrate (m3/h).
   * - ``maximum_flowrate``
     - ``989``
     - m3/h
     - Maximum UF pump flowrate (m3/h).
   * - ``nominal_recovery``
     - ``1``
     - dimensionless
     - Nominal UF recovery.
   * - ``minimum_uptime``
     - ``2``
     - time steps
     - Minimum number of time steps the UF pump stays on after startup.
   * - ``minimum_downtime``
     - ``2``
     - time steps
     - Minimum number of time steps the UF pump stays off after shutdown.
   * - ``startup_delay``
     - ``1``
     - time steps
     - Delay (time steps) between startup command and operation.
   * - ``allow_variable_recovery``
     - ``False``
     - None
     - Recovery is fixed at nominal value for this WRD UF model.
   * - ``surrogate_type``
     - ``quadratic_energy_intensity``
     - None
     - UF energy intensity surrogate form.
   * - ``surrogate_a``, ``surrogate_b``, ``surrogate_c``
     - ``1``, ``1``, ``1``
     - surrogate-dependent
     - Coefficients for the quadratic UF energy intensity surrogate.


RO
~~

The RO unit uses the ``WRD_ROParams`` dataclass in
``watertap.flowsheets.flex_desal.params``.

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default value
     - Units
     - Description
   * - ``num_ro_skids``
     - ``4``
     - None
     - Number of RO skids represented in the model.
   * - ``minimum_operating_skids``
     - ``2``
     - None
     - Minimum number of RO skids that must operate.
   * - ``allow_shutdown``
     - ``True``
     - None
     - Enables RO skid on/off commitment logic.
   * - ``minimum_flowrate``
     - ``0``
     - m3/h
     - Minimum RO skid flowrate when operating (m3/h).
   * - ``nominal_flowrate``
     - ``337.670``
     - m3/h
     - Nominal RO skid flowrate (m3/h).
   * - ``maximum_flowrate``
     - ``400``
     - m3/h
     - Maximum RO skid flowrate (m3/h).
   * - ``minimum_recovery``
     - ``0.88``
     - dimensionless
     - Minimum RO recovery.
   * - ``nominal_recovery``
     - ``0.92``
     - dimensionless
     - Nominal RO recovery.
   * - ``maximum_recovery``
     - ``0.925``
     - dimensionless
     - Maximum RO recovery.
   * - ``minimum_uptime``
     - ``2``
     - time steps
     - Minimum number of time steps a skid stays on after startup.
   * - ``minimum_downtime``
     - ``2``
     - time steps
     - Minimum number of time steps a skid stays off after shutdown.
   * - ``startup_delay``
     - ``1``
     - time steps
     - Delay (time steps) between startup command and operation.
   * - ``allow_variable_recovery``
     - ``False``
     - None
     - Recovery is bounded but not optimized as a free variable in default setup.
   * - ``surrogate_type``
     - ``constant_energy_intensity``
     - None
     - RO surrogate type selector for the WRD case.
   * - ``surrogate_file``
     - ``None``
     - path
     - Optional file path for a loaded RO energy surrogate.



Posttreatment (UV)
~~~~~~~~~~~~~~~~~~

The UV posttreatment unit uses the ``PosttreatmentParams`` dataclass in
``watertap.flowsheets.flex_desal.params``.

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default value
     - Units
     - Description
   * - ``energy_intensity``
     - ``0.41``
     - kWh/m3
     - Specific UV energy intensity.
   * - ``leakage_fraction``
     - ``0``
     - dimensionless
     - Fraction of inlet flow not recovered in UV posttreatment.
   * - ``chemical_cost``
     - ``None``
     - $/m3
     - Optional variable chemical cost ($/m3).

Brine Discharge
~~~~~~~~~~~~~~~~

The brine discharge unit uses the ``BrineDischargeParams`` dataclass in
``watertap.flowsheets.flex_desal.params``.

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default value
     - Units
     - Description
   * - ``energy_intensity``
     - ``0.1``
     - kWh/m3
     - Specific brine-discharge energy intensity.
   * - ``brine_cost``
     - ``None``
     - $/m3
     - Optional brine-disposal cost.


Common Variables
-----------------

There are several variables common across the unit models.

.. csv-table::
  :header: "Variable", "Symbol used in this document", "Used in equations"

  "Flowrate", ":math:`Q` (e.g., :math:`Q^{RO}_{t,i}`, :math:`Q^{UF}_{t,i}`)", "RO/UF positive and negative flow-change detection"
  "Recovery", ":math:`R`", "None in current Equations and Relationships table"
  "Operating mode", ":math:`u` (e.g., :math:`u^{brine}_{t}`)", "Total feed cost"
  "Startup", ":math:`SU`", "None in current Equations and Relationships table"
  "Shutdown", ":math:`s`", "Degree of flex"
  "Energy intensity", ":math:`EI`", "None in current Equations and Relationships table"
  "Flow cost", ":math:`c` (e.g., :math:`c^{intake}_{t}`, :math:`c^{brine}_{t}`)", "Total feed cost, Total brine-discharge cost, Total chemical cost"

Common Equations
----------------

The following equations are built for every unit model via ``_add_required_variables``
in ``watertap.flowsheets.flex_desal.unit_models``.

.. csv-table::
   :header: "Description", "Equation"

   "Mass balance", ":math:`Q^{feed} = Q^{product} + Q^{reject}`"
   "Product flowrate", ":math:`Q^{product} = Q^{feed} \cdot R`"
   "Power consumption", ":math:`P = EI \cdot Q^{product}`"

Symmetry Breaking Equations
---------------------------


Optional Helper Functions
-----------
In addition to the functions described in the first tutorial, the WRD implementation
uses several helper functions in ``watertap.flowsheets.flex_desal.wrd_ro_flowsheet``.
The following quantities are only built when specific helper functions are used.

* ``add_flow_costs(m)``
  Adds total feed, brine-discharge, and chemical-cost expressions over the time horizon.

  .. csv-table::
     :header: "Description", "Symbol", "Name", "Index", "Units"

     "Total feed cost", ":math:`C_{feed}`", "total_feed_cost", "None", ":math:`\$`"
     "Total brine-discharge cost", ":math:`C_{brine}`", "total_brine_cost", "None", ":math:`\$`"
     "Total chemical cost", ":math:`C_{chem}`", "total_chemical_cost", "None", ":math:`\$`"

* ``add_flow_changes_penalty(m)``
  Introduces binary variables and Big-M constraints to detect RO/UF flowrate changes
  between time steps, then adds a penalty expression for frequent changes.

  .. csv-table::
     :header: "Description", "Symbol", "Name", "Index", "Units"

      "RO flow-change indicator", ":math:`y^{RO}`", "flow_changed", "[t, i]", ":math:`\text{dimensionless}`"
      "UF flow-change indicator", ":math:`y^{UF}`", "uf_flow_changed", "[t, i]", ":math:`\text{dimensionless}`"
     "Total flow-change penalty", ":math:`C_{chg}`", "flow_changes_penalty", "None", ":math:`\$`"


* ``calculate_replacement_costs(m)``
  Computes a flexibility metric from shutdown behavior and uses that metric to form
  annualized replacement-cost expressions for configured replacement categories.

  .. csv-table::
     :header: "Description", "Symbol", "Name", "Index", "Units"

     "Degree of flex", ":math:`f`", "degree_of_flex", "None", ":math:`\text{dimensionless}`"
     "Total replacement cost", ":math:`C_{rep}`", "total_replacement_cost", "None", ":math:`\$`"

* ``calculate_flexibility_metrics(m, baseline_power, baseline_electricity_cost, baseline_replacement_cost)``
  Post-processes solved results to estimate flexibility quantities such as charge/discharge
  capacities and levelized value/cost of flexibility. These metrics are defined in Rao et al. [1].

  .. csv-table::
     :header: "Description", "Symbol", "Name", "Index", "Units"

     "Maximum power draw", ":math:`P_{max}`", "maximum_power", "None", ":math:`\text{kW}`"
     "Discharge energy capacity", ":math:`E_{dis}`", "energy_capacity", "None", ":math:`\text{kWh}`"
     "Discharge power capacity", ":math:`P_{dis}`", "power_capacity", "None", ":math:`\text{kW}`"
     "Levelized value/cost of flexibility", ":math:`LVOF`", "LVOF", "None", ":math:`\$/\text{kWh}`"

* ``begin_and_end_constraint(m)``
  Enforces cyclic operation by matching RO train 1 operating mode at the first and
  last time points.

* ``add_working_hours_constraint(m)``
  Adds constraints that prevent RO startup and shutdown events during configured
  nonworking hours.

* ``restrict_flexible_trains(m, num_flexible_trains)``
  Restricts flexibility to a selected subset of RO skids by fixing startup and
  shutdown decisions to zero for non-flexible skids.


Equations and Relationships
---------------------------

Each equation below is cross-referenced to the helper-function section where its
left-hand-side quantity is introduced.

.. csv-table::
  :header: "Description", "Defined in", "Equation"

  "Total feed cost", "``add_flow_costs(m)``", ":math:`C_{feed} = \Delta t \sum_{t} c^{intake}_{t}(1 - u^{brine}_{t})`"
  "Total brine-discharge cost", "``add_flow_costs(m)``", ":math:`C_{brine} = \Delta t \sum_{t} c^{brine}_{t}`"
  "Total chemical cost", "``add_flow_costs(m)``", ":math:`C_{chem} = \Delta t \left(\sum_{t} c^{intake,chem}_{t} + \sum_{t} c^{post,chem}_{t}\right)`"
  "RO positive flow-change detection", "``add_flow_changes_penalty(m)``", ":math:`M^{RO} y^{RO}_{t,i} \geq Q^{RO}_{t,i} - Q^{RO}_{t-1,i}`"
  "RO negative flow-change detection", "``add_flow_changes_penalty(m)``", ":math:`M^{RO} y^{RO}_{t,i} \geq Q^{RO}_{t-1,i} - Q^{RO}_{t,i}`"
  "UF positive flow-change detection", "``add_flow_changes_penalty(m)``", ":math:`M^{UF} y^{UF}_{t,i} \geq Q^{UF}_{t,i} - Q^{UF}_{t-1,i}`"
  "UF negative flow-change detection", "``add_flow_changes_penalty(m)``", ":math:`M^{UF} y^{UF}_{t,i} \geq Q^{UF}_{t-1,i} - Q^{UF}_{t,i}`"
  "Total flow-change penalty", "``add_flow_changes_penalty(m)``", ":math:`C_{chg} = 50\left(\sum_{t,i} y^{RO}_{t,i} + \sum_{t,i} y^{UF}_{t,i}\right)`"
  "Degree of flex", "``calculate_replacement_costs(m)``", ":math:`f = \frac{\sum_{t,i} s_{t,i}}{2 N_{days} N_{RO}}`"
  "Total replacement cost", "``calculate_replacement_costs(m)``", ":math:`C_{rep} = \sum_{k} \frac{C_{rep,k}}{L_k \left(1 - \phi_k f\right)} \frac{N_{months}}{12}`"
  "Maximum power draw", "``calculate_flexibility_metrics(...)``", ":math:`P_{max} = \max_{t} P^{grid}_{t}`"
  "Discharge energy capacity", "``calculate_flexibility_metrics(...)``", ":math:`E_{dis} = \Delta t \sum_{t} \max\left(0, P_{base} - P^{grid}_{t}\right)`"
  "Discharge power capacity", "``calculate_flexibility_metrics(...)``", ":math:`P_{dis} = E_{dis}/T_{dis}`"
  "Levelized value/cost of flexibility", "``calculate_flexibility_metrics(...)``", ":math:`LVOF = \frac{\left(C^{base}_{elec} - \left(C_{energy} + C_{demand} + C_{customer} - R_{DR}\right)\right) - \left(C^{base}_{rep} - C_{rep}\right)}{E_{dis}}`"
  "No shutdowns during nonworking hours", "``add_working_hours_constraint(m)``", ":math:`\sum_{i=1}^{N_{RO}} SD_{t,i} = 0, \quad \forall t \in \mathcal{T}_{nonwork}`"
  "No startups during nonworking hours", "``add_working_hours_constraint(m)``", ":math:`\sum_{i=1}^{N_{RO}} SU_{t,i} = 0, \quad \forall t \in \mathcal{T}_{nonwork}`"
  "Non-flexible skid startup fixed", "``restrict_flexible_trains(m, num_flexible_trains)``", ":math:`SU_{t,i} = 0, \quad \forall i \in \mathcal{I}_{nonflex}, \forall t`"
  "Non-flexible skid shutdown fixed", "``restrict_flexible_trains(m, num_flexible_trains)``", ":math:`SD_{t,i} = 0, \quad \forall i \in \mathcal{I}_{nonflex}, \forall t`"


Flowsheet Specifications
------------------------
TODO: Currently table is incomplete because I want feedback on whether this organizational structure is a good idea.
Also, these values aren't called / defined except in the tutoral. Not the flowsheet.
The first flowsheet represents the ???? Santa Barbra plant, which flexibly varies recovery. The key paramters are given in the table below.

.. csv-table::
  :header: "Description", "Value", "Units"
  
  "Number of RO skids", "4", "dimensionless"
  "Minimum RO recovery", "0.40", "dimensionless"
  "Maximum RO recovery", "0.52", "dimensionless"
  "Nominal RO recovery", "0.465", "dimensionless"
  "Nominal RO flowrate", "337.67", "m3/h"

The second flowsheet represents the Water Replenishment District (WRD) ARC facility in Pico Rivera, CA. 
The plant is modeled by 4 RO trains, 3 UF pumps and 1 UV unit. The energy intensity of each RO train is a function of flowrate and recovery. 
This is the key difference between this implementation and flex_ro_1. The operational limits and costing values used in this flowsheet are based on the WRD plant.
The default parameters for the WRD RO and UF unit models represent those limits. 

.. csv-table::
  :header: "Description", "Value", "Units"
  
  "Number of RO skids", "4", "dimensionless"
  "Number of UF pumps", "3", "dimensionless"
  "Minimum RO recovery", "0.88", "dimensionless"
  "Maximum RO recovery", "0.925", "dimensionless"
  "Minimum RO flowrate", "0", "m3/h"
  "Maximum RO flowrate", "400", "m3/h"

References
----------
[1] Rao, A. et al. "Valuing energy flexibility from water systems", 2024. Nature Water, 2, 10, pp. 1028-1037, https://www.nature.com/articles/s44221-024-00316-4

Eventual TO-DO: Add references to the cross-cutting and WRD papers 