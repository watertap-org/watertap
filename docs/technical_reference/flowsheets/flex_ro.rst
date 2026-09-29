Flexible Reverse Osmosis Flowsheets
===================================

Introduction
------------

The flowsheets ``flex_recovery.py`` and ``flex_recovery_and_flow.py`` represent two different RO plants.
Both flowsheets use the IDAES `Pricetaker model <https://github.com/IDAES/idaes-pse/blob/main/docs/reference_guides/apps/grid_integration/multiperiod/Price_Taker.rst>`_
to determine the cost-optimal operation based on treatment energy requirement, operational constraints, and variable grid electricity costs.
Each flowsheet characterizes the operational flexiblity of the respective plants.


File Structure
---------------

There are four files associated with each case study.
1) The flowsheet file (``flex_recovery.py`` or ``flex_recovery_and_flow.py``) contains the model build, constraints, and expressions used to construct the Pricetaker optimization problem.
2) The unit model file (``unit_models.py``) contains the unit model definitions for each treament step in the flowsheet.
3) The params file (``params.py``) contains the parameter definitions for each unit model and the overall FlexDesal model.
4) The utilities file (``utils.py``) includes a few additional functions.
 
Edits to the model typically invlove changes to the flowsheet, unit model, and parameters files.

Pricetaker Model Functions
--------------------------
To setup of the Pricetaker model, several helper function are used. These functions are part of the Pricetaker framework and are defined in IDAES.

* ``append_lmp_data``
  Appends the locational marginal price (LMP) time series to the PriceTaker model. This is used to calculate 
  the consumption or energy charge portion of the electricity costs.

  .. csv-table::
     :header: "Description", "Symbol", "Name", "Index", "Units"

     "Locational marginal price", ":math:`\lambda`", "lmp", "[t]", ":math:`\$/\text{kWh}`"

* ``build_multiperiod_model``
  Builds the full multiperiod PriceTaker model, attaching the plant flowsheet,
  operational constraints, and objective-function terms for every time period.
  This step links the WaterTAP desalination model to the IDAES price-taker
  framework and initializes the time-indexed operating model.

* ``update_operation_params``
  Updates the already-constructed PriceTaker model with scenario-specific operating
  data such as tariff rates, emissions intensity, and onsite generation factors.
  This allows the same model structure to be reused for different pricing or
  operational cases without rebuilding the full model.


Parameters
----------
Parameters from both flowsheets are described here. Note that not all are used in each flowsheet.

The top-level ``FlexDesalParams`` values include:

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default value
     - Units
     - Description
   * - ``start_date``
     - :math:`2022\text{-}07\text{-}05\ 00\text{:}00\text{:}00`
     - :math:`\text{timestamp}`
     - Start of the simulation horizon.
   * - ``end_date``
     - :math:`2022\text{-}07\text{-}06\ 00\text{:}00\text{:}00`
     - :math:`\text{timestamp}`
     - End of the simulation horizon.
   * - ``timestep_hours``
     - :math:`0.25`
     - :math:`\text{hour}`
     - Length of each simulation time step.
   * - ``product_water_price``
     - :math:`0`
     - :math:`\$/\text{m}^{3}`
     - Unit revenue for produced water.
   * - ``fixed_monthly_cost``
     - :math:`766000`
     - :math:`\$/\text{month}`
     - Fixed monthly customer/facility charge.
   * - ``customer_rate``
     - :math:`100`
     - :math:`\text{dimensionless}`
     - Customer-rate multiplier used in tariff calculations.
   * - ``constrain_to_baseline_production``
     - :math:`\text{False}`
     - :math:`\text{None}`
     - Whether to enforce baseline production tracking.
   * - ``curtailment_fraction``
     - :math:`0.0`
     - :math:`\text{dimensionless}`
     - Allowed fractional curtailment relative to baseline production.
   * - ``annual_production_AF``
     - :math:`3125`
     - :math:`\text{acre-ft/year}`
     - Annual production target used for absolute production constraints.
   * - ``production_constraint_to_objective``
     - :math:`\text{False}`
     - :math:`\text{None}`
     - Whether production compliance is enforced through objective penalization.
   * - ``production_constraint_penalty``
     - :math:`0.6`
     - :math:`\$/\text{m}^{3}`
     - Penalty weight applied when production target is incorporated in the objective.
   * - ``emissions_cost``
     - :math:`0`
     - :math:`\$/\text{kg}`
     - Cost assigned to emissions associated with grid electricity use.
   * - ``include_demand_response``
     - :math:`\text{False}`
     - :math:`\text{None}`
     - Enable demand-response price/revenue terms.
   * - ``include_battery``
     - :math:`\text{False}`
     - :math:`\text{None}`
     - Enable battery operation model.
   * - ``include_onsite_solar``
     - :math:`\text{False}`
     - :math:`\text{None}`
     - Enable onsite solar generation model.
   * - ``onsite_capacity``
     - :math:`0`
     - :math:`\text{kW}`
     - Installed onsite generation capacity.
   * - ``nonworking_hours``
     - :math:`[]`
     - :math:`\text{hour index list}`
     - Hours where startup/shutdown actions may be restricted.
   * - ``CAPEX_yr``
     - :math:`\text{None}`
     - :math:`\$/\text{year}`
     - Optional annualized CAPEX value for economic reporting.
   * - ``max_daily_shutdowns``
     - :math:`\text{None}`
     - :math:`\text{count/day}`
     - Optional limit on shutdown events over a daily rolling window.

Additional FlexDesal Parameter Details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The following optional parameters change the formulation during the model build (``build_multiperiod_model(m)``).

* ``include_demand_response``
  Enables the demand-response revenue term in the objective. When active, the
  plant can earn revenue by reducing grid power below a baseline reference level
  during event periods, and the model tracks the corresponding time-varying
  demand-response price.

  .. csv-table::
     :header: "Description", "Symbol", "Name", "Index", "Units"

     "Demand-response price", ":math:`\pi^{DR}_{t}`", "demand_response_price", "[t]", ":math:`\$/\text{kWh}`"
     "Baseline power reference", ":math:`P^{base}`", "baseline_power", "None", ":math:`\text{kW}`"
     "Demand-response revenue", ":math:`R_{DR}`", "demand_response_revenue", "[t]", ":math:`\$`"

* ``include_onsite_solar``
  Adds an onsite photovoltaic generation model and computes the power supplied to
  the plant from solar generation. Any unused generation is treated as excess
  solar power that can be curtailed rather than exported or stored.

  .. csv-table::
     :header: "Description", "Symbol", "Name", "Index", "Units"

     "Solar capacity factor", ":math:`CF_t`", "capacity_factor", "[t]", "dimensionless"
     "Onsite solar power used", ":math:`P^{solar}_{t}`", "power_utilized", "[t]", ":math:`\text{kW}`"
     "Excess solar power", ":math:`P^{excess}_{t}`", "excess_solar_power", "[t]", ":math:`\text{kW}`"

* ``include_battery(m)``
  Adds a battery storage model that can charge or discharge to shift electricity
  purchases across time. This supports temporal arbitrage and defers demand during
  higher-cost periods while conserving the plant's operating flexibility.

  .. csv-table::
     :header: "Description", "Symbol", "Name", "Index", "Units"

     "Battery charge rate", ":math:`P^{chg}_{t}`", "charge_rate", "[t]", ":math:`\text{kW}`"
     "Battery discharge rate", ":math:`P^{dis}_{t}`", "discharge_rate", "[t]", ":math:`\text{kW}`"
     "Battery energy state", ":math:`E_t`", "holdup", "[t]", ":math:`\text{kWh}`"


Unit Models
-----------
The existing flowsheets are built with the following unit models:

Intake
~~~~~~

The intake unit uses the ``IntakeParams`` dataclass:


.. list-table::
   :header-rows: 1

   * - Parameter
     - Default value
     - Units
     - Description
   * - ``energy_intensity``
     - :math:`0.157121734`
     - :math:`\text{kWh}/\text{m}^{3}`
     - Specific intake energy intensity.
   * - ``minimum_flowrate``
     - :math:`1063.5`
     - :math:`\text{m}^{3}/\text{h}`
     - Minimum intake flowrate.
   * - ``nominal_flowrate``
     - :math:`1063.5`
     - :math:`\text{m}^{3}/\text{h}`
     - Nominal intake flowrate.
   * - ``maximum_flowrate``
     - :math:`1063.5`
     - :math:`\text{m}^{3}/\text{h}`
     - Maximum intake flowrate.
   * - ``feed_cost``
     - :math:`\text{None}`
     - :math:`\$/\text{m}^{3}`
     - Optional feed-water cost.
   * - ``chemical_cost``
     - :math:`\text{None}`
     - :math:`\$/\text{m}^{3}`
     - Optional chemical cost.

Pretreament (Generic)
~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default value
     - Units
     - Description
   * - ``allow_shutdown``
     - :math:`\text{False}`
     - :math:`\text{None}`
     - Enables pretreatment unit on/off commitment logic.
   * - ``energy_intensity``
     - :math:`0.01`
     - :math:`\text{kWh}/\text{m}^{3}`
     - Specific pretreatment energy intensity.
   * - ``leakage_fraction``
     - :math:`0`
     - :math:`\text{dimensionless}`
     - Fraction of inlet flow not recovered in pretreatment.
   * - ``minimum_downtime``
     - :math:`0`
     - :math:`\text{num time steps}`
     - Minimum number of time steps the pretreatment unit remains off after shutdown.
   * - ``startup_delay``
     - :math:`0`
     - :math:`\text{num time steps}`
     - Delay between startup command and pretreatment operation.
   * - ``chemical_cost``
     - :math:`\text{None}`
     - :math:`\$/\text{m}^{3}`
     - Optional pretreatment chemical cost.



Ultrafiltration (UF)
~~~~~~~~~~~~~~~~~~~~

The UF unit uses the ``UFParams`` dataclass. This is an alternative pretreatment to the generic unit model above.

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default value
     - Units
     - Description
   * - ``num_uf_pumps``
     - :math:`4`
     - :math:`\text{None}`
     - Number of UF pumps represented in the model.
   * - ``minimum_operating_pumps``
     - :math:`1`
     - :math:`\text{None}`
     - Minimum number of UF pumps that must operate.
   * - ``allow_shutdown``
     - :math:`\text{True}`
     - :math:`\text{None}`
     - Enables UF unit on/off commitment logic.
   * - ``minimum_flowrate``
     - :math:`344`
     - :math:`\text{m}^{3}/\text{h}`
     - Minimum UF pump flowrate when operating (m3/h).
   * - ``nominal_flowrate``
     - :math:`900`
     - :math:`\text{m}^{3}/\text{h}`
     - Nominal UF pump flowrate (m3/h).
   * - ``maximum_flowrate``
     - :math:`989`
     - :math:`\text{m}^{3}/\text{h}`
     - Maximum UF pump flowrate (m3/h).
   * - ``nominal_recovery``
     - :math:`1`
     - :math:`\text{dimensionless}`
     - Nominal UF recovery.
   * - ``minimum_uptime``
     - :math:`2`
     - :math:`\text{time steps}`
     - Minimum number of time steps the UF pump stays on after startup.
   * - ``minimum_downtime``
     - :math:`2`
     - :math:`\text{time steps}`
     - Minimum number of time steps the UF pump stays off after shutdown.
   * - ``startup_delay``
     - :math:`1`
     - :math:`\text{time steps}`
     - Delay (time steps) between startup command and operation.
   * - ``allow_variable_recovery``
     - :math:`\text{False}`
     - :math:`\text{None}`
     - Recovery is fixed at nominal value for this WRD UF model.
   * - ``surrogate_type``
     - :math:`\text{quadratic_energy_intensity}`
     - :math:`\text{None}`
     - UF energy intensity surrogate form.
   * - ``surrogate_a``, ``surrogate_b``, ``surrogate_c``
     - :math:`1`, :math:`1`, :math:`1`
     - :math:`\text{surrogate-dependent}`
     - Coefficients for the quadratic UF energy intensity surrogate.


Reverse Omosis (RO)
~~~~~~~~~~~~~~~~~~~

The RO unit uses the ``ROParams`` dataclass.

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default value
     - Units
     - Description
   * - ``num_ro_skids``
     - :math:`4`
     - :math:`\text{num of skids}`
     - Number of RO skids represented in the model.
   * - ``minimum_operating_skids``
     - :math:`2`
     - :math:`\text{None}`
     - Minimum number of RO skids that must operate.
   * - ``allow_shutdown``
     - :math:`\text{True}`
     - :math:`\text{None}`
     - Enables RO skid on/off commitment logic.
   * - ``minimum_flowrate``
     - :math:`0`
     - :math:`\text{m}^{3}/\text{h}`
     - Minimum RO skid flowrate when operating (m3/h).
   * - ``nominal_flowrate``
     - :math:`337.670`
     - :math:`\text{m}^{3}/\text{h}`
     - Nominal RO skid flowrate (m3/h).
   * - ``maximum_flowrate``
     - :math:`400`
     - :math:`\text{m}^{3}/\text{h}`
     - Maximum RO skid flowrate (m3/h).
   * - ``minimum_recovery``
     - :math:`0.88`
     - :math:`\text{dimensionless}`
     - Minimum RO recovery.
   * - ``nominal_recovery``
     - :math:`0.92`
     - :math:`\text{dimensionless}`
     - Nominal RO recovery.
   * - ``maximum_recovery``
     - :math:`0.925`
     - :math:`\text{dimensionless}`
     - Maximum RO recovery.
   * - ``minimum_uptime``
     - :math:`2`
     - :math:`\text{time steps}`
     - Minimum number of time steps a skid stays on after startup.
   * - ``minimum_downtime``
     - :math:`2`
     - :math:`\text{time steps}`
     - Minimum number of time steps a skid stays off after shutdown.
   * - ``startup_delay``
     - :math:`1`
     - :math:`\text{time steps}`
     - Delay (time steps) between startup command and operation.
   * - ``max_num_skids_shutdown_per_timestep``
     - :math:`1`
     - :math:`\text{num of skids}`
     - Number of skids that can shutdown in a single time step
   * - ``allow_variable_recovery``
     - :math:`\text{False}`
     - :math:`\text{None}`
     - Recovery is bounded but not optimized as a free variable in default setup.
   * - ``surrogate_type``
     - :math:`\text{constant_energy_intensity}`
     - :math:`\text{None}`
     - RO surrogate type selector for the WRD case.
   * - ``surrogate_file``
     - :math:`\text{None}`
     - :math:`\text{path}`
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
     - :math:`0.41`
     - :math:`\text{kWh}/\text{m}^{3}`
     - Specific UV energy intensity.
   * - ``leakage_fraction``
     - :math:`0`
     - :math:`\text{dimensionless}`
     - Fraction of inlet flow not recovered in UV posttreatment.
   * - ``chemical_cost``
     - :math:`\text{None}`
     - :math:`\$/\text{m}^{3}`
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
     - :math:`0.1`
     - :math:`\text{kWh}/\text{m}^{3}`
     - Specific brine-discharge energy intensity.
   * - ``brine_cost``
     - :math:`\text{None}`
     - :math:`\$/\text{m}^{3}`
     - Optional brine-disposal cost.


Common Variables
-----------------

There are several variables common across the unit models.

.. csv-table::
  :header: "Variable", "Symbol", "Set"

  "Flowrate", ":math:`Q` (e.g., :math:`Q^{RO}_{t,i}`, :math:`Q^{UF}_{t,i}`)", "[t,i]"
  "Recovery", ":math:`R`", "[t,i]"
  "Operating mode", ":math:`u` (e.g., :math:`u^{brine}_{t}`)", "[t,i]"
  "Startup", ":math:`SU`", "[t,i]"
  "Shutdown", ":math:`s`", "[t,i]"
  "Energy intensity", ":math:`EI`", "[t,i]"
  "Flow cost", ":math:`c` (e.g., :math:`c^{intake}_{t}`, :math:`c^{brine}_{t}`)", "None"

Unit Model Common Equations
---------------------------

The following equations are built for every unit model via ``_add_required_variables``
in "unit_models".

.. csv-table::
   :header: "Description", "Equation"

   "Mass balance", ":math:`Q^{feed} = Q^{product} + Q^{reject}`"
   "Product flowrate", ":math:`Q^{product} = Q^{feed} \cdot R`"
   "Power consumption", ":math:`P = EI \cdot Q^{product}`"


Required Helper Functions
--------------------------
To complete setup, some WaterTAP defined function are used.

* ``add_demand_and_fixed_costs``
  Adds the horizon-level demand-charge variables and binds them to the per-period
  grid-power usage through lower-bound constraints. This helper is used to capture
  both fixed demand charges and variable demand charges, along with the fixed
  customer fee for the modeled billing horizon.

  .. csv-table::
     :header: "Description", "Symbol", "Name", "Index", "Units"

     "Fixed demand charge", ":math:`C_{demand}^{fixed}`", "fixed_demand_cost", "None", ":math:`\$`"
     "Variable demand charge", ":math:`C_{demand}^{var}`", "variable_demand_cost", "None", ":math:`\$`"
     "Fixed monthly customer cost", ":math:`C_{customer}`", "fixed_monthly_cost", "None", ":math:`\$`"

* ``constrain_water_production``
  Enforces the plant production requirement over the simulation horizon. Depending
  on the selected configuration, it either tracks a baseline production target
  (with optional curtailment) or enforces a fixed absolute target based on the
  annual production requirement.

  .. csv-table::
     :header: "Description", "Symbol", "Name", "Index", "Units"

     "Curtailment fraction", ":math:`\phi_{curtail}`", "curtailment_fraction", "None", ":math:`\text{dimensionless}`"
     "Baseline production", ":math:`W_{base}`", "baseline_production", "None", ":math:`\text{m}^3`"
     "Absolute production target", ":math:`W_{target}`", "production_target_abs", "None", ":math:`\text{m}^3`"


Optional Helper Functions
-------------------------
A number of helper functions are defined in the ``watertap.flowsheets.flex_desal.flex_recovery_and_flow.flex_recovery_and_flow_flowsheet``.
These helper function is build additional variables, constraints, and equations when called.

* ``add_flow_costs(m)``
  Adds total feed, brine-discharge, and chemical-cost expressions.

  .. csv-table::
     :header: "Description", "Symbol", "Name", "Index", "Units"

     "Total feed cost", ":math:`C_{feed}`", "total_feed_cost", "None", ":math:`\$`"
     "Total brine-discharge cost", ":math:`C_{brine}`", "total_brine_cost", "None", ":math:`\$`"
     "Total chemical cost", ":math:`C_{chem}`", "total_chemical_cost", "None", ":math:`\$`"

* ``add_flow_changes_penalty(m)``
  Adds binary indicators for RO and UF flow changes between consecutive time
  steps and penalizes such changes in the objective. This discourages frequent
  ramping and stabilizes operational patterns.

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
  Post-processes solved results to estimate flexibility metrics such as charge/discharge
  capacities and levelized value or cost of flexibility. These metrics are defined in Rao et al. [1].

  .. csv-table::
     :header: "Description", "Symbol", "Name", "Index", "Units"

     "Maximum power draw", ":math:`P_{max}`", "maximum_power", "None", ":math:`\text{kW}`"
     "Discharge energy capacity", ":math:`E_{dis}`", "energy_capacity", "None", ":math:`\text{kWh}`"
     "Discharge power capacity", ":math:`P_{dis}`", "power_capacity", "None", ":math:`\text{kW}`"
     "Levelized value/cost of flexibility", ":math:`LVOF`", "LVOF", "None", ":math:`\$/\text{kWh}`"

* ``begin_and_end_constraint(m)``
  Enforces a cyclic operating-state condition by matching RO train 1 operating mode at the first and
  last time points.

* ``add_working_hours_constraint(m)``
  Prevents RO startup and shutdown actions during configured nonworking hours. "Working hours" are based on operator availability
  to enact chnages to operating state of the system.

* ``restrict_flexible_trains(m, num_flexible_trains)``
  Fixes startup and shutdown decisions to zero for the non-flexible RO trains so
  that only a selected subset of skids can participate in the flexible operation
  strategy.


Equations and Relationships
---------------------------

Each equation below is cross-referenced to the helper-function section where its
left-hand-side quantity is introduced.

.. csv-table::
  :header: "Description", "Defined in", "Equation"

  "Grid power balance", "``build_desal_flowsheet(...)``", ":math:`P^{grid}_{t} = P^{process}_{t} - P^{solar}_{t} + P^{chg}_{t} - P^{dis}_{t} + P^{excess}_{t}`"
  "Demand-response revenue", "``add_operational_cost_expressions(...)``", ":math:`R^{DR}_{t} = \pi^{DR}_{t} \left(P^{base} - P^{grid}_{t}\right) \Delta t`"
  "Energy charge", "``add_operational_cost_expressions(...)``", ":math:`C_{energy} = \Delta t \sum_{t} \lambda_{t} P^{grid}_{t}`"
  "Emissions cost", "``add_operational_cost_expressions(...)``", ":math:`C_{emissions} = \Delta t \sum_{t} \alpha_{t} P^{grid}_{t}`"
  "Fixed demand charge", "``add_demand_and_fixed_costs(m)``", ":math:`C_{demand,t}^{fixed} \geq r_{demand,t}^{fixed} P^{grid}_{t} N_{months}`"
  "Variable demand charge", "``add_demand_and_fixed_costs(m)``", ":math:`C_{demand,t}^{var} \geq r_{demand,t}^{var} P^{grid}_{t} N_{months}`"
  "Fixed monthly customer cost", "``add_demand_and_fixed_costs(m)``", ":math:`C_{customer} = c_{customer}^{monthly} N_{months}`"
  "Baseline production target", "``constrain_water_production(m, baseline_production)``", ":math:`W_{prod} \geq W_{base} \left(1 - \phi_{curtail}\right)`"
  "Absolute production target", "``constrain_water_production(m)``", ":math:`W_{prod} \geq W_{target}`"
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


Flowsheet Specifications in Tutorials
-------------------------------------
NOTE: These values aren't called / defined except in the tutoral. Not the flowsheet.
The first flowsheet represents the Santa Barbra plant, which flexibly varies recovery.

.. csv-table::
  :header: "Description", "Value", "Units"
  
  "Number of RO skids", ":math:`4`", ":math:`\text{dimensionless}`"
  "Minimum RO recovery", ":math:`0.40`", ":math:`\text{dimensionless}`"
  "Maximum RO recovery", ":math:`0.52`", ":math:`\text{dimensionless}`"
  "Nominal RO recovery", ":math:`0.465`", ":math:`\text{dimensionless}`"
  "Nominal RO flowrate", ":math:`337.67`", ":math:`\text{m}^{3}/\text{h}`"

The second flowsheet represents the Water Replenishment District (WRD) ARC facility in Pico Rivera, CA. 
The plant is modeled by 4 RO trains, 3 UF pumps and 1 UV unit. The energy intensity of each RO train is a function of flowrate and recovery. 
This is the key difference between this implementation and flex_ro_1. The operational limits and costing values used in this flowsheet are based on the WRD plant.
The default parameters for the WRD RO and UF unit models represent those limits. 

.. csv-table::
  :header: "Description", "Value", "Units"
  
  "Number of RO skids", ":math:`4`", ":math:`\text{dimensionless}`"
  "Number of UF pumps", ":math:`3`", ":math:`\text{dimensionless}`"
  "Minimum RO recovery", ":math:`0.88`", ":math:`\text{dimensionless}`"
  "Maximum RO recovery", ":math:`0.925`", ":math:`\text{dimensionless}`"
  "Minimum RO flowrate", ":math:`0`", ":math:`\text{m}^{3}/\text{h}`"
  "Maximum RO flowrate", ":math:`400`", ":math:`\text{m}^{3}/\text{h}`"

References
----------
[1] Rao, A. et al. "Valuing energy flexibility from water systems", 2024. Nature Water, 2, 10, pp. 1028-1037, https://www.nature.com/articles/s44221-024-00316-4

Eventual TO-DO: Add references to the cross-cutting and WRD papers 
