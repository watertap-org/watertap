Electrocoagulation Costing Method
==================================


Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., ``m.fs.costing.electrocoagulation``) 
when applying the ``cost_electrocoagulation`` costing method in the ``watertap_costing_package``:


.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units", "Notes"

   "Reactor capital cost base parameter", ":math:`A_r`", ``reactor_capital_cost_base``, 11500, ":math:`\text{USD}_{2020}`", "Parameters from Table 2.1 for "Agitated Reactor" in Smith (2005)"
   "Reactor capital cost exponent", ":math:`b_r`", ``reactor_capital_cost_exponent``, 0.45, ":math:`\text{dimensionless}`", "Parameters from Table 2.1 for "Agitated Reactor" in Smith (2005)"
   "Reactor capital cost material coefficient", ":math:`x_{r,m}`", ``reactor_material_coeff``, 1.0, ":math:`\text{dimensionless}`", "1 for carbon steel; 3.4 for stainless steel; 0.062 for PVC"
   "Reactor capital cost safety factor", ":math:`SF_r`", ``reactor_capital_safety_factor``, 2.5, ":math:`\text{dimensionless}`", "Developed with feedback from industry experts"
   "Power supply capital cost equation slope", ":math:`A_p`", ``power_supply_capital_slope``, 0.51972, ":math:`\text{USD}_{2020}\text{ W}^{-1}`", "DC power supply + transformer + electrical connection base cost, developed from magna-power.com"
   "Flocculator capital cost base parameter", ":math:`A_f`", ``floc_capital_cost_base``, 1075700, ":math:`\text{USD}_{2007}`", "Figure 5.5.22 in McGivney & Kawamura (2008); refit to power equation"
   "Flocculator capital cost equation exponent", ":math:`b_f`", ``floc_capital_cost_exponent``, -0.95139, ":math:`\text{dimensionless}`", "Figure 5.5.22 in McGivney & Kawamura (2008); refit to power equation"
   "Sludge handling cost", ":math:`c_{sh}`", ``sludge_handling_cost``, 0.0, ":math:`\text{USD}\text{ kg}^{-1}`", "Cost of sludge handling is assumed to be zero by default"
   "Electrode material cost", ":math:`c_{mat}`", ``electrode_material_cost``, 2, ":math:`\text{USD}_{2021}\text{ kg}^{-1}`", "Cost per kg for electrode material; 2.23 for Al, 3.41 for Fe"
   "Electrode material cost safety factor", ":math:`SF_{mat}`", ``electrode_material_cost_safety_factor``, 2.0, ":math:`\text{dimensionless}`", "Developed with feedback from industry experts"


Costing Method Variables
++++++++++++++++++++++++

The following variables are constructed on the unit block (e.g., ``m.fs.unit.costing``) when 
applying the ``cost_electrocoagulation`` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"
   
   "Capital cost of reactor", ":math:`C_{r}`", ``capital_cost_reactor``, None, ":math:`\text{USD}`",
   "Capital cost of electrodes", ":math:`C_{e}`", ``capital_cost_electrodes``, None, ":math:`\text{USD}`",
   "Capital cost of power supply", ":math:`C_{p}`", ``capital_cost_power_supply``, None, ":math:`\text{USD}`",
   "Capital cost of floc reactor", ":math:`C_{f}`", ``capital_cost_floc_reactor``, None, ":math:`\text{USD}`",
   "Annual cost of sludge management", ":math:`C_{sh}`", ``annual_sludge_management``, None, ":math:`\text{USD year}^{-1}`",


Capital Cost Calculations
+++++++++++++++++++++++++

Capital costs for electrocoagulation are the summation of the capital cost of the reactor, electrodes, power supply, and flocculator.

.. math::
    C_{total} = C_r + C_e + C_p + C_f

The capital cost of the reactor is calculated according to:

.. math::
    C_r = \left( A_r  V_{r}^{b_r} \right) \left( x_{r,m} \right) \left( SF_r \right)

The capital cost of the electrodes is calculated from the mass of the electrodes:

.. math::
    C_e = \left( c_{mat} m_{electrode} \right) SF_{mat}

The capital cost of the power supply is determined from the power required for the electrocoagulation process:

.. math::
    C_p = A_p P_{tot}

The flocculator capital cost is a function of the flocculator volume:

.. math::
    C_f = A_f V_{floc}^{b_f}

Operating Cost Calculations
++++++++++++++++++++++++++++

Operating costs for electrocoagulation are the summation of the electrode replacement, 
electricity required, and the annual cost of sludge management.

Electricity costs are calculated with the power demand :math:`P_{tot}` on an annual basis:

.. math::
    C_{elec} = P_{tot} c_{elec}

Electrode replacement costs are a function of the dose of coagulant, volumetric flow (on an annual basis), 
and the cost of the electrode material:

.. math::
    C_e = \left( D_c q_{liq} c_{mat} m_{elec} \right) SF_{mat}

And the annual cost of sludge management is from the total annual mass flow of all non-water components from the `byproduct` port 
on the electrocoagulation unit model:

.. math::
    C_{sh} = \left(  \sum_{j} M_j \right) c_{sh}

Note: due to the uncertainty in the cost of sludge management, this cost is assumed to be zero by default (i.e., :math:`c_{sh} = 0`).
The user is encouraged to provide their own value for this parameter if they desire to include it in the costing calculations.

References
++++++++++

| W. McGivney and S. Kawamura (2008)
| Cost Estimating Manual for Water Treatment Facilities
| DOI: 10.1002/9780470260036

| Power supply cost estimation from magna-power.com
| Linear equation fit to SL and MT series cost data