How to build a Flexible Reverse Osmosis Flowsheet
=================================================

Introduction
------------

This document describes the steps to build a multiperiod model to optimize the
operation of a flexible reverse osmosis (RO) desalination plant. The model is
built using the IDAES Pricetaker framework and is designed to account for varying
electricity prices, operational constraints, and the energy intensity of the treatment
process.

The workflow typically includes the following steps:

1) Create electricity price signal files.

  - time series should match the start, end, and time step parameters of the model.

2) Define the operational boundaries for the RO plant. This typically includes the following:

  - Product water flowrate
  - Recovery

  Other bounds may also be relevant, such as:

  - Maximum flowrate limits through pumps or membranes
  - Minimum flowrate limits through pumps or membranes
  - Pressure limits
  - Product water quality

3) Develop a surrogate model for energy intensity.

  The approach used in the `Flexible RO Tutorial <watertap/tutorials/flex_ro/flex_recovery_and_flow/flex_recovery_and_flow.ipynb>`
  was to create a flowsheet with WaterTAP unit models and perform a parameter sweep over different flowrates and recoveries. Then, 
  the resulting data can be used to develop a surrogate model using PYSMO.

  A data-driven approach can be applied using operational data directly to form a surrogate.
  This option is limited when there are few data points throughout the operating range or when
  most data represent a single nominal operating point.
  

4) Build a flowsheet compatible with the PriceTaker model.

  - The flowsheet is used to determine energy consumption and to represent the operating limitations. See the example in `Flexible Recovery<watertap/flowsheets/flex_desal/flex_recovery_flowsheet.py>` and `Flexible Recovery and Flow<watertap/flowsheets/flex_desal/flex_recovery_and_flow_flowsheet.py>` for reference.
  

5) Build custom unit models, if needed.

  - Existing unit models for RO, UF, intake, posttreatment etc. can be used directly or adapted. 
  - The key requirements for a unit model are the recovery and energy intensity.

6) Adapt parameters, if needed.

  - The parameters class defines valid inputs and must be updated if defaults
    change.
  - This may include the flexible desalination model or any of the individual
    unit models.

7) Create Script to build the multiperiod model.

  - The script will build the multiperiod model. The steps can be seen in the `Flexible RO Tutorial <watertap/tutorials/flex_ro/flex_recovery_and_flow/flex_recovery_and_flow.ipynb>`.


8) Apply functions for constraints.

  - Existing functions are cataloged in `Flexible RO Documentation <watertap/docs/technical_reference/flowsheets/flex_ro.rst>`
  - Create new functions that represent the unique constraints of the system to be modeled.

9) Set other cost values, if required.

  - Feed, brine, and chemical costs
  - Replacement costs
  - Capital costs
  - These are passed as optional parameters to the unit models. 

10) Solve the PriceTaker problem formulation with optimization software.

  - The problem formulation is a MINLP, so the default WaterTAP solver is not
    appropriate. Instead, use Gurobi, SCIP, BARON, or other suitable MINLP solvers.
