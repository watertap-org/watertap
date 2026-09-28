How to build a Flexible Reverse Osmosis Flowsheet
============================

Introduction
------------
This document describes the steps to build a multiperiod model to optimize the operation of a flexible reverse osmosis (RO) desalination plant. The model is built using the WaterTAP framework and is designed to take into account varying electricity prices, operational constraints, and energy intensity of the RO process.    

*This may not need to be it's own document but somewhere all steps should be listed-
It would make more sense for this to be a how-to guide I think.*

- Create Electricity Price signals file
    - Need both the energy rates with matching a matching time step to the model
- Defining the operational boundaries
This typically would include define the following:
    - Product water flowrate
    - Recovery
But these bounds could be defined by other factors such as:
    - Maximum flowrate limits through pumps or membranes
    - Minimum flowrate limits through pumps or membranes
    - Pressure limits 
    - Product water quality

- Developing a surrogate model for energy intensity. This  PySMO Options for this include:
    - Creating a flowsheet with WaterTAP unit models, performing a param sweep across the , like in the second flex_ro tutorial
    - Using operational data to directly develop a surrogate model. This option is limited if there are few data points throughout the operational range / if all data represents one nominal operating point.

- Build custom flowsheet
    - The flowsheet will be used by the Pricetaker model to determine energy consumption
- Build custom unit models (if needed)
    - Existing unit models can be used or adapted. Existing models are cataloged in the other flex ro document
    - The requirements are the recovery and energy intensity (seen in the unit model file)
- Adapt parameters (if needed)
    - The parameters class is used to define valid inputs and must be updated if any changes are made to defaults
    - This could include the flex desal model or any of the individual unit models
- Apply functions for constraints 
    - Existing functions are cataloged in the other flex ro document.
    - create new ones as needed. Diff. plants have unique constraints.
- Set cost values
    - Feed, brine, chemical costs
    - Replacement costs
    - Capital Costs
    - All described inthe other flex ro document
- Run Pricetaker optimization
    - The problem formulation is MINLP, meaning the default waterTAP solver cannot be used. Instead, one of the solvers in the tutorial should be selected.