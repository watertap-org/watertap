Stoichiometric Reactor Costing Method
======================================

Currently, the costing method is implemented for lime and soda ash softening and acidification which only include
the capital cost of building the reactor. The capital cost of lime soda ash is a function of 
total reagent mass being added to the softening process and is only valid when both precipitant and reagents are provided.
While acid addition capital cost is only constructed if only reagents are provided. Acid addition costing is 
base on folume flow of acid per day. 
(Please refer to the `stoichiometric reactor documentation <https://watertap.readthedocs.io/en/stable/technical_reference/unit_models/stoichiometric_reactor.html>`_ for details on dissolution and precipitation reaction configurations). 

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., `m.fs.costing.stoichiometric_reactor`) when applying the `cost_stoichiometric_reactor` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units"

   "Softening capital cost", ":math:`C_{softening}`", "``capital_cost_softening``", "374.9", ":math:`\text{USD/(lb/day)}`"
   "Acid addition capital cost", ":math:`C_{acid}`", "``capital_cost_acid_addition``", "127.8", ":math:`\text{USD/(gallon/day)}`"

Costing Method Variables
++++++++++++++++++++++++

There are not unique costing variables constructed for stoichiometric reactor

Capital Cost Calculations
+++++++++++++++++++++++++

(if user only includes reagents).(Refer to `stoichiometric reactor documentation <https://watertap.readthedocs.io/en/stable/technical_reference/unit_models/stoichiometric_reactor.html>`_) 

    .. math::

        CC_{softening}=C_{softening}*\sum{M_{reagent}}

        CC_{acidification}=C_{acid}*\sum{Q_{reagent}}


 where M_{reagent} is mass flow all reagent added in softening processes, and Q_{reagent} is flow volume of chemical added. 

Operating Cost Calculations
+++++++++++++++++++++++++++

There are no operational costs provided, user needs to manually cost the reagent addition per mass/flow of reagent as in example below for 
softening operation with soda ash and lime addition. 
.. code-block::

   # build the unit model 
    reagents = {
         "Na2CO3": {
               "mw": 105.99 * pyunits.g / pyunits.mol,
               "dissolution_stoichiometric": {"Na_+": 2, "HCO3_-": 1},
         },
         "CaO": {
               "mw": 56.0774 * pyunits.g / pyunits.mol,
               "dissolution_stoichiometric": {"Ca_2+": 1, "H2O": 1},
         },
      }
    m.fs.chemical_addition = StoichiometricReactor(
         property_package=m.fs.properties,
         reagent=reagents,
      )
    # The user must the specify how much reagent to add
    m.fs.chemical_addition.reagent_dose["Na2CO3"].fix(1e-3)
    m.fs.chemical_addition.reagent_dose["CaO"].fix(1e-3)

    # specify the costs for lime (CaO)
    blk.lime_cost = Param(
        initialize=0.13,
        units=m.fs.costing.base_currency / pyunits.kg,
        mutable=True,
    )
    # specify the costs for soda ash (Na2CO3)
    blk.soda_ash_cost = Param(
        initialize=0.13,
        units=m.fs.costing.base_currency / pyunits.kg,
        mutable=True,
    )
    # Register the cost for each chemical being added
    m.fs.costing.register_flow_type("lime_cost", blk.lime_cost )
    m.fs.costing.register_flow_type("soda_ash_cost", blk.soda_ash_cost )
    
    # Register the mass or volume flow for each chemical being added
    m.fs.costing.cost_flow(
        blk.chemical_precipitator.flow_mass_reagent["CaO"],
        "lime_cost",
    )
    m.fs.costing.cost_flow(
        blk.chemical_precipitator.flow_mass_reagent["Na2CO3"],
        "soda_ash_cost",
    )
 
Code Documentation
------------------

* :mod:`watertap.costing.unit_models.stoichiometric_reactor`

References
----------
O. Amusat, A. Atia, A. Dudchenko, T. Bartholomew, “Modeling Framework for Cost Optimization of Process-Scale Desalination Systems with Mineral Scaling and Precipitation”, ACS ES&T Engr, (2024)