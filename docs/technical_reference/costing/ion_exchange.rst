Ion Exchange Costing Method
============================

Costing Method Parameters
+++++++++++++++++++++++++

The following parameters are constructed for the unit on the FlowsheetCostingBlock (e.g., `m.fs.costing.ion_exchange`) when applying the `cost_ion_exchange` costing method in the ``watertap_costing_package``:


.. csv-table::
   :header: "Description", "Symbol", "Parameter Name", "Default Value", "Units", "Notes"

   "Anion exchange resin cost", ":math:`c_{res}`", "``anion_exchange_resin_cost``", "205", ":math:`\text{USD}_{2020}\text{/ft}^{3}`", "Assumes strong base polystyrenic gel-type Type II. From EPA-WBS cost model."
   "Cation exchange resin cost", ":math:`c_{res}`", "``cation_exchange_resin_cost``", "153", ":math:`\text{USD}_{2020}\text{/ft}^{3}`", "Assumes strong acid polystyrenic gel-type. From EPA-WBS cost model."
   "Regenerant dose per volume of resin", ":math:`D_{regen}`", "``regen_dose``", "300", ":math:`\text{kg/}\text{m}^{3}`", "Mass of regenerant chemical per cubic meter of resin volume"
   "Ion exchange column cost equation A coeff", ":math:`C_{col,A}`", "``vessel_A_coeff``", "1596.499", ":math:`\text{USD}_{2020}`", "Carbon steel w/ stainless steel internals. From EPA-WBS cost model."
   "Ion exchange column cost equation B coeff", ":math:`C_{col,b}`", "``vessel_b_coeff``", "0.459496", ":math:`\text{dimensionless}`", "Carbon steel w/ stainless steel internals. From EPA-WBS cost model."
   "Backwash/rinse tank cost equation A coeff", ":math:`C_{bw,A}`", "``backwash_tank_A_coeff``", "308.9371", ":math:`\text{USD}_{2020}`", "Steel tank. From EPA-WBS cost model."
   "Backwash/rinse tank cost equation B coeff", ":math:`C_{bw,b}`", "``backwash_tank_b_coeff``", "0.501467", ":math:`\text{dimensionless}`", "Steel tank. From EPA-WBS cost model."
   "Regeneration solution tank cost equation A coeff", ":math:`C_{regen,A}`", "``regen_tank_A_coeff``", "57.02158", ":math:`\text{USD}_{2020}`", "Stainless steel tank. From EPA-WBS cost model."
   "Regeneration solution tank cost equation B coeff", ":math:`C_{regen,b}`", "``regen_tank_b_coeff``", "0.729325", ":math:`\text{dimensionless}`", "Stainless steel tank. From EPA-WBS cost model."
   "Fraction of resin replaced per year", ":math:`f_{res}`", "``annual_resin_replacement_factor``", "0.05", ":math:`1/\text{yr}`", "Estimated 4-5% per year. From EPA-WBS cost model."
   "Minimum hazardous waste disposal cost", ":math:`f_{haz,min}`", "``hazardous_min_cost``", "3240", ":math:`\text{USD}_{2020}\text{/yr}`", "Minimum cost per hazardous waste shipment. From EPA-WBS cost model."
   "Unit cost for hazardous waste resin disposal", ":math:`f_{haz,res}`", "``hazardous_resin_disposal``", "347.10", ":math:`\text{USD}_{2020}\text{/ton}`", "From EPA-WBS cost model."
   "Unit cost for hazardous waste regeneration solution disposal", ":math:`f_{haz,regen}`", "``hazardous_regen_disposal``", "3.64", ":math:`\text{USD}_{2020}\text{/gal}`", "From EPA-WBS cost model."
   "Number of cycles the regenerant can be reused before disposal", ":math:`f_{recycle}`", "``regen_recycle``", "1", ":math:`\text{dimensionless}`", "Can optionally be set by the user to investigate more efficient regen regimes."
   "Costing factor to account for total installed cost installation of equipment", ":math:`f_{TIC}`", "``total_installed_cost_factor``", "1.65", ":math:`\text{dimensionless}`", "Costing factor to account for total installed cost of equipment"
   "Unit cost of NaCl", ":math:`c_{regen}`", "``costing.nacl``", "0.09", ":math:`\text{USD}_{2020}\text{/kg}`", "Assumes solid NaCl. From CatCost v 1.0.4"
   "Unit cost of HCl", ":math:`c_{regen}`", "``costing.hcl``", "0.17", ":math:`\text{USD}_{2020}\text{/kg}`", "Assumes 37% solution HCl. From CatCost v 1.0.4"
   "Unit cost of NaOH", ":math:`c_{regen}`", "``costing.naoh``", "0.59", ":math:`\text{USD}_{2020}\text{/kg}`", "Assumes 30% solution NaOH. From iDST"
   "Unit cost of Methanol (MeOH)", ":math:`c_{regen}`", "``costing.meoh``", "3.395", ":math:`\text{USD}_{2008}\text{/kg}`", "Assumes 100% pure MeOH. From ICIS"

Costing Method Variables
++++++++++++++++++++++++

The following variables are constructed on the unit block (e.g., m.fs.unit.costing) when applying the `cost_ion_exchange` costing method in the ``watertap_costing_package``:

.. csv-table::
   :header: "Description", "Symbol", "Variable Name", "Index", "Units"

   "Density of regenerant solution", ":math:`\rho_{regen}`", "``regen_soln_dens``", "None", ":math:`\text{kg/}\text{m}^{3}`"
   "Regenerant dose required for regeneration per volume of resin [kg regenerant/m3 resin]", ":math:`D_{regen}`", "``regen_dose``", "None", ":math:`\text{kg/}\text{m}^{3}`"
   "Capital cost for one vessel", ":math:`C_{col}`", "``capital_cost_vessel``", "None", ":math:`\text{USD}`"
   "Capital cost for resin for one vessel", ":math:`C_{resin}`", "``capital_cost_resin``", "None", ":math:`\text{USD}`"
   "Capital cost for regeneration solution tank", ":math:`C_{regen}`", "``capital_cost_regen_tank``", "None", ":math:`\text{USD}`"
   "Capital cost for backwash + rinse solution tank", ":math:`C_{bw}`", "``capital_cost_backwash_tank``", "None", ":math:`\text{USD}`"
   "Operating cost for hazardous waste disposal", ":math:`D_{regen}`", "``operating_cost_hazardous``", "None", ":math:`\text{USD/}\text{yr}`"
   "Regeneration solution flow", ":math:`\dot{v}_{regen}`", "``flow_mass_regen_soln``", "None", ":math:`\text{kg/}\text{yr}`"
   "Total pumping power required", ":math:`P_{tot}`", "``total_pumping_power``", "None", ":math:`\text{kW}`"

Capital Cost Calculations
+++++++++++++++++++++++++

Capital costs for ion exchange in the ``watertap_costing_package`` are the summation of the 
total cost of the resin, columns, backwashing tank, and regeneration solution tank:

Resin is costed based on the total volume of resin required for the system, where :math:`c_{res}` is the cost per volume of resin (either cation or anion exchange resin):

.. math::
    C_{resin} = V_{res,tot} c_{res}

Vessel cost as a function of volume was fit to a power function to determine capital cost of each column:

.. math::
    C_{col} = C_{col,A} V_{col}^{C_{col,b}}
   

The backwashing tank is assumed to include backwash and rinsing volumes. The total volume of this tank is:

.. math::
    V_{bw} = Q_{bw} t_{bw} + Q_{rinse} t_{rinse}

Backwashing tank cost as a function of volume was fit to a power function to determine capital cost of the backwashing tank:

.. math::
    C_{bw} = C_{bw,A} V_{bw}^{C_{bw,b}}
   
Regeneration tank cost as a function of volume was fit to a power function to determine capital cost of the regeneration tank:

.. math::
    C_{regen} = C_{regen,A} V_{regen}^{C_{regen,b}}

And the total capital cost for the ion exchange system is the summation of these:

.. math::
    C_{tot} = ((C_{resin} + C_{col}) (n_{op} + n_{red}) + C_{bw} + C_{regen}) f_{TIC}

A total installed cost (:math:`f_{TIC}`) factor of 1.65 is applied to account for installation costs. 

.. note::
    If using ``single_use`` option for ``regenerant`` configuration keyword, the capital for the regeneration tank is zero.

 
Operating Cost Calculations
+++++++++++++++++++++++++++


The operating costs for ion exchange includes the annual resin replacement cost, regeneration solution flow, energy consumption for booster pumps, 
and any hazardous waste handling costs.

Generally, the largest operating cost is the cost of the regeneration solution. The type of regeneration solution used is set via the 
optional model configuration keyword ``regenerant``. Costing data is available for the following regenerant chemicals:

* NaCl
* HCl
* NaOH
* MeOH

If the user does not provide a value for this option, the model defaults to a NaCl regeneration solution. The dose of regenerant needed
is set by the parameter ``regen_dose`` in kg regenerant per cubic meter of resin volume. The mass flow of regenerant solution [kg/yr] is:

.. math::
    \dot{m}_{regen} = \frac{D_{regen} V_{res} (n_{op} + n_{red})}{t_{cycle} f_{recycle}}

Annual resin replacement cost is:

.. math::
    C_{op,res} = V_{res} (n_{op} + n_{red}) f_{res} c_{res}

If the spent resin and regenerant contains hazardous material, the user designates this by the model configuration keyword ``hazardous_waste``. If set to ``True``, hazardous
disposal costs are calculated as a function of the annual mass of resin replaced and regenerant consumed:

.. math::
    C_{op,haz} = f_{haz,min} + \bigg( M_{res} (n_{op} + n_{red}) f_{res} \bigg)  f_{haz,res} + \dot{v}_{regen} f_{haz,regen}

Where :math:`M_{res}` is the resin mass for a single bed and :math:`\dot{v}_{regen}` is the volumetric flow of regenerant solution. If ``hazardous_waste`` is set to ``False``,
:math:`C_{op,haz} = 0`

The total energy consumed by the unit is the summation of the power required for each of the booster pump, backwashing pump, regeneration pump, and rinsing pump. Each is scaled 
by the total time required for each step:

.. math::
    P_{tot} = \cfrac{P_{main} t_{break} + P_{bw} t_{bw} + P_{regen} t_{regen} + P_{rinse} t_{rinse}}{t_{cycle}} 

If the user chooses ``single_use`` for the ``regenerant`` configuration keyword, there is no cost for regeneration solution:

.. math::
    \dot{m}_{regen} = \dot{v}_{regen} = 0

Instead, the model assumes the entire volume of resin for the operational columns is replaced at the end of each service cycle by calculating the 
volumetric "flow" of resin:

.. math::
    \dot{v}_{resin} = \frac{V_{res, tot}}{t_{break}} 

And then operational cost of replacing the entire bed is:

.. math::
    C_{op,res} = \dot{v}_{resin} c_{res}

If ``hazardous_waste`` is set to ``True``, the hazardous waste disposal costs are: 

.. math::
    C_{op,haz} = f_{haz,min} + ( \dot{v}_{resin} \rho_{b} n_{op})  f_{haz,res}

Otherwise, :math:`C_{op,haz} = 0` as before. 

Lastly, the total energy consumed by the unit for ``single_use`` configuration includes the booster pump, backwashing pump, and rinsing pump:

.. math::
    P_{tot} = \cfrac{P_{main} t_{break} + P_{bw} t_{bw} + P_{rinse} t_{rinse}}{t_{cycle}} 

Code Documentation
------------------

* :mod:`watertap.costing.unit_models.ion_exchange`

References
----------
| United States Environmental Protection Agency. (2021). Work Breakdown Structure-Based Cost Models
| https://www.epa.gov/sdwa/drinking-water-treatment-technology-unit-cost-models

| CatCost https://catcost.chemcatbio.org/
| v 1.0.4 available here: https://datahub.chemcatbio.org/dataset/catcost-v1-0-4

| Integrated Decision Support Tool (i-DST) 
| https://idst.mines.edu/