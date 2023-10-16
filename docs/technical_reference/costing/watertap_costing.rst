.. _watertap_costing:

WaterTAP Costing Package
========================

.. currentmodule:: watertap.costing.watertap_costing

The WaterTAP costing module contains the costing package typically used for high-fidelity models, but it can also be utilized with zero-order models. Technoeconomic data is set on the package or utilizing global parameters added by the unit model costing method.

Usage
-----

The WatertTAPCosting class contains all the variables and constraints needed to cost a unit model. It also inherits the functionality of the :ref:`WaterTAPCostingBlockData`.

The code below shows an outline of how the WatertTAPCosing class is intended to be used to cost unit model.

.. testcode::

    import pyomo.environ as pyo
    import idaes.core as idc
    from watertap.costing import (
        WaterTAPCosting,
        register_costing_parameter_block,
        make_capital_cost_var,
    )

    def build_my_unit_model_param_block(blk):
        """
        This function builds the global parameters for MyUnitModel.

        This function should also register needed flows using the
        blk.parent_block().register_flow_type method on the costing package.
        """
        blk.fixed_capital_cost = pyo.Var(
            initialize=42,
            doc="Fixed capital cost for all of my units",
            units=pyo.units.USD_2020,
        )

    # This decorator ensures that the function
    # `build_my_unit_model_param_block` is only
    # added to the costing package once.
    # It registers it as a sub-block with the
    # name `my_unit`.
    @register_costing_parameter_block(
        build_rule=build_my_unit_model_param_block,
        parameter_block_name="my_unit",
    )
    def cost_my_unit_model(blk):
        """
        Cost an instance of MyUnitModel
        """
        # creates the `capital_cost` Var
        make_capital_cost_var(blk)

        # here we reference the `fixed_capital_cost` parameter
        # automatically added by the `register_costing_parameter_block`
        # decorator.
        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.costing_package.my_unit.fixed_capital_cost,
            name="fixed capital cost constraint",
        )

    @idc.declare_process_block_class("MyUnitModel")
    class MyUnitModelData(idc.UnitModelBlockData):

        @property
        def default_costing_method(self):
            # could point to a static method on
            # this class, could be function in 
            # a different module even
            return cost_my_unit_model

    m = pyo.ConcreteModel()
    m.fs = idc.FlowsheetBlock(dynamic=False)
    m.fs.costing = WaterTAPCosting()

    m.fs.my_unit_1 = MyUnitModel()
    
    # The `default_costing_method_attribute` on the
    # unit model is checked, and the function
    # `cost_my_unit_model` returned then build the costing block.
    # This method also adds the `my_unit` global parameter block,
    # so the global costing parameter m.fs.costing.my_unit.fixed_capital_cost
    # is the same for all instances of MyUnitModel.
    m.fs.my_unit_1.costing = idc.UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
    )

    m.fs.my_unit_2 = MyUnitModel()

    # Here everythin as before, but the global parameter block
    # m.fs.costing.my_unit is not re-built. 
    m.fs.my_unit_2.costing = idc.UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
    )


Costing WaterTAP Models
-----------------------

The WaterTAPCostingData class includes variables and constraints necessary to calculate process-wide costs:

=============================================  ====================  ===============================================  ==============================================================================
                 Cost                               Variable                 Name                                         Description
=============================================  ====================  ===============================================  ==============================================================================
Total capital cost                              :math:`C_{ca,tot}`    ``total_capital_cost``                          Total capital cost
Unit capital cost                               :math:`C_{ca,u}`      ``aggregate_capital_cost``                      Unit processes capital cost
Total operating cost                            :math:`C_{op,tot}`    ``total_operating_cost``                        Total operating cost for unit process
Maintenance, labor, and chemical cost           :math:`C_{mlc}`       ``maintenance_labor_chemical_operating_cost``   Maintenance cost for unit process
Total annualized cost                           :math:`C_{annual}`    ``total_annualized_costs``                      Total cost on a annualized basis
=============================================  ====================  ===============================================  ==============================================================================

Calculations for each of these costs are presented below.

Costing Index and Technoeconomic Factors
----------------------------------------

Default costing indices are provide with the WaterTAP Costing Package, but the user is free to modify these for their needs.

WaterTAP uses the CE (Chemical Engineering) Cost Index to help account for the time-value of investments and are used in the capital
and operating cost calculations. Unit process capital costs are adjusted to the year of the case study. The default year is 2018.

Technoeconomic factors used to calculate various system metrics, capital, and operating costs are presented in the table below:

=============================================  ====================  =======================================  ===============  ==============================================================================
                 Cost factor                     Variable                 Name                                 Default Value    Description
=============================================  ====================  =======================================  ===============  ==============================================================================
Plant capacity utilization factor                 :math:`f_{util}`    ``utilization_factor``                   90%             Percentage of year plant is operating
Total investment factor                           :math:`f_{toti}`    ``factor_total_investment``              2.0             Total investment factor (investment cost / equipment cost)
Maintenance-labor-chemical factor                 :math:`f_{mlc}`     ``factor_maintenance_labor_chemical``    0.03            Maintenance, labor, and chemical factor (fraction of investment cost / year)
Captial annualization factor                      :math:`f_{caf}`     ``factor_capital_annualization``         0.1             Capital annualization factor (fraction of investment cost / year) 
Capital recovery factor                           :math:`f_{crf}`     ``capital_recovery_factor``              0.1             Identical to `factor_capital_annualization`
=============================================  ====================  =======================================  ===============  ==============================================================================

Costing Calculations
--------------------

Total annulized cost is a simple function of the annualized capital cost and the annualized operating cost:

    .. math::
    
        C_{annual} = f_{crf} C_{ca,tot} + C_{op,tot}

The total capital cost is a simple factor of the sum of the unit model capital costs:

    .. math::

        C_{ca,tot} = f_{toti} C_{ca,u}

The operating cost :math:`C_{op,tot}` is the sum of the maintence, labor, and chemical operating costs, :math:`C_{mlc}`, the total fixed operating costs from the unit models, :math:`C_{fop,u}`, the total variable operating cost from the unit models, :math:`C_{vop,u}`, and the sum of the flow costs, :math:`C_{flow,tot}`:

   .. math::

        C_{op,totl} = C_{mlc} + C_{fop,u} + C_{vop,u} + f_{util} C_{flow,tot}

Where the maintenance-labor-chemical operating cost :math:`C_{mlc}` is defined as:

   .. math::

        C_{mlc} = f_{mlc} C_{ca,tot}

Other aggregates, like levelized cost of water (LCOW), are provided through the :ref:`WaterTAPCostingBlockData`: :ref:`technical_reference/costing/costing_base:Aggregates Metrics`. 

Class Documentation
-------------------

* :class:`WaterTAPCostingData`
