###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

from collections import MutableMapping

import pyomo.environ as pyo

from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core import declare_process_block_class
from idaes.core.base.costing_base import (
    FlowsheetCostingBlockData,
    register_idaes_currency_units,
)

from idaes.models.unit_models import Mixer

from watertap.unit_models import (
    ReverseOsmosis0D,
    ReverseOsmosis1D,
    NanoFiltration0D,
    NanofiltrationZO,
    PressureExchanger,
    Crystallization,
    Ultraviolet0D,
    Pump,
    EnergyRecoveryDevice,
    Electrodialysis0D,
    Electrodialysis1D,
    IonExchange0D,
    GAC,
)

from .units.crystallizer import cost_crystallizer
from .units.electrodialysis import cost_electrodialysis
from .units.energy_recovery_device import cost_energy_recovery_device
from .units.gac import cost_gac
from .units.ion_exchange import cost_ion_exchange
from .units.nanofiltration import cost_nanofiltration
from .units.mixer import cost_mixer
from .units.pressure_exchanger import cost_pressure_exchanger
from .units.pump import cost_pump
from .units.reverse_osmosis import cost_reverse_osmosis
from .units.uv_aop import cost_uv_aop


class _DefinedFlowsDict(MutableMapping, dict):
    # use dict methods
    __getitem__ = dict.__getitem__
    __iter__ = dict.__iter__
    __len__ = dict.__len__

    def __setitem__(self, key, value):
        if key in self:
            raise KeyError(f"{key} has already been defined as a flow")
        dict.__setitem__(self, key, value)

    def __delitem__(self, key):
        raise KeyError("defined flows cannot be removed")


@declare_process_block_class("WaterTAPCosting")
class WaterTAPCostingData(FlowsheetCostingBlockData):
    # Define default mapping of costing methods to unit models
    unit_mapping = {
        Mixer: cost_mixer,
        Pump: cost_pump,
        EnergyRecoveryDevice: cost_energy_recovery_device,
        PressureExchanger: cost_pressure_exchanger,
        ReverseOsmosis0D: cost_reverse_osmosis,
        ReverseOsmosis1D: cost_reverse_osmosis,
        NanoFiltration0D: cost_nanofiltration,
        NanofiltrationZO: cost_nanofiltration,
        Crystallization: cost_crystallizer,
        Ultraviolet0D: cost_uv_aop,
        Electrodialysis0D: cost_electrodialysis,
        Electrodialysis1D: cost_electrodialysis,
        IonExchange0D: cost_ion_exchange,
        GAC: cost_gac,
    }

    def build(self):
        super().build()
        self._registered_LCOWs = {}

    def build_global_params(self):

        # Register currency and conversion rates based on CE Index
        register_idaes_currency_units()

        # Set the base year for all costs
        self.base_currency = pyo.units.USD_2018
        # Set a base period for all operating costs
        self.base_period = pyo.units.year

        # Define standard material flows and costs
        # The WaterTAP costing package creates flows
        # in a lazy fashion, the first time `cost_flow`
        # is called for a flow. The `_DefinedFlowsDict`
        # prevents defining more than one flow with
        # the same name.
        self.defined_flows = _DefinedFlowsDict()

        # Build flowsheet level costing components
        # These are the global parameters
        self.utilization_factor = pyo.Var(
            initialize=0.9,
            doc="Plant capacity utilization [fraction of uptime]",
            units=pyo.units.dimensionless,
        )
        self.factor_total_investment = pyo.Var(
            initialize=2,
            doc="Total investment factor [investment cost/equipment cost]",
            units=pyo.units.dimensionless,
        )
        self.factor_maintenance_labor_chemical = pyo.Var(
            initialize=0.03,
            doc="Maintenance-labor-chemical factor [fraction of investment cost/year]",
            units=pyo.units.year**-1,
        )
        self.factor_capital_annualization = pyo.Var(
            initialize=0.1,
            doc="Capital annualization factor [fraction of investment cost/year]",
            units=pyo.units.year**-1,
        )

        self.electricity_base_cost = pyo.Param(
            mutable=True,
            initialize=0.07,
            doc="Electricity cost",
            units=pyo.units.USD_2018 / pyo.units.kWh,
        )
        self.defined_flows["electricity"] = self.electricity_base_cost

        self.electrical_carbon_intensity = pyo.Param(
            mutable=True,
            initialize=0.475,
            doc="Grid carbon intensity [kgCO2_eq/kWh]",
            units=pyo.units.kg / pyo.units.kWh,
        )

        # Crystallizer operating cost information from literature
        self.steam_unit_cost = pyo.Var(
            initialize=0.004,
            units=pyo.units.USD_2018 / (pyo.units.meter**3),
            doc="Steam cost, Panagopoulos (2019)",
        )
        self.defined_flows["steam"] = self.steam_unit_cost

        def build_naocl_cost_param_block(blk):

            blk.cost = pyo.Param(
                initialize=0.23,
                doc="NaOCl cost",
                units=pyo.units.USD_2018 / pyo.units.kg,
            )
            blk.purity = pyo.Param(
                mutable=True,
                initialize=0.15,
                doc="NaOCl purity",
                units=pyo.units.dimensionless,
            )

        self.naocl = pyo.Block(rule=build_naocl_cost_param_block)
        self.defined_flows["NaOCl"] = self.naocl.cost / self.naocl.purity

        def build_caoh2_cost_param_block(blk):
            blk.cost = pyo.Param(
                mutable=True,
                initialize=0.12,
                doc="CaOH2 cost",
                units=pyo.units.USD_2018 / pyo.units.kg,
            )
            blk.purity = pyo.Param(
                mutable=True,
                initialize=1,
                doc="CaOH2 purity",
                units=pyo.units.dimensionless,
            )

        self.caoh2 = pyo.Block(rule=build_caoh2_cost_param_block)
        self.defined_flows["CaOH2"] = self.caoh2.cost / self.caoh2.purity

        def build_hcl_cost_param_block(blk):

            blk.cost = pyo.Param(
                mutable=True,
                initialize=0.17,
                doc="HCl cost",  # for 37% sol'n - CatCost v 1.0.4
                units=pyo.units.USD_2020 / pyo.units.kg,
            )
            blk.purity = pyo.Param(
                mutable=True,
                initialize=0.37,
                doc="HCl purity",
                units=pyo.units.dimensionless,
            )

        self.hcl = pyo.Block(rule=build_hcl_cost_param_block)
        self.defined_flows["HCl"] = self.hcl.cost / self.hcl.purity

        def build_naoh_cost_param_block(blk):

            blk.cost = pyo.Param(
                mutable=True,
                initialize=0.59,
                doc="NaOH cost",  # for 30% sol'n - iDST
                units=pyo.units.USD_2020 / pyo.units.kg,
            )

            blk.purity = pyo.Param(
                mutable=True,
                initialize=0.30,
                doc="NaOH purity",
                units=pyo.units.dimensionless,
            )

        self.naoh = pyo.Block(rule=build_naoh_cost_param_block)
        self.defined_flows["NaOH"] = self.naoh.cost / self.naoh.purity

        def build_meoh_cost_param_block(blk):
            # MeOH = Methanol
            blk.cost = pyo.Param(
                mutable=True,
                initialize=3.395,
                doc="MeOH cost",  # for 100% purity - ICIS
                units=pyo.units.USD_2008 / pyo.units.kg,
            )

            blk.purity = pyo.Param(
                mutable=True,
                initialize=1,
                doc="MeOH purity",
                units=pyo.units.dimensionless,
            )

        self.meoh = pyo.Block(rule=build_meoh_cost_param_block)
        self.defined_flows["MeOH"] = self.meoh.cost / self.meoh.purity

        def build_nacl_cost_param_block(blk):

            blk.cost = pyo.Param(
                mutable=True,
                initialize=0.09,
                doc="NaCl cost",  # for solid, 100% purity - CatCost
                units=pyo.units.USD_2020 / pyo.units.kg,
            )

            blk.purity = pyo.Param(
                mutable=True,
                initialize=1,
                doc="NaCl purity",
                units=pyo.units.dimensionless,
            )

        self.nacl = pyo.Block(rule=build_nacl_cost_param_block)
        self.defined_flows["NaCl"] = self.nacl.cost / self.nacl.purity

        # fix the parameters
        for var in self.component_objects(pyo.Var, descend_into=True):
            var.fix()

    def cost_flow(self, flow_expr, flow_type):
        """
        This method registers a given flow component (Var or expression) for
        costing. All flows are required to be bounded to be non-negative (i.e.
        a lower bound equal to or greater than 0).

        Args:
            flow_expr: Pyomo Var or expression that represents a material flow
                that should be included in the process costing. Units are
                expected to be on a per time basis.
            flow_type: string identifying the material this flow represents.
                This string must be available to the FlowsheetCostingBlock
                as a known flow type.

        Raises:
            ValueError if flow_type is not recognized.
            TypeError if flow_expr is an indexed Var.
        """
        if flow_type not in self.defined_flows:
            raise ValueError(
                f"{flow_type} is not a recognized flow type. Please check "
                "your spelling and that the flow type has been available to"
                " the FlowsheetCostingBlock."
            )
        if flow_type not in self.flow_types:
            self.register_flow_type(flow_type, self.defined_flows[flow_type])
        super().cost_flow(flow_expr, flow_type)

    def build_process_costs(self):
        self.total_capital_cost = pyo.Expression(
            expr=self.aggregate_capital_cost, doc="Total capital cost"
        )
        self.total_investment_cost = pyo.Var(
            initialize=1e3,
            domain=pyo.NonNegativeReals,
            doc="Total investment cost",
            units=self.base_currency,
        )
        self.maintenance_labor_chemical_operating_cost = pyo.Var(
            initialize=1e3,
            domain=pyo.NonNegativeReals,
            doc="Maintenance-labor-chemical operating cost",
            units=self.base_currency / self.base_period,
        )
        self.total_operating_cost = pyo.Var(
            initialize=1e3,
            domain=pyo.NonNegativeReals,
            doc="Total operating cost",
            units=self.base_currency / self.base_period,
        )

        self.total_investment_cost_constraint = pyo.Constraint(
            expr=self.total_investment_cost
            == self.factor_total_investment * self.total_capital_cost
        )
        self.maintenance_labor_chemical_operating_cost_constraint = pyo.Constraint(
            expr=self.maintenance_labor_chemical_operating_cost
            == self.factor_maintenance_labor_chemical * self.total_investment_cost
        )

        self.total_operating_cost_constraint = pyo.Constraint(
            expr=self.total_operating_cost
            == self.maintenance_labor_chemical_operating_cost
            + self.aggregate_fixed_operating_cost
            + self.aggregate_variable_operating_cost
            + sum(self.aggregate_flow_costs.values()) * self.utilization_factor
        )

    def initialize_build(self):
        calculate_variable_from_constraint(
            self.total_investment_cost, self.total_investment_cost_constraint
        )
        calculate_variable_from_constraint(
            self.maintenance_labor_chemical_operating_cost,
            self.maintenance_labor_chemical_operating_cost_constraint,
        )
        calculate_variable_from_constraint(
            self.total_operating_cost, self.total_operating_cost_constraint
        )

        for var, con in self._registered_LCOWs.values():
            calculate_variable_from_constraint(var, con)

    def add_LCOW(self, flow_rate, name="LCOW"):
        """
        Add Levelized Cost of Water (LCOW) to costing block.
        Args:
            flow_rate - flow rate of water (volumetric) to be used in
                        calculating LCOW
            name (optional) - name for the LCOW variable (default: LCOW)
        """

        LCOW = pyo.Var(
            doc=f"Levelized Cost of Water based on flow {flow_rate.name}",
            units=self.base_currency / pyo.units.m**3,
        )
        self.add_component(name, LCOW)

        LCOW_constraint = pyo.Constraint(
            expr=LCOW
            == (
                self.total_investment_cost * self.factor_capital_annualization
                + self.total_operating_cost
            )
            / (
                pyo.units.convert(
                    flow_rate, to_units=pyo.units.m**3 / self.base_period
                )
                * self.utilization_factor
            ),
            doc=f"Constraint for Levelized Cost of Water based on flow {flow_rate.name}",
        )
        self.add_component(name + "_constraint", LCOW_constraint)

        self._registered_LCOWs[name] = (LCOW, LCOW_constraint)

    def add_annual_water_production(self, flow_rate, name="annual_water_production"):
        """
        Add annual water production to costing block.
        Args:
            flow_rate - flow rate of water (volumetric) to be used in
                        calculating annual water production
            name (optional) - name for the annual water productionvariable
                              Expression (default: annual_water_production)
        """
        self.add_component(
            name,
            pyo.Expression(
                expr=(
                    pyo.units.convert(
                        flow_rate, to_units=pyo.units.m**3 / self.base_period
                    )
                    * self.utilization_factor
                ),
                doc=f"Annual water production based on flow {flow_rate.name}",
            ),
        )

    def add_specific_energy_consumption(
        self, flow_rate, name="specific_energy_consumption"
    ):
        """
        Add specific energy consumption (kWh/m**3) to costing block.
        Args:
            flow_rate - flow rate of water (volumetric) to be used in
                        calculating specific energy consumption
            name (optional) - the name of the Expression for the specific
                              energy consumption (default: specific_energy_consumption)
        """

        self.add_component(
            name,
            pyo.Expression(
                expr=self.aggregate_flow_electricity
                / pyo.units.convert(
                    flow_rate, to_units=pyo.units.m**3 / pyo.units.hr
                ),
                doc=f"Specific energy consumption based on flow {flow_rate.name}",
            ),
        )

    def add_specific_electrical_carbon_intensity(
        self, flow_rate, name="specific_electrical_carbon_intensity"
    ):
        """
        Add specific electrical carbon intensity (kg_CO2eq/m**3) to costing block.
        Args:
            flow_rate - flow rate of water (volumetric) to be used in
                        calculating specific electrical carbon intensity
            name (optional) - the name of the Expression for the specific
                              energy consumption (default: specific_electrical_carbon_intensity)
        """

        self.add_component(
            name,
            pyo.Expression(
                expr=self.aggregate_flow_electricity
                * self.electrical_carbon_intensity
                / pyo.units.convert(
                    flow_rate, to_units=pyo.units.m**3 / pyo.units.hr
                ),
                doc=f"Specific electrical carbon intensity based on flow {flow_rate.name}",
            ),
        )
