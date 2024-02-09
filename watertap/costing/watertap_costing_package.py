#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pyomo.environ as pyo

from idaes.core import declare_process_block_class
from idaes.core.base.costing_base import register_idaes_currency_units

from watertap.costing.costing_base import WaterTAPCostingBlockData


@declare_process_block_class("WaterTAPCosting")
class WaterTAPCostingData(WaterTAPCostingBlockData):
    def build_global_params(self):
        # Register currency and conversion rates based on CE Index
        register_idaes_currency_units()

        self._build_common_global_params()

        # Set the base year for all costs
        self.base_currency = pyo.units.USD_2018
        # Set a base period for all operating costs
        self.base_period = pyo.units.year

        # Build flowsheet level costing components
        # These are the global parameters
        self.factor_total_investment = pyo.Var(
            initialize=1.0,
            doc="Total investment factor [investment cost/equipment cost]",
            units=pyo.units.dimensionless,
        )
        self.factor_maintenance_labor_chemical = pyo.Var(
            initialize=0.03,
            doc="Maintenance-labor-chemical factor [fraction of equipment cost/year]",
            units=pyo.units.year**-1,
        )

        # fix the parameters
        self.fix_all_vars()

        # keep wacc floating by default
        # so that factor_capital_annualization
        # can be fixed
        self.wacc.unfix()
