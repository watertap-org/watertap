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

from enum import Enum

import pyomo.environ as pyo

from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core import declare_process_block_class
from idaes.core.base.costing_base import (
    FlowsheetCostingBlockData,
    register_idaes_currency_units,
)
from idaes.core.util.constants import Constants
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.misc import StrEnum
from idaes.core.util.math import smooth_min

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


class ROType(StrEnum):
    standard = "standard"
    high_pressure = "high_pressure"


class PumpType(StrEnum):
    low_pressure = "low_pressure"
    high_pressure = "high_pressure"


class EnergyRecoveryDeviceType(StrEnum):
    pressure_exchanger = "pressure_exchanger"


class MixerType(StrEnum):
    default = "default"
    NaOCl = "NaOCl"
    CaOH2 = "CaOH2"


class CrystallizerCostType(StrEnum):
    default = "default"
    mass_basis = "mass_basis"
    volume_basis = "volume_basis"


@declare_process_block_class("WaterTAPCosting")
class WaterTAPCostingData(FlowsheetCostingBlockData):
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
        # is called for a flow. The `available_flows`
        # are all the flows which can be costed with
        # the WaterTAP costing package.
        self.available_flows = {}

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
            units=self.base_currency / pyo.units.kWh,
        )
        self.available_flows["electricity"] = self.electricity_base_cost

        self.electrical_carbon_intensity = pyo.Param(
            mutable=True,
            initialize=0.475,
            doc="Grid carbon intensity [kgCO2_eq/kWh]",
            units=pyo.units.kg / pyo.units.kWh,
        )

        def build_reverse_osmosis_cost_param_block(blk):

            costing = blk.parent_block()

            blk.factor_membrane_replacement = pyo.Var(
                initialize=0.2,
                doc="Membrane replacement factor [fraction of membrane replaced/year]",
                units=pyo.units.year**-1,
            )
            blk.membrane_cost = pyo.Var(
                initialize=30,
                doc="Membrane cost",
                units=costing.base_currency / (pyo.units.meter**2),
            )
            blk.high_pressure_membrane_cost = pyo.Var(
                initialize=75,
                doc="Membrane cost",
                units=costing.base_currency / (pyo.units.meter**2),
            )

        self.reverse_osmosis = pyo.Block(rule=build_reverse_osmosis_cost_param_block)

        def build_nanofiltration_cost_param_block(blk):

            costing = blk.parent_block()

            blk.factor_membrane_replacement = pyo.Var(
                initialize=0.2,
                doc="Membrane replacement factor [fraction of membrane replaced/year]",
                units=pyo.units.year**-1,
            )
            blk.membrane_cost = pyo.Var(
                initialize=15,
                doc="Membrane cost",
                units=costing.base_currency / (pyo.units.meter**2),
            )

        self.nanofiltration = pyo.Block(rule=build_nanofiltration_cost_param_block)

        def build_uv_cost_param_block(blk):

            costing = blk.parent_block()

            blk.factor_lamp_replacement = pyo.Var(
                initialize=0.33278,
                doc="UV replacement factor accounting for lamps, sleeves, ballasts and sensors [fraction of uv replaced/year]",
                units=pyo.units.year**-1,
            )
            blk.reactor_cost = pyo.Var(
                initialize=202.346,
                doc="UV reactor cost",
                units=costing.base_currency / (pyo.units.m**3 / pyo.units.hr),
            )
            blk.lamp_cost = pyo.Var(
                initialize=235.5,
                doc="UV lamps, sleeves, ballasts and sensors cost",
                units=costing.base_currency / pyo.units.kW,
            )

        self.ultraviolet = pyo.Block(rule=build_uv_cost_param_block)

        def build_high_pressure_pump_cost_param_block(blk):

            costing = blk.parent_block()

            blk.cost = pyo.Var(
                initialize=53 / 1e5 * 3600,
                doc="High pressure pump cost",
                units=costing.base_currency / pyo.units.watt,
            )

        self.high_pressure_pump = pyo.Block(
            rule=build_high_pressure_pump_cost_param_block
        )

        def build_low_pressure_pump_cost_param_block(blk):

            costing = blk.parent_block()

            blk.cost = pyo.Var(
                initialize=889,
                doc="Low pressure pump cost",
                units=costing.base_currency / (pyo.units.liter / pyo.units.second),
            )

        self.low_pressure_pump = pyo.Block(
            rule=build_low_pressure_pump_cost_param_block
        )

        def build_energy_recovery_device_cost_param_block(blk):

            costing = blk.parent_block()

            blk.pressure_exchanger_cost = pyo.Var(
                initialize=535,
                doc="Pressure exchanger cost",
                units=costing.base_currency / (pyo.units.meter**3 / pyo.units.hours),
            )

        self.energy_recovery_device = pyo.Block(
            rule=build_energy_recovery_device_cost_param_block
        )

        def build_pressure_exchanger_cost_param_block(blk):

            costing = blk.parent_block()

            blk.cost = pyo.Var(
                initialize=535,
                doc="Pressure exchanger cost",
                units=costing.base_currency / (pyo.units.meter**3 / pyo.units.hours),
            )

        self.pressure_exchanger = pyo.Block(
            rule=build_pressure_exchanger_cost_param_block
        )

        def build_mixer_cost_param_block(blk):

            costing = blk.parent_block()

            blk.unit_cost = pyo.Var(
                initialize=361,
                doc="Mixer cost",
                units=costing.base_currency / (pyo.units.liters / pyo.units.second),
            )

        self.mixer = pyo.Block(rule=build_mixer_cost_param_block)

        def build_naocl_cost_param_block(blk):

            costing = blk.parent_block()

            blk.mixer_unit_cost = pyo.Var(
                initialize=5.08,
                doc="NaOCl mixer cost",
                units=costing.base_currency / (pyo.units.m**3 / pyo.units.day),
            )
            blk.cost = pyo.Param(
                initialize=0.23,
                doc="NaOCl cost",
                units=costing.base_currency / pyo.units.kg,
            )
            blk.purity = pyo.Param(
                mutable=True,
                initialize=0.15,
                doc="NaOCl purity",
                units=pyo.units.dimensionless,
            )
            costing.available_flows["NaOCl"] = blk.cost / blk.purity

        self.naocl = pyo.Block(rule=build_naocl_cost_param_block)

        def build_caoh2_cost_param_block(blk):

            costing = blk.parent_block()

            blk.mixer_unit_cost = pyo.Var(
                initialize=792.8 * 2.20462,
                doc="Ca(OH)2 mixer cost",
                units=costing.base_currency / (pyo.units.kg / pyo.units.day),
            )
            blk.cost = pyo.Param(
                mutable=True,
                initialize=0.12,
                doc="CaOH2 cost",
                units=costing.base_currency / pyo.units.kg,
            )
            blk.purity = pyo.Param(
                mutable=True,
                initialize=1,
                doc="CaOH2 purity",
                units=pyo.units.dimensionless,
            )
            costing.available_flows["CaOH2"] = blk.cost / blk.purity

        self.caoh2 = pyo.Block(rule=build_caoh2_cost_param_block)

        def build_hcl_cost_param_block(blk):

            costing = blk.parent_block()

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
            costing.available_flows["HCl"] = blk.cost / blk.purity

        self.hcl = pyo.Block(rule=build_hcl_cost_param_block)

        def build_naoh_cost_param_block(blk):

            costing = blk.parent_block()

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
            costing.available_flows["NaOH"] = blk.cost / blk.purity

        self.naoh = pyo.Block(rule=build_naoh_cost_param_block)

        def build_meoh_cost_param_block(blk):
            # MeOH = Methanol

            costing = blk.parent_block()

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
            costing.available_flows["MeOH"] = blk.cost / blk.purity

        self.meoh = pyo.Block(rule=build_meoh_cost_param_block)

        def build_nacl_cost_param_block(blk):

            costing = blk.parent_block()

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
            costing.available_flows["NaCl"] = blk.cost / blk.purity

        self.nacl = pyo.Block(rule=build_nacl_cost_param_block)

        def build_electrodialysis_cost_param_block(blk):

            costing = blk.parent_block()

            blk.cem_membrane_cost = pyo.Var(
                initialize=43,
                doc="Cost of CEM membrane used in Electrodialysis ($/CEM/area)",
                units=costing.base_currency / (pyo.units.meter**2),
            )
            blk.aem_membrane_cost = pyo.Var(
                initialize=43,
                doc="Cost of AEM membrane used in Electrodialysis ($/AEM/area)",
                units=costing.base_currency / (pyo.units.meter**2),
            )
            blk.flowspacer_cost = pyo.Var(
                initialize=3,
                doc="Cost of the spacers used in Electrodialysis ($/spacer/area)",
                units=costing.base_currency / (pyo.units.meter**2),
            )
            blk.factor_membrane_housing_replacement = pyo.Var(
                initialize=0.2,
                doc="Membrane housing replacement factor for CEM, AEM, and spacer replacements [fraction of membrane replaced/year]",
                units=pyo.units.year**-1,
            )
            blk.electrode_cost = pyo.Var(
                initialize=2000,
                doc="Cost of the electrodes used in Electrodialysis ($/electrode/area)",
                units=costing.base_currency / (pyo.units.meter**2),
            )
            blk.factor_electrode_replacement = pyo.Var(
                initialize=0.02,
                doc="Electrode replacements [fraction of electrode replaced/year]",
                units=pyo.units.year**-1,
            )

        self.electrodialysis = pyo.Block(rule=build_electrodialysis_cost_param_block)

        def build_crystallizer_cost_param_block(blk):

            blk.steam_pressure = pyo.Var(
                initialize=3,
                units=pyo.units.bar,
                doc="Steam pressure (gauge) for crystallizer heating: 3 bar default based on Dutta example",
            )

            blk.efficiency_pump = pyo.Var(
                initialize=0.7,
                units=pyo.units.dimensionless,
                doc="Crystallizer pump efficiency - assumed",
            )

            blk.pump_head_height = pyo.Var(
                initialize=1,
                units=pyo.units.m,
                doc="Crystallizer pump head height -  assumed, unvalidated",
            )

            # Crystallizer operating cost information from literature
            blk.fob_unit_cost = pyo.Var(
                initialize=675000,
                doc="Forced circulation crystallizer reference free-on-board cost (Woods, 2007)",
                units=pyo.units.USD_2007,
            )

            blk.ref_capacity = pyo.Var(
                initialize=1,
                doc="Forced circulation crystallizer reference crystal capacity (Woods, 2007)",
                units=pyo.units.kg / pyo.units.s,
            )

            blk.ref_exponent = pyo.Var(
                initialize=0.53,
                doc="Forced circulation crystallizer cost exponent factor (Woods, 2007)",
                units=pyo.units.dimensionless,
            )

            blk.iec_percent = pyo.Var(
                initialize=1.43,
                doc="Forced circulation crystallizer installed equipment cost (Diab and Gerogiorgis, 2017)",
                units=pyo.units.dimensionless,
            )

            blk.volume_cost = pyo.Var(
                initialize=16320,
                doc="Forced circulation crystallizer cost per volume (Yusuf et al., 2019)",
                units=pyo.units.USD_2007,  ## TODO: Needs confirmation, but data is from Perry apparently
            )

            blk.vol_basis_exponent = pyo.Var(
                initialize=0.47,
                doc="Forced circulation crystallizer volume-based cost exponent (Yusuf et al., 2019)",
                units=pyo.units.dimensionless,
            )

        self.crystallizer = pyo.Block(rule=build_crystallizer_cost_param_block)

        # Crystallizer operating cost information from literature
        self.steam_unit_cost = pyo.Var(
            initialize=0.004,
            units=pyo.units.USD_2018 / (pyo.units.meter**3),
            doc="Steam cost, Panagopoulos (2019)",
        )
        self.available_flows["steam"] = self.steam_unit_cost

        def build_ion_exhange_cost_param_block(blk):
            blk.anion_exchange_resin_cost = pyo.Var(
                initialize=205,
                units=pyo.units.USD_2020 / pyo.units.ft**3,
                doc="Anion exchange resin cost per cubic ft. Assumes strong base polystyrenic gel-type Type II. From EPA-WBS cost model.",
            )
            blk.cation_exchange_resin_cost = pyo.Var(
                initialize=153,
                units=pyo.units.USD_2020 / pyo.units.ft**3,
                doc="Cation exchange resin cost per cubic ft. Assumes strong acid polystyrenic gel-type. From EPA-WBS cost model.",
            )
            # Ion exchange pressure vessels costed with 3rd order polynomial:
            #   pv_cost = A * col_vol^3 + B * col_vol^2 + C * col_vol + intercept

            blk.vessel_intercept = pyo.Var(
                initialize=10010.86,
                units=pyo.units.USD_2020,
                doc="Ion exchange pressure vessel cost equation - intercept, Carbon steel w/ plastic internals",
            )
            blk.vessel_A_coeff = pyo.Var(
                initialize=6e-9,
                units=pyo.units.USD_2020 / pyo.units.gal**3,
                doc="Ion exchange pressure vessel cost equation - A coeff., Carbon steel w/ plastic internals",
            )
            blk.vessel_B_coeff = pyo.Var(
                initialize=-2.284e-4,
                units=pyo.units.USD_2020 / pyo.units.gal**2,
                doc="Ion exchange pressure vessel cost equation - B coeff., Carbon steel w/ plastic internals",
            )
            blk.vessel_C_coeff = pyo.Var(
                initialize=8.3472,
                units=pyo.units.USD_2020 / pyo.units.gal,
                doc="Ion exchange pressure vessel cost equation - C coeff., Carbon steel w/ plastic internals",
            )
            # Ion exchange pressure vessels costed with 3rd order polynomial:
            #   pv_cost = A * col_vol^3 + B * col_vol^2 + C * col_vol + intercept

            blk.backwash_tank_A_coeff = pyo.Var(
                initialize=1e-9,
                units=pyo.units.USD_2020 / pyo.units.gal**3,
                doc="Ion exchange backwash tank cost equation - A coeff., Fiberglass tank",
            )
            blk.backwash_tank_B_coeff = pyo.Var(
                initialize=-5.8587e-05,
                units=pyo.units.USD_2020 / pyo.units.gal**2,
                doc="Ion exchange backwash tank cost equation - B coeff., Fiberglass tank",
            )
            blk.backwash_tank_C_coeff = pyo.Var(
                initialize=2.2911,
                units=pyo.units.USD_2020 / pyo.units.gal,
                doc="Ion exchange backwash tank cost equation - C coeff., Fiberglass tank",
            )
            blk.backwash_tank_intercept = pyo.Var(
                initialize=4717.255,
                units=pyo.units.dimensionless,
                doc="Ion exchange backwash tank cost equation - exponent, Fiberglass tank",
            )
            # Ion exchange regeneration solution tank costed with 2nd order polynomial:
            #   regen_tank_cost = A * tank_vol^2 + B * tank_vol + intercept

            blk.regen_tank_intercept = pyo.Var(
                initialize=4408.327,
                units=pyo.units.USD_2020 / pyo.units.gal,
                doc="Ion exchange regen tank cost equation - C coeff. Stainless steel",
            )
            blk.regen_tank_A_coeff = pyo.Var(
                initialize=-3.258e-5,
                units=pyo.units.USD_2020 / pyo.units.gal**2,
                doc="Ion exchange regen tank cost equation - A coeff. Stainless steel",
            )
            blk.regen_tank_B_coeff = pyo.Var(
                initialize=3.846,
                units=pyo.units.USD_2020 / pyo.units.gal,
                doc="Ion exchange regen tank cost equation - B coeff. Stainless steel",
            )
            blk.annual_resin_replacement_factor = pyo.Var(
                initialize=0.05,
                units=pyo.units.year**-1,
                doc="Fraction of ion excange resin replaced per year, 4-5% of bed volume - EPA",
            )
            blk.hazardous_min_cost = pyo.Var(
                initialize=3240,
                units=pyo.units.USD_2020,
                doc="Min cost per hazardous waste shipment - EPA",
            )
            blk.hazardous_resin_disposal = pyo.Var(
                initialize=347.10,
                units=pyo.units.USD_2020 * pyo.units.ton**-1,
                doc="Hazardous resin disposal cost - EPA",
            )
            blk.hazardous_regen_disposal = pyo.Var(
                initialize=3.64,
                units=pyo.units.USD_2020 * pyo.units.gal**-1,
                doc="Hazardous liquid disposal cost - EPA",
            )

        self.ion_exchange = pyo.Block(rule=build_ion_exhange_cost_param_block)

        def build_gac_cost_param_block(blk):

            blk.num_contactors_op = pyo.Var(
                initialize=1,
                units=pyo.units.dimensionless,
                doc="Number of GAC contactors in operation in parallel",
            )

            blk.num_contactors_redundant = pyo.Var(
                initialize=1,
                units=pyo.units.dimensionless,
                doc="Number of off-line redundant GAC contactors in parallel",
            )

            blk.contactor_cost_coeff_0 = pyo.Var(
                initialize=10010.9,
                units=pyo.units.USD_2020,
                doc="GAC contactor polynomial cost coefficient 0",
            )

            blk.contactor_cost_coeff_1 = pyo.Var(
                initialize=2204.95,
                units=pyo.units.USD_2020 * (pyo.units.m**3) ** -1,
                doc="GAC contactor polynomial cost coefficient 1",
            )

            blk.contactor_cost_coeff_2 = pyo.Var(
                initialize=-15.9378,
                units=pyo.units.USD_2020 * (pyo.units.m**3) ** -2,
                doc="GAC contactor polynomial cost coefficient 2",
            )

            blk.contactor_cost_coeff_3 = pyo.Var(
                initialize=0.110592,
                units=pyo.units.USD_2020 * (pyo.units.m**3) ** -3,
                doc="GAC contactor polynomial cost coefficient 3",
            )

            blk.bed_mass_max_ref = pyo.Var(
                initialize=18143.7,
                units=pyo.units.kg,
                doc="Reference maximum value of GAC mass needed for initial charge where "
                "economy of scale no longer discounts the unit price",
            )

            blk.adsorbent_unit_cost_coeff = pyo.Var(
                initialize=4.58342,
                units=pyo.units.USD_2020 * pyo.units.kg**-1,
                doc="GAC adsorbent exponential cost pre-exponential coefficient",
            )

            blk.adsorbent_unit_cost_exp_coeff = pyo.Var(
                initialize=-1.25311e-5,
                units=pyo.units.kg**-1,
                doc="GAC adsorbent exponential cost parameter coefficient",
            )

            blk.other_cost_coeff = pyo.Var(
                initialize=16660.7,
                units=pyo.units.USD_2020,
                doc="GAC other cost power law coefficient",
            )

            blk.other_cost_exp = pyo.Var(
                initialize=0.552207,
                units=pyo.units.dimensionless,
                doc="GAC other cost power law exponent",
            )

            blk.regen_frac = pyo.Var(
                initialize=0.70,
                units=pyo.units.dimensionless,
                doc="Fraction of spent GAC adsorbent that can be regenerated for reuse",
            )

            blk.regen_unit_cost = pyo.Var(
                initialize=4.28352,
                units=pyo.units.USD_2020 * pyo.units.kg**-1,
                doc="Unit cost to regenerate spent GAC adsorbent by an offsite regeneration facility",
            )

            blk.makeup_unit_cost = pyo.Var(
                initialize=4.58223,
                units=pyo.units.USD_2020 * pyo.units.kg**-1,
                doc="Unit cost to makeup spent GAC adsorbent with fresh adsorbent",
            )

        self.gac = pyo.Block(rule=build_gac_cost_param_block)

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
        if flow_type not in self.available_flows:
            raise ValueError(
                f"{flow_type} is not a recognized flow type. Please check "
                "your spelling and that the flow type has been available to"
                " the FlowsheetCostingBlock."
            )
        if flow_type not in self.flow_types:
            self.register_flow_type(flow_type, self.available_flows[flow_type])
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

    # Define costing methods supported by package

    @staticmethod
    def cost_uv_aop(blk, cost_electricity_flow=True):
        """
        UV-AOP costing method
        """
        cost_uv_aop_bundle(
            blk,
            blk.costing_package.ultraviolet.reactor_cost,
            blk.costing_package.ultraviolet.lamp_cost,
            blk.costing_package.ultraviolet.factor_lamp_replacement,
        )

        t0 = blk.flowsheet().time.first()
        if cost_electricity_flow:
            blk.costing_package.cost_flow(
                pyo.units.convert(
                    blk.unit_model.electricity_demand[t0],
                    to_units=pyo.units.kW,
                ),
                "electricity",
            )

    @staticmethod
    def cost_nanofiltration(blk):
        """
        Nanofiltration costing method

        TODO: describe equations
        """
        cost_membrane(
            blk,
            blk.costing_package.nanofiltration.membrane_cost,
            blk.costing_package.nanofiltration.factor_membrane_replacement,
        )

    @staticmethod
    def cost_reverse_osmosis(blk, ro_type=ROType.standard):
        """
        Reverse osmosis costing method

        TODO: describe equations

        Args:
            ro_type: ROType Enum indicating reverse osmosis type,
                default = ROType.standard
        """
        if ro_type == ROType.standard:
            membrane_cost = blk.costing_package.reverse_osmosis.membrane_cost
        elif ro_type == ROType.high_pressure:
            membrane_cost = (
                blk.costing_package.reverse_osmosis.high_pressure_membrane_cost
            )
        else:
            raise ConfigurationError(
                f"{blk.unit_model.name} received invalid argument for ro_type:"
                f" {ro_type}. Argument must be a member of the ROType Enum."
            )
        cost_membrane(
            blk,
            membrane_cost,
            blk.costing_package.reverse_osmosis.factor_membrane_replacement,
        )

    @staticmethod
    def cost_pump(blk, pump_type=PumpType.high_pressure, cost_electricity_flow=True):
        """
        Pump costing method

        TODO: describe equations

        Args:
            pump_type: PumpType Enum indicating pump type,
                default = PumpType.high_pressure

            cost_electricity_flow: bool, if True, the Pump's work_mechanical will be
                converted to kW and costed as an electricity, default = True
        """
        if pump_type == PumpType.high_pressure:
            WaterTAPCostingData.cost_high_pressure_pump(blk, cost_electricity_flow)
        elif pump_type == PumpType.low_pressure:
            WaterTAPCostingData.cost_low_pressure_pump(blk, cost_electricity_flow)
        else:
            raise ConfigurationError(
                f"{blk.unit_model.name} received invalid argument for pump_type:"
                f" {pump_type}. Argument must be a member of the PumpType Enum."
            )

    @staticmethod
    def cost_energy_recovery_device(
        blk,
        energy_recovery_device_type=EnergyRecoveryDeviceType.pressure_exchanger,
        cost_electricity_flow=True,
    ):
        """
        Energy recovery device costing method

        TODO: describe equations

        Args:
            energy_recovery_device_type: EnergyRecoveryDeviceType Enum indicating ERD type,
                default = EnergyRecoveryDeviceType.pressure_exchanger.

            cost_electricity_flow: bool, if True, the ERD's work_mechanical will
                be converted to kW and costed as an electricity default = True
        """
        if energy_recovery_device_type == EnergyRecoveryDeviceType.pressure_exchanger:
            WaterTAPCostingData.cost_pressure_exchanger_erd(blk)
        else:
            raise ConfigurationError(
                f"{blk.unit_model.name} received invalid argument for energy_recovery_device_type:"
                f" {energy_recovery_device_type}. Argument must be a member of the EnergyRecoveryDeviceType Enum."
            )

    @staticmethod
    def cost_high_pressure_pump(blk, cost_electricity_flow=True):
        """
        High pressure pump costing method

        `TODO: describe equations`

        Args:
            cost_electricity_flow - bool, if True, the Pump's work_mechanical will
                                    be converted to kW and costed as an electricity
                                    default = True
        """
        t0 = blk.flowsheet().time.first()
        make_capital_cost_var(blk)
        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost
            == blk.costing_package.high_pressure_pump.cost
            * pyo.units.convert(blk.unit_model.work_mechanical[t0], pyo.units.W)
        )
        if cost_electricity_flow:
            blk.costing_package.cost_flow(
                pyo.units.convert(
                    blk.unit_model.work_mechanical[t0], to_units=pyo.units.kW
                ),
                "electricity",
            )

    @staticmethod
    def cost_low_pressure_pump(blk, cost_electricity_flow=True):
        """
        Low pressure pump costing method

        TODO: describe equations

        Args:
            cost_electricity_flow - bool, if True, the Pump's work_mechanical will
                                    be converted to kW and costed as an electricity
                                    default = True
        """
        t0 = blk.flowsheet().time.first()
        cost_by_flow_volume(
            blk,
            blk.costing_package.low_pressure_pump.cost,
            pyo.units.convert(
                blk.unit_model.control_volume.properties_in[t0].flow_vol,
                (pyo.units.m**3 / pyo.units.s),
            ),
        )
        if cost_electricity_flow:
            blk.costing_package.cost_flow(
                pyo.units.convert(
                    blk.unit_model.work_mechanical[t0], to_units=pyo.units.kW
                ),
                "electricity",
            )

    @staticmethod
    def cost_pressure_exchanger_erd(blk, cost_electricity_flow=True):
        """
        ERD pressure exchanger costing method

        TODO: describe equations

        Args:
            cost_electricity_flow - bool, if True, the ERD's work_mechanical will
                                    be converted to kW and costed as an electricity
                                    default = True
        """
        t0 = blk.flowsheet().time.first()
        cost_by_flow_volume(
            blk,
            blk.costing_package.energy_recovery_device.pressure_exchanger_cost,
            pyo.units.convert(
                blk.unit_model.control_volume.properties_in[t0].flow_vol,
                (pyo.units.meter**3 / pyo.units.hours),
            ),
        )
        if cost_electricity_flow:
            blk.costing_package.cost_flow(
                pyo.units.convert(
                    blk.unit_model.work_mechanical[t0], to_units=pyo.units.kW
                ),
                "electricity",
            )

    @staticmethod
    def cost_pressure_exchanger(blk):
        """
        Pressure exchanger costing method

        TODO: describe equations
        """
        cost_by_flow_volume(
            blk,
            blk.costing_package.pressure_exchanger.cost,
            pyo.units.convert(
                blk.unit_model.low_pressure_side.properties_in[0].flow_vol,
                (pyo.units.meter**3 / pyo.units.hours),
            ),
        )

    @staticmethod
    def cost_mixer(blk, mixer_type=MixerType.default, **kwargs):
        """
        Mixer costing method

        TODO: describe equations

        Args:
            mixer_type: MixerType Enum indicating mixer type,
                default = MixerType.default

            `**kwargs`: Additional keywords for the MixerType, e.g., NaOCl
                and CaOH2 mixers expect the `dosing_rate` keyword
                argument.
        """
        if mixer_type == MixerType.default:
            WaterTAPCostingData.cost_default_mixer(blk, **kwargs)
        elif mixer_type == MixerType.NaOCl:
            WaterTAPCostingData.cost_naocl_mixer(blk, **kwargs)
        elif mixer_type == MixerType.CaOH2:
            WaterTAPCostingData.cost_caoh2_mixer(blk, **kwargs)
        else:
            raise ConfigurationError(
                f"{blk.unit_model.name} received invalid argument for mixer_type:"
                f" {mixer_type}. Argument must be a member of the MixerType Enum."
            )

    @staticmethod
    def cost_default_mixer(blk):
        """
        Default mixer costing method

        TODO: describe equations
        """
        cost_by_flow_volume(
            blk,
            blk.costing_package.mixer.unit_cost,
            pyo.units.convert(
                blk.unit_model.mixed_state[0].flow_vol,
                pyo.units.liter / pyo.units.second,
            ),
        )

    @staticmethod
    def cost_naocl_mixer(blk, dosing_rate):
        """
        NaOCl mixer costing method

        TODO: describe equations

        Args:
            dosing_rate: An expression in [mass/time] for NaOCl dosage
        """
        cost_by_flow_volume(
            blk,
            blk.costing_package.naocl.mixer_unit_cost,
            pyo.units.convert(
                blk.unit_model.inlet_stream_state[0].flow_vol,
                pyo.units.m**3 / pyo.units.day,
            ),
        )
        blk.costing_package.cost_flow(
            pyo.units.convert(dosing_rate, pyo.units.kg / pyo.units.s), "NaOCl"
        )

    @staticmethod
    def cost_caoh2_mixer(blk, dosing_rate):
        """
        CaOH2 mixer costing method

        TODO: describe equations

        Args:
            dosing_rate: An expression in [mass/time] for CaOH2 dosage
        """
        stream = blk.unit_model.lime_stream
        blk.lime_kg_per_day = pyo.Expression(
            expr=pyo.units.convert(
                dosing_rate,
                pyo.units.kg / pyo.units.day,
            )
        )
        cost_by_flow_volume(
            blk,
            blk.costing_package.caoh2.mixer_unit_cost
            / blk.costing_package.factor_total_investment,
            blk.lime_kg_per_day,
        )
        blk.costing_package.cost_flow(
            pyo.units.convert(dosing_rate, pyo.units.kg / pyo.units.s), "CaOH2"
        )

    @staticmethod
    def cost_electrodialysis(blk, cost_electricity_flow=True):
        """
        Function for costing the Electrodialysis unit

        Args:
            cost_electricity_flow - Option for including the costing of electricity
        """
        t0 = blk.flowsheet().time.first()

        membrane_cost = (
            blk.costing_package.electrodialysis.cem_membrane_cost
            + blk.costing_package.electrodialysis.aem_membrane_cost
        )
        spacer_cost = 2.0 * blk.costing_package.electrodialysis.flowspacer_cost
        electrode_cost = 2.0 * blk.costing_package.electrodialysis.electrode_cost

        cost_electrodialysis_stack(
            blk,
            membrane_cost,
            spacer_cost,
            blk.costing_package.electrodialysis.factor_membrane_housing_replacement,
            electrode_cost,
            blk.costing_package.electrodialysis.factor_electrode_replacement,
        )

        # Changed this to grab power from performance table which is identified
        # by same key regardless of whether the Electrodialysis unit is 0D or 1D
        if cost_electricity_flow:
            blk.costing_package.cost_flow(
                pyo.units.convert(
                    blk.unit_model.get_power_electrical(t0),
                    to_units=pyo.units.kW,
                ),
                "electricity",
            )

    @staticmethod
    def cost_crystallizer(blk, cost_type=CrystallizerCostType.default):
        """
        Function for costing the FC crystallizer by the mass flow of produced crystals.
        The operating cost model assumes that heat is supplied via condensation of saturated steam (see Dutta et al.)

        Args:
            cost_type - Option for crystallizer cost function type - volume or mass basis
        """
        if (
            cost_type == CrystallizerCostType.default
            or cost_type == CrystallizerCostType.mass_basis
        ):
            WaterTAPCostingData.cost_crystallizer_by_crystal_mass(blk)
        elif cost_type == CrystallizerCostType.volume_basis:
            WaterTAPCostingData.cost_crystallizer_by_volume(blk)
        else:
            raise ConfigurationError(
                f"{blk.unit_model.name} received invalid argument for cost_type:"
                f" {cost_type}. Argument must be a member of the CrystallizerCostType Enum."
            )

        blk.costing_package.cost_flow(
            pyo.units.convert(
                (
                    blk.unit_model.magma_circulation_flow_vol
                    * blk.unit_model.dens_mass_slurry
                    * Constants.acceleration_gravity
                    * blk.costing_package.crystallizer.pump_head_height
                    / blk.costing_package.crystallizer.efficiency_pump
                ),
                to_units=pyo.units.kW,
            ),
            "electricity",
        )

        blk.costing_package.cost_flow(
            pyo.units.convert(
                (
                    blk.unit_model.work_mechanical[0]
                    / WaterTAPCostingData._compute_steam_properties(blk)
                ),
                to_units=pyo.units.m**3 / pyo.units.s,
            ),
            "steam",
        )

    @staticmethod
    def cost_crystallizer_by_crystal_mass(blk):
        """
        Mass-based capital cost for FC crystallizer
        """
        make_capital_cost_var(blk)
        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost
            == pyo.units.convert(
                (
                    blk.costing_package.crystallizer.iec_percent
                    * blk.costing_package.crystallizer.fob_unit_cost
                    * (
                        sum(
                            blk.unit_model.solids.flow_mass_phase_comp[0, "Sol", j]
                            for j in blk.unit_model.config.property_package.solute_set
                        )
                        / blk.costing_package.crystallizer.ref_capacity
                    )
                    ** blk.costing_package.crystallizer.ref_exponent
                ),
                to_units=blk.costing_package.base_currency,
            )
        )

    @staticmethod
    def cost_crystallizer_by_volume(blk):
        """
        Volume-based capital cost for FC crystallizer
        """
        make_capital_cost_var(blk)
        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost
            == pyo.units.convert(
                (
                    blk.costing_package.crystallizer.volume_cost
                    * (
                        (
                            pyo.units.convert(
                                blk.unit_model.volume_suspension
                                * (
                                    blk.unit_model.height_crystallizer
                                    / blk.unit_model.height_slurry
                                ),
                                to_units=(pyo.units.ft) ** 3,
                            )
                        )
                        / pyo.units.ft**3
                    )
                    ** blk.costing_package.crystallizer.vol_basis_exponent
                ),
                to_units=blk.costing_package.base_currency,
            )
        )

    @staticmethod
    def cost_ion_exchange(blk):
        """
        Volume-based capital cost for Ion Exchange
        """
        make_capital_cost_var(blk)
        make_fixed_operating_cost_var(blk)
        # Conversions to use units from cost equations in reference
        col_vol_gal = pyo.units.convert(
            blk.unit_model.col_vol_per, to_units=pyo.units.gal
        )
        bed_vol_ft3 = pyo.units.convert(
            blk.unit_model.bed_vol, to_units=pyo.units.ft**3
        )
        bed_mass_ton = pyo.units.convert(
            blk.unit_model.bed_vol * blk.unit_model.resin_bulk_dens,
            to_units=pyo.units.ton,
        )
        bw_tank_vol = pyo.units.convert(
            (
                blk.unit_model.bw_flow * blk.unit_model.t_bw
                + blk.unit_model.rinse_flow * blk.unit_model.t_rinse
            ),
            to_units=pyo.units.gal,
        )
        regen_tank_vol = pyo.units.convert(
            (
                blk.unit_model.properties_regen[0].flow_vol_phase["Liq"]
                * blk.unit_model.t_regen
            ),
            to_units=pyo.units.gal,
        )
        ix_type = blk.unit_model.ion_exchange_type
        TIC = 1.65
        blk.capital_cost_vessel = pyo.Var(
            initialize=1e5,
            domain=pyo.NonNegativeReals,
            units=blk.costing_package.base_currency,
            doc="Capital cost for one vessel",
        )
        blk.capital_cost_resin = pyo.Var(
            initialize=1e5,
            domain=pyo.NonNegativeReals,
            units=blk.costing_package.base_currency,
            doc="Capital cost for resin for one vessel",
        )
        blk.capital_cost_regen_tank = pyo.Var(
            initialize=1e5,
            domain=pyo.NonNegativeReals,
            units=blk.costing_package.base_currency,
            doc="Capital cost for regeneration solution tank",
        )
        blk.capital_cost_backwash_tank = pyo.Var(
            initialize=1e5,
            domain=pyo.NonNegativeReals,
            units=blk.costing_package.base_currency,
            doc="Capital cost for backwash + rinse solution tank",
        )
        blk.operating_cost_hazardous = pyo.Var(
            initialize=1e5,
            domain=pyo.NonNegativeReals,
            units=blk.costing_package.base_currency,
            doc="Operating cost for hazardous waste disposal",
        )

        ion_exchange_params = blk.costing_package.ion_exchange

        # TODO: add way to have other regen chemicals
        if ix_type == "cation":
            resin_cost = ion_exchange_params.cation_exchange_resin_cost

        elif ix_type == "anion":
            resin_cost = ion_exchange_params.anion_exchange_resin_cost

        elif ix_type == "mixed":
            raise ConfigurationError(
                "IonExchangeOD model for IonExchangeType.mixed has not been implemented yet."
            )

        blk.capital_cost_vessel_constraint = pyo.Constraint(
            expr=blk.capital_cost_vessel
            == ion_exchange_params.vessel_intercept
            + ion_exchange_params.vessel_A_coeff * col_vol_gal**3
            + ion_exchange_params.vessel_B_coeff * col_vol_gal**2
            + ion_exchange_params.vessel_C_coeff * col_vol_gal
        )
        blk.capital_cost_resin_constraint = pyo.Constraint(
            expr=blk.capital_cost_resin == resin_cost * bed_vol_ft3
        )
        blk.capital_cost_backwash_tank_constraint = pyo.Constraint(
            expr=blk.capital_cost_backwash_tank
            == ion_exchange_params.backwash_tank_intercept
            + ion_exchange_params.backwash_tank_A_coeff * bw_tank_vol**3
            + ion_exchange_params.backwash_tank_B_coeff * bw_tank_vol**2
            + ion_exchange_params.backwash_tank_C_coeff * bw_tank_vol
        )
        blk.capital_cost_regen_tank_constraint = pyo.Constraint(
            expr=blk.capital_cost_regen_tank
            == ion_exchange_params.regen_tank_intercept
            + ion_exchange_params.regen_tank_A_coeff * regen_tank_vol**2
            + ion_exchange_params.regen_tank_B_coeff * regen_tank_vol
        )
        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost
            == (
                (blk.capital_cost_vessel + blk.capital_cost_resin)
                * (blk.unit_model.number_columns + blk.unit_model.number_columns_redund)
                + blk.capital_cost_backwash_tank
                + blk.capital_cost_regen_tank
            )
            * TIC
        )
        if blk.unit_model.config.hazardous_waste:
            blk.operating_cost_hazardous_constraint = pyo.Constraint(
                expr=blk.operating_cost_hazardous
                == (
                    bw_tank_vol * ion_exchange_params.hazardous_regen_disposal
                    + bed_mass_ton * ion_exchange_params.hazardous_resin_disposal
                )
                * ion_exchange_params.annual_resin_replacement_factor
                + ion_exchange_params.hazardous_min_cost
            )
        else:
            blk.operating_cost_hazardous.fix(0)
        blk.fixed_operating_cost_constraint = pyo.Constraint(
            expr=blk.fixed_operating_cost
            == (
                (
                    bed_vol_ft3
                    * ion_exchange_params.annual_resin_replacement_factor
                    * resin_cost
                )
                + blk.operating_cost_hazardous
            )
        )

        regen_soln_flow = (
            (
                blk.unit_model.regen_dose
                * blk.unit_model.bed_vol
                * (blk.unit_model.number_columns + blk.unit_model.number_columns_redund)
            )
            / (blk.unit_model.t_cycle)
        ) / blk.unit_model.regen_recycle

        electricity_flow = (
            blk.unit_model.main_pump_power
            + blk.unit_model.regen_pump_power
            + blk.unit_model.bw_pump_power
            + blk.unit_model.rinse_pump_power
        )

        blk.costing_package.cost_flow(electricity_flow, "electricity")
        blk.costing_package.cost_flow(regen_soln_flow, blk.unit_model.regen_chem)

    def cost_gac(blk):
        """
        3 equation capital cost estimation for GAC systems with: (i), contactor/pressure vessel cost by polynomial
        as a function of individual contactor volume; (ii), initial charge of GAC adsorbent cost by exponential as a
        function of required mass of GAC adsorbent; and (iii), other process costs (vessels, pipes, instrumentation, and
        controls) calculated by power law as a function of total contactor(s) volume. Operating costs calculated as the
        required makeup and regeneration of GAC adsorbent. Energy for backwash and booster pumps considered negligible
        compared to regeneration costs
        """
        make_capital_cost_var(blk)
        blk.contactor_cost = pyo.Var(
            initialize=1e5,
            domain=pyo.NonNegativeReals,
            units=blk.costing_package.base_currency,
            doc="Unit contactor(s) capital cost",
        )
        blk.bed_mass_gac_ref = pyo.Var(
            initialize=4,
            domain=pyo.NonNegativeReals,
            units=pyo.units.kg,
            doc="Reference value of GAC mass needed for initial charge where "
            "economy of scale no longer discounts the unit price",
        )
        blk.adsorbent_unit_cost = pyo.Var(
            initialize=2,
            domain=pyo.NonNegativeReals,
            units=blk.costing_package.base_currency * pyo.units.kg**-1,
            doc="GAC adsorbent cost per unit mass",
        )
        blk.adsorbent_cost = pyo.Var(
            initialize=1e5,
            domain=pyo.NonNegativeReals,
            units=blk.costing_package.base_currency,
            doc="Unit adsorbent capital cost",
        )
        blk.other_process_cost = pyo.Var(
            initialize=1e5,
            domain=pyo.NonNegativeReals,
            units=blk.costing_package.base_currency,
            doc="Unit other process capital cost",
        )

        blk.contactor_cost_constraint = pyo.Constraint(
            expr=blk.contactor_cost
            == (
                blk.costing_package.gac.num_contactors_op
                + blk.costing_package.gac.num_contactors_redundant
            )
            * pyo.units.convert(
                (
                    blk.costing_package.gac.contactor_cost_coeff_3
                    * (
                        blk.unit_model.bed_volume
                        / blk.costing_package.gac.num_contactors_op
                    )
                    ** 3
                    + blk.costing_package.gac.contactor_cost_coeff_2
                    * (
                        blk.unit_model.bed_volume
                        / blk.costing_package.gac.num_contactors_op
                    )
                    ** 2
                    + blk.costing_package.gac.contactor_cost_coeff_1
                    * (
                        blk.unit_model.bed_volume
                        / blk.costing_package.gac.num_contactors_op
                    )
                    ** 1
                    + blk.costing_package.gac.contactor_cost_coeff_0
                ),
                to_units=blk.costing_package.base_currency,
            )
        )
        blk.bed_mass_gac_ref_constraint = pyo.Constraint(
            expr=blk.bed_mass_gac_ref
            == smooth_min(
                blk.costing_package.gac.bed_mass_max_ref / pyo.units.kg,
                pyo.units.convert(blk.unit_model.bed_mass_gac, to_units=pyo.units.kg)
                / pyo.units.kg,
            )
            * pyo.units.kg
        )
        blk.adsorbent_unit_cost_constraint = pyo.Constraint(
            expr=blk.adsorbent_unit_cost
            == pyo.units.convert(
                blk.costing_package.gac.adsorbent_unit_cost_coeff
                * pyo.exp(
                    blk.bed_mass_gac_ref
                    * blk.costing_package.gac.adsorbent_unit_cost_exp_coeff
                ),
                to_units=blk.costing_package.base_currency * pyo.units.kg**-1,
            )
        )
        blk.adsorbent_cost_constraint = pyo.Constraint(
            expr=blk.adsorbent_cost
            == blk.adsorbent_unit_cost * blk.unit_model.bed_mass_gac
        )
        blk.other_process_cost_constraint = pyo.Constraint(
            expr=blk.other_process_cost
            == pyo.units.convert(
                (
                    blk.costing_package.gac.other_cost_coeff
                    * ((pyo.units.m**3) ** -blk.costing_package.gac.other_cost_exp)
                    * (
                        (
                            blk.costing_package.gac.num_contactors_op
                            + blk.costing_package.gac.num_contactors_redundant
                        )
                        * (
                            blk.unit_model.bed_volume
                            / blk.costing_package.gac.num_contactors_op
                        )
                    )
                    ** blk.costing_package.gac.other_cost_exp
                ),
                to_units=blk.costing_package.base_currency,
            )
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost
            == blk.contactor_cost + blk.adsorbent_cost + blk.other_process_cost
        )

        make_fixed_operating_cost_var(blk)
        blk.gac_regen_cost = pyo.Var(
            initialize=1e5,
            domain=pyo.NonNegativeReals,
            units=blk.costing_package.base_currency / blk.costing_package.base_period,
            doc="Cost to regenerate spent GAC adsorbent by an offsite regeneration facility",
        )
        blk.gac_makeup_cost = pyo.Var(
            initialize=1e5,
            domain=pyo.NonNegativeReals,
            units=blk.costing_package.base_currency / blk.costing_package.base_period,
            doc="Cost to makeup spent GAC adsorbent with fresh adsorbent",
        )

        blk.gac_regen_cost_constraint = pyo.Constraint(
            expr=blk.gac_regen_cost
            == pyo.units.convert(
                (
                    blk.costing_package.gac.regen_unit_cost
                    * (
                        blk.costing_package.gac.regen_frac
                        * blk.unit_model.gac_mass_replace_rate
                    )
                ),
                to_units=blk.costing_package.base_currency
                / blk.costing_package.base_period,
            )
        )
        blk.gac_makeup_cost_constraint = pyo.Constraint(
            expr=blk.gac_makeup_cost
            == pyo.units.convert(
                (
                    blk.costing_package.gac.makeup_unit_cost
                    * (
                        (1 - blk.costing_package.gac.regen_frac)
                        * blk.unit_model.gac_mass_replace_rate
                    )
                ),
                to_units=blk.costing_package.base_currency
                / blk.costing_package.base_period,
            )
        )
        blk.fixed_operating_cost_constraint = pyo.Constraint(
            expr=blk.fixed_operating_cost == blk.gac_regen_cost + blk.gac_makeup_cost
        )

    def _compute_steam_properties(blk):
        """
        Function for computing saturated steam properties for thermal heating estimation.

        Args:
            pressure_sat:   Steam gauge pressure in bar

        Out:
            Steam thermal capacity (latent heat of condensation * density) in kJ/m3
        """
        pressure_sat = blk.costing_package.crystallizer.steam_pressure
        # 1. Compute saturation temperature of steam: computed from El-Dessouky expression
        tsat_constants = [
            42.6776 * pyo.units.K,
            -3892.7 * pyo.units.K,
            1000 * pyo.units.kPa,
            -9.48654 * pyo.units.dimensionless,
        ]
        psat = (
            pyo.units.convert(pressure_sat, to_units=pyo.units.kPa)
            + 101.325 * pyo.units.kPa
        )
        temperature_sat = tsat_constants[0] + tsat_constants[1] / (
            pyo.log(psat / tsat_constants[2]) + tsat_constants[3]
        )

        # 2. Compute latent heat of condensation/vaporization: computed from Sharqawy expression
        t = temperature_sat - 273.15 * pyo.units.K
        enth_mass_units = pyo.units.J / pyo.units.kg
        t_inv_units = pyo.units.K**-1
        dh_constants = [
            2.501e6 * enth_mass_units,
            -2.369e3 * enth_mass_units * t_inv_units**1,
            2.678e-1 * enth_mass_units * t_inv_units**2,
            -8.103e-3 * enth_mass_units * t_inv_units**3,
            -2.079e-5 * enth_mass_units * t_inv_units**4,
        ]
        dh_vap = (
            dh_constants[0]
            + dh_constants[1] * t
            + dh_constants[2] * t**2
            + dh_constants[3] * t**3
            + dh_constants[4] * t**4
        )
        dh_vap = pyo.units.convert(dh_vap, to_units=pyo.units.kJ / pyo.units.kg)

        # 3. Compute specific volume: computed from Affandi expression (Eq 5)
        t_critical = 647.096 * pyo.units.K
        t_red = temperature_sat / t_critical  # Reduced temperature
        sp_vol_constants = [
            -7.75883 * pyo.units.dimensionless,
            3.23753 * pyo.units.dimensionless,
            2.05755 * pyo.units.dimensionless,
            -0.06052 * pyo.units.dimensionless,
            0.00529 * pyo.units.dimensionless,
        ]
        log_sp_vol = (
            sp_vol_constants[0]
            + sp_vol_constants[1] * (pyo.log(1 / t_red)) ** 0.4
            + sp_vol_constants[2] / (t_red**2)
            + sp_vol_constants[3] / (t_red**4)
            + sp_vol_constants[4] / (t_red**5)
        )
        sp_vol = pyo.exp(log_sp_vol) * pyo.units.m**3 / pyo.units.kg

        # 4. Return specific energy: density * latent heat
        return dh_vap / sp_vol


# Define default mapping of costing methods to unit models
WaterTAPCostingData.unit_mapping = {
    Mixer: WaterTAPCostingData.cost_mixer,
    Pump: WaterTAPCostingData.cost_pump,
    EnergyRecoveryDevice: WaterTAPCostingData.cost_energy_recovery_device,
    PressureExchanger: WaterTAPCostingData.cost_pressure_exchanger,
    ReverseOsmosis0D: WaterTAPCostingData.cost_reverse_osmosis,
    ReverseOsmosis1D: WaterTAPCostingData.cost_reverse_osmosis,
    NanoFiltration0D: WaterTAPCostingData.cost_nanofiltration,
    NanofiltrationZO: WaterTAPCostingData.cost_nanofiltration,
    Crystallization: WaterTAPCostingData.cost_crystallizer,
    Ultraviolet0D: WaterTAPCostingData.cost_uv_aop,
    Electrodialysis0D: WaterTAPCostingData.cost_electrodialysis,
    Electrodialysis1D: WaterTAPCostingData.cost_electrodialysis,
    IonExchange0D: WaterTAPCostingData.cost_ion_exchange,
    GAC: WaterTAPCostingData.cost_gac,
}


def make_capital_cost_var(blk):
    blk.capital_cost = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Unit capital cost",
    )


def make_fixed_operating_cost_var(blk):
    blk.fixed_operating_cost = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency / blk.costing_package.base_period,
        doc="Unit fixed operating cost",
    )


def cost_membrane(blk, membrane_cost, factor_membrane_replacement):
    """
    Generic function for costing a membrane. Assumes the unit_model
    has an `area` variable or parameter.

    Args:
        membrane_cost - The cost of the membrane in currency per area
        factor_membrane_replacement - Membrane replacement factor
                                      [fraction of membrane replaced/year]
    """

    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)
    blk.membrane_cost = pyo.Expression(expr=membrane_cost)
    blk.factor_membrane_replacement = pyo.Expression(expr=factor_membrane_replacement)

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost == blk.membrane_cost * blk.unit_model.area
    )
    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost
        == blk.factor_membrane_replacement * blk.membrane_cost * blk.unit_model.area
    )


def cost_electrodialysis_stack(
    blk,
    membrane_cost,
    spacer_cost,
    membrane_replacement_factor,
    electrode_cost,
    electrode_replacement_factor,
):
    """
    Generic function for costing the stack in an electrodialysis unit.
    Assumes the unit_model has a `cell_pair_num`, `cell_width`, and `cell_length`
    set of variables used to size the total membrane area.

    Args:
        membrane_cost - The total cost of the CEM and AEM per cell pair in currency per area

        spacer_cost - The total cost of the spacers per cell pair in currency per area

        membrane_replacement_factor - Replacement factor for membranes and spacers
                                      [fraction of membranes/spacers replaced/year]

        electrode_cost - The total cost of electrodes in a given stack in currency per area

        electrode_replacement_factor - Replacement factor for electrodes
                                        [fraction of electrodes replaced/year]
    """
    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    blk.membrane_cost = pyo.Expression(expr=membrane_cost)
    blk.membrane_replacement_factor = pyo.Expression(expr=membrane_replacement_factor)
    blk.spacer_cost = pyo.Expression(expr=spacer_cost)
    blk.electrode_cost = pyo.Expression(expr=electrode_cost)
    blk.electrode_replacement_factor = pyo.Expression(expr=electrode_replacement_factor)

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == (blk.membrane_cost + blk.spacer_cost)
        * (
            blk.unit_model.cell_pair_num
            * blk.unit_model.cell_width
            * blk.unit_model.cell_length
        )
        + blk.electrode_cost * (blk.unit_model.cell_width * blk.unit_model.cell_length)
    )
    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost
        == blk.membrane_replacement_factor
        * (blk.membrane_cost + blk.spacer_cost)
        * (
            blk.unit_model.cell_pair_num
            * blk.unit_model.cell_width
            * blk.unit_model.cell_length
        )
        + blk.electrode_replacement_factor
        * blk.electrode_cost
        * (blk.unit_model.cell_width * blk.unit_model.cell_length)
    )


def cost_by_flow_volume(blk, flow_cost, flow_to_cost):
    """
    Generic function for costing by flow volume.

    Args:
        flow_cost - The cost of the device in [currency]/([volume]/[time])
        flow_to_cost - The flow costed in [volume]/[time]
    """
    make_capital_cost_var(blk)
    blk.flow_cost = pyo.Expression(expr=flow_cost)
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost == blk.flow_cost * flow_to_cost
    )


def cost_uv_aop_bundle(blk, reactor_cost, lamp_cost, factor_lamp_replacement):
    """
    Generic function for costing a UV system.

    Args:
        reactor_cost - The cost of UV reactor in [currency]/[volume]
        lamp_cost - The costs of the lamps, sleeves, ballasts and sensors in [currency]/[kW]
    """
    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)
    blk.reactor_cost = pyo.Expression(expr=reactor_cost)
    blk.lamp_cost = pyo.Expression(expr=lamp_cost)
    blk.factor_lamp_replacement = pyo.Expression(expr=factor_lamp_replacement)

    flow_in = pyo.units.convert(
        blk.unit_model.control_volume.properties_in[0].flow_vol,
        to_units=pyo.units.m**3 / pyo.units.hr,
    )

    electricity_demand = pyo.units.convert(
        blk.unit_model.electricity_demand[0], to_units=pyo.units.kW
    )

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.reactor_cost * flow_in + blk.lamp_cost * electricity_demand
    )
    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost
        == blk.factor_lamp_replacement * blk.lamp_cost * electricity_demand
    )
