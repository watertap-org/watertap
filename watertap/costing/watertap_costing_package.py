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

from idaes.models.unit_models import Mixer

from watertap.unit_models import (
    ReverseOsmosis0D,
    ReverseOsmosis1D,
    NanoFiltration0D,
    NanofiltrationZO,
    PressureExchanger,
    Crystallization,
    Pump,
    EnergyRecoveryDevice,
    Electrodialysis0D,
    Electrodialysis1D,
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

        # Build flowsheet level costing components
        # This is package specific
        self.load_factor = pyo.Var(
            initialize=0.9,
            doc="Load factor [fraction of uptime]",
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
        self.factor_membrane_replacement = pyo.Var(
            initialize=0.2,
            doc="Membrane replacement factor [fraction of membrane replaced/year]",
            units=pyo.units.year**-1,
        )
        self.reverse_osmosis_membrane_cost = pyo.Var(
            initialize=30,
            doc="Membrane cost",
            units=self.base_currency / (pyo.units.meter**2),
        )
        self.reverse_osmosis_high_pressure_membrane_cost = pyo.Var(
            initialize=75,
            doc="Membrane cost",
            units=self.base_currency / (pyo.units.meter**2),
        )
        self.nanofiltration_membrane_cost = pyo.Var(
            initialize=15,
            doc="Membrane cost",
            units=self.base_currency / (pyo.units.meter**2),
        )
        self.high_pressure_pump_cost = pyo.Var(
            initialize=53 / 1e5 * 3600,
            doc="High pressure pump cost",
            units=self.base_currency / pyo.units.watt,
        )
        self.low_pressure_pump_cost = pyo.Var(
            initialize=889,
            doc="Low pressure pump cost",
            units=self.base_currency / (pyo.units.liter / pyo.units.second),
        )
        self.erd_pressure_exchanger_cost = pyo.Var(
            initialize=535,
            doc="Pressure exchanger cost",
            units=self.base_currency / (pyo.units.meter**3 / pyo.units.hours),
        )
        self.pressure_exchanger_cost = pyo.Var(
            initialize=535,
            doc="Pressure exchanger cost",
            units=self.base_currency / (pyo.units.meter**3 / pyo.units.hours),
        )
        self.energy_recovery_device_linear_coefficient = pyo.Var(
            initialize=3134.7,
            doc="Energy recovery device linear coefficient",
            units=self.base_currency,
        )
        self.energy_recovery_device_exponent = pyo.Var(
            initialize=0.58,
            doc="Energy recovery device exponent",
            units=pyo.units.dimensionless,
        )
        self.mixer_unit_cost = pyo.Var(
            initialize=361,
            doc="Mixer cost",
            units=self.base_currency / (pyo.units.liters / pyo.units.second),
        )
        self.naocl_mixer_unit_cost = pyo.Var(
            initialize=5.08,
            doc="NaOCl mixer cost",
            units=self.base_currency / (pyo.units.m**3 / pyo.units.day),
        )
        self.caoh2_mixer_unit_cost = pyo.Var(
            initialize=792.8 * 2.20462,
            doc="Ca(OH)2 mixer cost",
            units=self.base_currency / (pyo.units.kg / pyo.units.day),
        )

        self.electricity_base_cost = pyo.Param(
            mutable=True,
            initialize=0.07,
            doc="Electricity cost",
            units=self.base_currency / pyo.units.kWh,
        )
        self.naocl_cost = pyo.Param(
            initialize=0.23, doc="NaOCl cost", units=self.base_currency / pyo.units.kg
        )
        self.naocl_purity = pyo.Param(
            mutable=True,
            initialize=0.15,
            doc="NaOCl purity",
            units=pyo.units.dimensionless,
        )
        self.caoh2_cost = pyo.Param(
            mutable=True,
            initialize=0.12,
            doc="CaOH2 cost",
            units=self.base_currency / pyo.units.kg,
        )
        self.caoh2_purity = pyo.Param(
            mutable=True,
            initialize=1,
            doc="CaOH2 purity",
            units=pyo.units.dimensionless,
        )

        self.electrodialysis_cem_membrane_cost = pyo.Var(
            initialize=43,
            doc="Cost of CEM membrane used in Electrodialysis ($/CEM/area)",
            units=self.base_currency / (pyo.units.meter**2),
        )
        self.electrodialysis_aem_membrane_cost = pyo.Var(
            initialize=43,
            doc="Cost of AEM membrane used in Electrodialysis ($/AEM/area)",
            units=self.base_currency / (pyo.units.meter**2),
        )
        self.electrodialysis_flowspacer_cost = pyo.Var(
            initialize=3,
            doc="Cost of the spacers used in Electrodialysis ($/spacer/area)",
            units=self.base_currency / (pyo.units.meter**2),
        )
        self.factor_electrodialysis_membrane_housing_replacement = pyo.Var(
            initialize=0.2,
            doc="Membrane housing replacement factor for CEM, AEM, and spacer replacements [fraction of membrane replaced/year]",
            units=pyo.units.year**-1,
        )
        self.electrodialysis_electrode_cost = pyo.Var(
            initialize=2000,
            doc="Cost of the electrodes used in Electrodialysis ($/electrode/area)",
            units=self.base_currency / (pyo.units.meter**2),
        )
        self.factor_electrodialysis_electrode_replacement = pyo.Var(
            initialize=0.02,
            doc="Electrode replacements [fraction of electrode replaced/year]",
            units=pyo.units.year**-1,
        )

        self.fc_crystallizer_fob_unit_cost = pyo.Var(
            initialize=675000,
            doc="Forced circulation crystallizer reference free-on-board cost (Woods, 2007)",
            units=pyo.units.USD_2007,
        )

        self.fc_crystallizer_ref_capacity = pyo.Var(
            initialize=1,
            doc="Forced circulation crystallizer reference crystal capacity (Woods, 2007)",
            units=pyo.units.kg / pyo.units.s,
        )

        self.fc_crystallizer_ref_exponent = pyo.Var(
            initialize=0.53,
            doc="Forced circulation crystallizer cost exponent factor (Woods, 2007)",
            units=pyo.units.dimensionless,
        )

        self.fc_crystallizer_iec_percent = pyo.Var(
            initialize=1.43,
            doc="Forced circulation crystallizer installed equipment cost (Diab and Gerogiorgis, 2017)",
            units=pyo.units.dimensionless,
        )

        self.fc_crystallizer_volume_cost = pyo.Var(
            initialize=16320,
            doc="Forced circulation crystallizer cost per volume (Yusuf et al., 2019)",
            units=pyo.units.USD_2007,  ## TODO: Needs confirmation, but data is from Perry apparently
        )

        self.fc_crystallizer_vol_basis_exponent = pyo.Var(
            initialize=0.47,
            doc="Forced circulation crystallizer volume-based cost exponent (Yusuf et al., 2019)",
            units=pyo.units.dimensionless,
        )

        # Crystallizer operating cost information from literature
        self.steam_unit_cost = pyo.Var(
            initialize=0.004,
            units=pyo.units.USD_2018 / (pyo.units.meter**3),
            doc="Steam cost, Panagopoulos (2019)",
        )

        self.crystallizer_steam_pressure = pyo.Var(
            initialize=3,
            units=pyo.units.bar,
            doc="Steam pressure (gauge) for crystallizer heating: 3 bar default based on Dutta example",
        )

        self.crystallizer_efficiency_pump = pyo.Var(
            initialize=0.7,
            units=pyo.units.dimensionless,
            doc="Crystallizer pump efficiency - assumed",
        )

        self.crystallizer_pump_head_height = pyo.Var(
            initialize=1,
            units=pyo.units.m,
            doc="Crystallizer pump head height -  assumed, unvalidated",
        )

        # fix the parameters
        for var in self.component_objects(pyo.Var):
            var.fix()

        # Define standard material flows and costs
        self.defined_flows["electricity"] = self.electricity_base_cost
        self.defined_flows["NaOCl"] = self.naocl_cost / self.naocl_purity
        self.defined_flows["CaOH2"] = self.caoh2_cost / self.caoh2_purity
        self.defined_flows["steam"] = self.steam_unit_cost

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
            + sum(self.aggregate_flow_costs.values()) * self.load_factor
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
                * self.load_factor
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
                    * self.load_factor
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

    # Define costing methods supported by package
    @staticmethod
    def cost_nanofiltration(blk):
        """
        Nanofiltration costing method

        TODO: describe equations
        """
        cost_membrane(
            blk,
            blk.costing_package.nanofiltration_membrane_cost,
            blk.costing_package.factor_membrane_replacement,
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
            membrane_cost = blk.costing_package.reverse_osmosis_membrane_cost
        elif ro_type == ROType.high_pressure:
            membrane_cost = (
                blk.costing_package.reverse_osmosis_high_pressure_membrane_cost
            )
        else:
            raise ConfigurationError(
                f"{blk.unit_model.name} received invalid argument for ro_type:"
                f" {ro_type}. Argument must be a member of the ROType Enum."
            )
        cost_membrane(
            blk, membrane_cost, blk.costing_package.factor_membrane_replacement
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
            == blk.costing_package.high_pressure_pump_cost
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
            blk.costing_package.low_pressure_pump_cost,
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
            blk.costing_package.erd_pressure_exchanger_cost,
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
            blk.costing_package.pressure_exchanger_cost,
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
            blk.costing_package.mixer_unit_cost,
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
            blk.costing_package.naocl_mixer_unit_cost,
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
            blk.costing_package.caoh2_mixer_unit_cost
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
            blk.costing_package.electrodialysis_cem_membrane_cost
            + blk.costing_package.electrodialysis_aem_membrane_cost
        )
        spacer_cost = 2.0 * blk.costing_package.electrodialysis_flowspacer_cost
        electrode_cost = 2.0 * blk.costing_package.electrodialysis_electrode_cost

        cost_electrodialysis_stack(
            blk,
            membrane_cost,
            spacer_cost,
            blk.costing_package.factor_electrodialysis_membrane_housing_replacement,
            electrode_cost,
            blk.costing_package.factor_electrodialysis_electrode_replacement,
        )

        # Changed this to grab power from performance table which is identified
        #   by same key regardless of whether the Electrodialysis unit is 0D or 1D
        if cost_electricity_flow:
            blk.costing_package.cost_flow(
                pyo.units.convert(
                    blk.unit_model._get_performance_contents(t0)["vars"][
                        "Total electrical power consumption(Watt)"
                    ],
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
                    * blk.costing_package.crystallizer_pump_head_height
                    / blk.costing_package.crystallizer_efficiency_pump
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
                    blk.costing_package.fc_crystallizer_iec_percent
                    * blk.costing_package.fc_crystallizer_fob_unit_cost
                    * (
                        sum(
                            blk.unit_model.solids.flow_mass_phase_comp[0, "Sol", j]
                            for j in blk.unit_model.config.property_package.solute_set
                        )
                        / blk.costing_package.fc_crystallizer_ref_capacity
                    )
                    ** blk.costing_package.fc_crystallizer_ref_exponent
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
                    blk.costing_package.fc_crystallizer_volume_cost
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
                    ** blk.costing_package.fc_crystallizer_vol_basis_exponent
                ),
                to_units=blk.costing_package.base_currency,
            )
        )

    def _compute_steam_properties(blk):
        """
        Function for computing saturated steam properties for thermal heating estimation.

        Args:
            pressure_sat:   Steam gauge pressure in bar

        Out:
            Steam thermal capacity (latent heat of condensation * density) in kJ/m3
        """
        pressure_sat = blk.costing_package.crystallizer_steam_pressure
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
    Electrodialysis0D: WaterTAPCostingData.cost_electrodialysis,
    Electrodialysis1D: WaterTAPCostingData.cost_electrodialysis,
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
        flow_cost - The cost of the pump in [currency]/([volume]/[time])
        flow_to_cost - The flow costed in [volume]/[time]
    """
    make_capital_cost_var(blk)
    blk.flow_cost = pyo.Expression(expr=flow_cost)
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost == blk.flow_cost * flow_to_cost
    )
