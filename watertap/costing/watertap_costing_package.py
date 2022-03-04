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
from pyomo.environ import (
    Block, Constraint, Expression, Var, Param, Reals, NonNegativeReals, units as pyunits)
from idaes.core.util.exceptions import BurntToast, ConfigurationError
from idaes.core.util.misc import register_units_of_measurement
from idaes.generic_models.costing import CostingPackageBase

from watertap.unit_models import (
        ReverseOsmosis0D,
        ReverseOsmosis1D,
        NanoFiltration0D,
        NanoFiltrationZO,
        PressureExchanger,
        Pump,
        )

class ROType(str, Enum):
    standard = "standard"
    high_pressure = "high_pressure"
    
    def __str__(self):
        return self.value

class PumpType(str, Enum):
    low_pressure = "low_pressure"
    high_pressure = "high_pressure"
    pressure_exchangers = "pressure_exchanger"
    
    def __str__(self):
        return self.value

class WaterTAPCostingPackage(CostingPackageBase):

    # TODO: this is in SSLW -- we need an equivalent
    # Register currency and conversion rates based on CE Index
    register_units_of_measurement("USD_500", "[currency]")  # base USD @ CEI 500
    register_units_of_measurement("USD2010", "500/550.8 * USD_500")
    register_units_of_measurement("USD2011", "500/585.7 * USD_500")
    register_units_of_measurement("USD2012", "500/584.6 * USD_500")
    register_units_of_measurement("USD2013", "500/567.3 * USD_500")
    register_units_of_measurement("USD2014", "500/576.1* USD_500")
    register_units_of_measurement("USD2015", "500/556.8 * USD_500")
    register_units_of_measurement("USD2016", "500/541.7 * USD_500")
    register_units_of_measurement("USD2017", "500/567.5 * USD_500")
    register_units_of_measurement("USD2018", "500/671.1 * USD_500")
    register_units_of_measurement("USD2019", "500/680.0 * USD_500")
    register_units_of_measurement("USD_394", "500/394 * USD_500")  # required for pump costing

    # Set the base year for all costs
    base_currency = pyunits.USD2018
    # Set a base period for all operating costs
    base_period = pyunits.year

    # Define standard material flows and costs
    defined_flows = {
        "electricity": 0.07 * pyunits.USD2018 / pyunits.kWh}

    @staticmethod
    def build_global_params(blk):
        # Build flowsheet level costing components
        # This is package specific
        blk.load_factor = Param(
                mutable=True,
                initialize=0.9,
                doc='Load factor [fraction of uptime]',
                units=pyunits.dimensionless)
        blk.factor_total_investment = Param(
                mutable=True,
                initialize=2,
                doc='Total investment factor [investment cost/equipment cost]',
                units=pyunits.dimensionless)
        blk.factor_maintenance_labor_chemical = Param(
                mutable=True,
                initialize=0.03,
                doc='Maintenance-labor-chemical factor [fraction of investment cost/year]',
                units=pyunits.year**-1)
        blk.factor_capital_annualization = Param(
                mutable=True,
                initialize=0.1,
                doc='Capital annualization factor [fraction of investment cost/year]',
                units=pyunits.year**-1)
        blk.factor_membrane_replacement = Param(
                mutable=True,
                initialize=0.2,
                doc='Membrane replacement factor [fraction of membrane replaced/year]',
                units=pyunits.year**-1)
        blk.reverse_osmosis_membrane_cost = Param(
                mutable=True,
                initialize=30,
                doc='Membrane cost [$/m2]',
                units=pyunits.USD_500/(pyunits.meter**2))
        blk.reverse_osmosis_high_pressure_membrane_cost = Param(
                mutable=True,
                initialize=75,
                doc='Membrane cost [$/m2]',
                units=pyunits.USD_500/(pyunits.meter**2))
        blk.nanofiltration_membrane_cost = Param(
                mutable=True,
                initialize=15,
                doc='Membrane cost [$/m2]',
                units=pyunits.USD_500/(pyunits.meter**2))
        blk.high_pressure_pump_cost = Param(
                mutable=True,
                initialize=53 / 1e5 * 3600,
                doc='High pressure pump cost [$/W]',
                units=pyunits.USD_500/pyunits.watt)
        blk.low_pressure_pump_cost = Param(
                mutable=True,
                initialize=889,
                doc='Low pressure pump cost [$/(Liter/second)]',
                units=pyunits.USD_500/(pyunits.liter/pyunits.second))
        blk.pump_pressure_exchanger_cost = Param(
                mutable=True,
                initialize=535,
                doc='Pressure exchanger cost [$/(m3/h)]',
                units=pyunits.USD_500/(pyunits.meter**3/pyunits.hours))
        blk.pressure_exchanger_cost = Param(
                mutable=True,
                initialize=535,
                doc='Pressure exchanger cost [$/(m3/h)]',
                units=pyunits.USD_500/(pyunits.meter**3/pyunits.hours))

    @staticmethod
    def build_process_costs(blk):
        # Build flowsheet level total costs
        blk.total_capital_cost = Expression(expr = blk.factor_total_investment * blk.aggregate_capital_cost )
        blk.total_fixed_operating_cost = Expression(expr = \
                 blk.aggregate_fixed_operating_cost + blk.factor_maintenance_labor_chemical * blk.total_capital_cost)

        blk.total_variable_operating_cost = Expression(expr = blk.aggregate_variable_operating_cost)
        blk.total_flow_costs = Expression(expr = blk.aggregate_flow_costs)

    @staticmethod
    def cost_nanofiltration(blk):
        """
        Nanofiltration costing method

        TODO: describe equations
        """
        _make_captial_cost_var(blk)
        _make_fixed_operating_cost_var(blk)

        fcb = blk.config.flowsheet_costing_block
        _cost_membrane(blk, fcb.nanofiltration_membrane_cost, fcb.factor_membrane_replacement)

    # Define costing methods supported by package
    @staticmethod
    def cost_reverse_osmosis(blk, ro_type=ROType.standard):
        """
        Reverse osmosis costing method

        TODO: describe equations

        Args:
            ro_type - ROType Enum indicating reverse osmosis type,
                      default = ROType.standard
        """
        # Validate arguments
        if ro_type not in ROType:
            raise ConfigurationError(
                f"{blk.unit_model.name} received invalid argument for ro_type:"
                f" {ro_type}. Argument must be a member of the ROType Enum.")

        _make_captial_cost_var(blk)
        _make_fixed_operating_cost_var(blk)

        fcb = blk.config.flowsheet_costing_block
        if ro_type == ROType.standard:
            membrane_cost = fcb.reverse_osmosis_membrane_cost
        elif ro_type == ROType.high_pressure:
            membrane_cost = fac.reverse_osmosis_high_pressure_membrane_cost
        else:
            raise BurntToast(f"Unrecognized ro_type: {ro_type}")
        _cost_membrane(blk, membrane_cost, fcb.factor_membrane_replacement)

    @staticmethod
    def cost_pump(blk, pump_type=PumpType.high_pressure):
        """
        Pump costing method

        TODO: describe equations

        Args:
            pump_type - PumpType Enum indicating pump type,
                        default = PumpType.high_pressure
        """
        if pump_type not in PumpType:
            raise ConfigurationError(
                f"{blk.unit_model.name} received invalid argument for pump_type:"
                f" {pump_type}. Argument must be a member of the PumpType Enum.")

        if pump_type=PumpType.high_pressure:
            WaterTAPCostingPackage.cost_high_pressure_pump(blk)
        elif pump_type=PumpType.low_pressure:
            WaterTAPCostingPackage.cost_low_pressure_pump(blk)
        elif pump_type=PumpType.pressure_exchanger:
            WaterTAPCostingPackage.cost_pressure_exchanger_pump(blk)
        else:
            raise BurntToast(f"Unrecognized pump_type: {pump_type}")

    @staticmethod
    def cost_high_pressure_pump(blk):
        """
        High pressure pump costing method

        TODO: describe equations
        """
        _make_captial_cost_var(blk)

        fcb = blk.config.flowsheet_costing_block
        blk.eq_capital_cost = Constraint(expr = \
                capital_cost == fcb.high_pressure_pump_cost * blk.unit_model.work_mechanical[0])

    @staticmethod
    def cost_low_pressure_pump(blk):
        """
        High pressure pump costing method

        TODO: describe equations
        """
        _make_captial_cost_var(blk)
        fcb = blk.config.flowsheet_costing_block
        cost_by_flow_volume(blk, pyunits.conver(fcb.low_pressure_pump_cost, pyunits.USD_500/(pyunits.m**3/pyunits.s)),
                blk.unit_model.control_volume.properties_in[0].flow_vol)

    @staticmethod
    def cost_pressure_exchanger_pump(blk):
        """
        Pump pressure exchanger costing method

        TODO: describe equations
        """
        _make_captial_cost_var(blk)
        fcb = blk.config.flowsheet_costing_block
        cost_by_flow_volume(blk, fcb.pump_pressure_exchanger_cost, blk.unit_model.control_volume.properties_in[0].flow_vol)

    @staticmethod
    def cost_pressure_exchanger(blk):
        """
        Pressure exchanger costing method

        TODO: describe equations
        """
        _make_captial_cost_var(blk)
        fcb = blk.config.flowsheet_costing_block
        cost_by_flow_volume(blk, fcb.pressure_exchanger_cost, blk.unit_model.low_pressure_side.properties_in[0].flow_vol)


    ## TODO; Mixer and Separator

    # Define default mapping of costing methods to unit models
    unit_mapping = {
        Pump: cost_pump,  # TODO: this may be carrying over a limitation of IDAES
        PressureExchanger: cost_pressure_exchanger,
        ReverseOsmosis0D: cost_reverse_osmosis,
        ReverseOsmosis1D: cost_reverse_osmosis,
        NanoFiltration0D: cost_nanofiltration,
        NanoFiltrationZO: cost_nanofiltration,
        }

def _make_captial_cost_var(blk):
    blk.capital_cost = Var(initialize=1e5,
                           domain=NonNegativeReals,
                           units=pyunits.USD2018,
                           doc="Unit capital cost")

def _make_fixed_operating_cost_var(blk):
    blk.fixed_operating_cost = Var(initialize=1e5,
                                   domain=NonNegativeReals,
                                   units=pyunits.USD2018/pyunits.year,
                                   doc="Unit fixed operating cost")

def cost_membrane(blk, membrane_cost, factor_membrane_replacement):
    """
    Generic function for costing a membrane. Assumes the unit_model
    has an `area` variable or parameter.

    Args:
        membrane_cost - The cost of the membrane in currency per area
        factor_membrane_replacement - Membrane replacement factor
                                      [fraction of membrane replaced/year]

    """
    blk.membrane_cost = Expression(expr=membrane_cost)
    blk.factor_membrane_replacement = Expression(expr=factor_membrane_replacement)

    blk.eq_capital_cost = Constraint(expr = \
            blk.capital_cost == blk.membrane_cost * blk.unit_model.area)
    blk.eq_fixed_operating_cost = Constraint(expr = \
            blk.fixed_operating_cost == blk.factor_membrane_replacement * blk.membrane_cost * blk.unit_model.area)

def cost_pump_by_flow_volume(blk, flow_cost, flow_to_cost):
    """
    Generic function for costing a pump by flow volume.

    Args:
        flow_cost - The cost of the pump in [currency]/([volume]/[time])
        flow_to_cost - The flow costed in [volume]/[time]
    """
    blk.flow_cost = Expression(expr=flow_cost)
    blk.eq_capital_cost = Constraint(expr = \
            capital_cost = blk.flow_cost * flow_to_cost)
