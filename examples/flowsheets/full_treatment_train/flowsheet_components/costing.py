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

from pyomo.environ import (
    Block,
    ConcreteModel,
    Constraint,
    Expression,
    Var,
    Param,
    value,
    TransformationFactory,
    units as pyunits,
)
from idaes.core import UnitModelCostingBlock
from watertap.costing import (
    WaterTAPCosting,
    PumpType,
    EnergyRecoveryDeviceType,
    MixerType,
    ROType,
)

from watertap.examples.flowsheets.full_treatment_train.flowsheet_components import (
    feed_block,
)
from watertap.examples.flowsheets.full_treatment_train.model_components import (
    unit_separator,
    unit_0DRO,
    unit_1DRO,
    property_models,
)

from watertap.examples.flowsheets.full_treatment_train.flowsheet_components.desalination import (
    build_desalination,
    solve_desalination,
    scale_desalination,
    initialize_desalination,
    display_desalination,
)

import idaes.core.util.scaling as iscale


def build_costing(m, costing_package=WaterTAPCosting, **kwargs):
    """Add costing to a given flowsheet

    Args:
        m: model
        costing_package : FlowsheetCostingBlock
    """

    # call get_costing for each unit model
    m.fs.costing = costing_package()
    # the full_treatment_train uses a lower than default value
    # for factor_maintenance_labor_chemical
    m.fs.costing.factor_maintenance_labor_chemical.fix(0.02)
    crf = m.fs.costing.factor_capital_annualization

    # Nanofiltration
    if hasattr(m.fs, "NF"):
        if kwargs["NF_type"] == "ZO":
            m.fs.NF.costing = UnitModelCostingBlock(
                flowsheet_costing_block=m.fs.costing
            )
        elif kwargs["NF_type"] == "Sep":
            raise NotImplementedError(
                "get_costing will not be implemented for the NF separator model."
            )
    if hasattr(m.fs, "pump_NF"):
        m.fs.pump_NF.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={"pump_type": PumpType.low_pressure},
        )

    # Reverse Osmosis
    if hasattr(m.fs, "RO"):
        if kwargs["RO_type"] == "0D" or kwargs["RO_type"] == "1D":
            m.fs.RO.costing = UnitModelCostingBlock(
                flowsheet_costing_block=m.fs.costing
            )
        elif kwargs["RO_type"] == "Sep":
            raise NotImplementedError(
                "get_costing will not be implemented for the RO separator model."
            )

    # Stage 2 RO
    if hasattr(m.fs, "RO2"):
        m.fs.RO2.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={"ro_type": ROType.high_pressure},
        )

    # Pump
    if hasattr(m.fs, "pump_RO"):
        m.fs.pump_RO.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={"pump_type": PumpType.high_pressure},
        )

    # Stage 2 pump
    if hasattr(m.fs, "pump_RO2"):
        m.fs.pump_RO2.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={"pump_type": PumpType.high_pressure},
        )

    # ERD
    if hasattr(m.fs, "ERD"):
        m.fs.ERD.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={
                "energy_recovery_device_type": EnergyRecoveryDeviceType.pressure_exchanger
            },
        )

    # Pretreatment
    if hasattr(m.fs, "stoich_softening_mixer_unit"):
        m.fs.stoich_softening_mixer_unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={
                "mixer_type": MixerType.CaOH2,
                "dosing_rate": m.fs.stoich_softening_mixer_unit.dosing_rate,
            },
        )

    # Post-treatment
    if hasattr(m.fs, "ideal_naocl_mixer_unit"):
        # print('FOUND CHLORINATION UNIT')
        m.fs.ideal_naocl_mixer_unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={
                "mixer_type": MixerType.NaOCl,
                "dosing_rate": m.fs.ideal_naocl_mixer_unit.dosing_rate,
            },
        )

    if hasattr(m.fs, "mixer_permeate"):
        m.fs.mixer_permeate.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={"mixer_type": MixerType.default},
        )

    # call get_system_costing for whole flowsheet
    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(m.fs.treated_flow_vol)
    m.fs.costing.add_LCOW(m.fs.treated_flow_vol)

    if hasattr(m.fs, "stoich_softening_mixer_unit"):
        m.fs.lime_softening_unit_capex = Expression(
            expr=m.fs.stoich_softening_mixer_unit.costing.capital_cost
            / m.fs.costing.annual_water_production
            * crf
        )
    else:
        m.fs.lime_softening_unit_capex = Expression(expr=0)

    if hasattr(m.fs, "ideal_naocl_mixer_unit"):
        m.fs.chlorination_unit_capex = Expression(
            expr=m.fs.ideal_naocl_mixer_unit.costing.capital_cost
            / m.fs.costing.annual_water_production
            * crf
        )
    else:
        m.fs.chlorination_unit_capex = Expression(expr=0)

    # apply scaling to cost variables and constraints
    scale_costing(m)


def scale_costing(self):
    for b_unit in self.component_objects(Block, descend_into=True):
        if hasattr(b_unit, "costing"):
            base = b_unit.costing
            for var in base.component_objects(Var):
                if iscale.get_scaling_factor(var) is None:
                    iscale.set_scaling_factor(var, 1e-3)
            for con in base.component_data_objects(Constraint):
                iscale.constraint_scaling_transform(con, 1e-3)


def display_costing(m):
    crf = m.fs.costing.factor_capital_annualization
    if not hasattr(m.fs, "pump_RO2"):
        m.fs.pump_RO2 = Block()
        m.fs.pump_RO2.costing = Block()
        m.fs.pump_RO2.costing.fixed_operating_cost = Param(initialize=0)
    if not hasattr(m.fs, "NF"):
        m.fs.NF = Block()
        m.fs.NF.costing = Block()
        m.fs.NF.costing.fixed_operating_cost = Param(initialize=0)
    if not hasattr(m.fs, "RO2"):
        m.fs.RO2 = Block()
        m.fs.RO2.costing = Block()
        m.fs.RO2.costing.fixed_operating_cost = Param(initialize=0)

    # UNITS FOR ALL COST COMPONENTS [=] $/m3 of permeate water produced
    cost_dict = {
        "LCOW": m.fs.costing.LCOW,  # Total LCOW
        "Total CAPEX": m.fs.costing.total_investment_cost
        * crf
        / m.fs.costing.annual_water_production,  # Direct + Indirect CAPEX
        "Direct CAPEX": m.fs.costing.total_capital_cost
        * crf
        / m.fs.costing.annual_water_production,  # Direct CAPEX for all system components
        "Indirect CAPEX": (
            m.fs.costing.total_investment_cost - m.fs.costing.total_capital_cost
        )
        * crf
        / m.fs.costing.annual_water_production,  # Indirect CAPEX for miscellaneous items
        "Total OPEX": m.fs.costing.total_operating_cost
        / m.fs.costing.annual_water_production,  # Total OPEX
        "Labor & Maintenance Costs": m.fs.costing.maintenance_labor_chemical_operating_cost
        / m.fs.costing.annual_water_production,
        "Total Electricity Cost": m.fs.costing.aggregate_flow_costs["electricity"]
        / m.fs.costing.annual_water_production
        if "electricity" in m.fs.costing.aggregate_flow_costs
        else 0.0,
        "Total Membrane Replacement Cost": (
            m.fs.NF.costing.fixed_operating_cost
            + m.fs.RO.costing.fixed_operating_cost
            + m.fs.RO2.costing.fixed_operating_cost
        )
        / m.fs.costing.annual_water_production,
        "Lime softener CAPEX": m.fs.lime_softening_unit_capex,
        "Lime softener OPEX": m.fs.costing.aggregate_flow_costs["CaOH2"]
        / m.fs.costing.annual_water_production
        if "CaOH2" in m.fs.costing.aggregate_flow_costs
        else 0.0,
        "Chlorination CAPEX": m.fs.chlorination_unit_capex,
        "Chlorination OPEX": m.fs.costing.aggregate_flow_costs["NaOCl"]
        / m.fs.costing.annual_water_production
        if "NaOCl" in m.fs.costing.aggregate_flow_costs
        else 0.0,
    }

    for item, val in cost_dict.items():
        print(f"{item} = {value(val)}")

    return cost_dict
