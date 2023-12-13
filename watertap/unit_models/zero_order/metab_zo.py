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
"""
This module contains a zero-order representation of a METAB bioreactor with simple reactions
(i.e., conversion fractions for key reagents and conversion ratios for other reactive species).
"""

import pyomo.environ as pyo
from pyomo.environ import units as pyunits, Var
from idaes.core import declare_process_block_class
from watertap.core import build_sido_reactive, ZeroOrderBaseData

# Some more information about this module
__author__ = "Tim Bartholomew"


@declare_process_block_class("MetabZO")
class MetabZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a METAB bioreactor
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "metab"
        build_sido_reactive(self)

        self._gas_comp = self.config.process_subtype

        # unit variables
        self.volume = Var(
            initialize=1, bounds=(0, None), units=pyunits.m**3, doc="Reactor volume"
        )
        self.hydraulic_retention_time = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.hr,
            doc="Hydraulic residence time",
        )
        self._fixed_perf_vars.append(self.hydraulic_retention_time)

        @self.Constraint(
            doc="Constraint for reactor volume based on hydraulic residence time"
        )
        def eq_reactor_volume(b):
            return b.volume == (
                pyunits.convert(
                    b.get_inlet_flow(0), to_units=pyunits.m**3 / pyunits.hour
                )
                * b.hydraulic_retention_time
            )

        # energy consumption
        self.electricity = Var(
            self.flowsheet().time,
            initialize=1,
            bounds=(0, None),
            units=pyunits.kW,
            doc="Electricity demand of unit",
        )
        self.heat = Var(
            self.flowsheet().time,
            initialize=1,
            bounds=(0, None),
            units=pyunits.kW,
            doc="Thermal demand of unit",
        )
        self.energy_electric_mixer_vol = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.kW / pyunits.m**3,
            doc="Electricity intensity of mixer with respect to reactor volume",
        )
        self._fixed_perf_vars.append(self.energy_electric_mixer_vol)
        self.energy_electric_vacuum_flow_vol_byproduct = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.kW / (pyunits.kg / pyunits.hr),
            doc="Electricity intensity of vacuum pump with respect to product gas flow",
        )
        self._fixed_perf_vars.append(self.energy_electric_vacuum_flow_vol_byproduct)
        self.energy_thermal_flow_vol_inlet = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.kJ / pyunits.m**3,
            doc="Thermal energy intensity of reactor with respect to inlet volumetric flowrate",
        )
        self._fixed_perf_vars.append(self.energy_thermal_flow_vol_inlet)

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for electricity consumption based on " "feed flowrate.",
        )
        def electricity_consumption(b, t):
            return b.electricity[t] == (
                b.energy_electric_mixer_vol * b.volume
                + b.energy_electric_vacuum_flow_vol_byproduct
                * pyunits.convert(
                    b.properties_byproduct[t].flow_mass_comp[b._gas_comp],
                    to_units=pyunits.kg / pyunits.hr,
                )
            )

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for heat demand based on " "feed flowrate.",
        )
        def heat_demand(b, t):
            return b.heat[t] == (
                b.energy_thermal_flow_vol_inlet
                * pyunits.convert(
                    b.get_inlet_flow(t), to_units=pyunits.m**3 / pyunits.s
                )
            )

        self._perf_var_dict["Electricity Demand"] = self.electricity
        self._perf_var_dict["Thermal Energy Demand"] = self.heat

    @property
    def default_costing_method(self):
        return self.cost_metab

    @staticmethod
    def cost_metab(blk):
        """
        General method for costing the metab reactor. Capital cost
        is based on the cost of reactor, mixer, METAB beads, membrane,
        and vacuum pump.
        This method also registers the electricity demand as a costed flow.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        (
            reactor_cost,
            mixer_cost,
            bead_bulk_density,
            bead_cost,
            bead_replacement_factor,
            membrane_sidestream_fraction,
            membrane_specific_size,
            membrane_cost,
            vacuum_cost,
        ) = blk.unit_model._get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "reactor_cost",
                "mixer_cost",
                "bead_bulk_density",
                "bead_cost",
                "bead_replacement_factor",
                "membrane_sidestream_fraction",
                "membrane_specific_size",
                "membrane_cost",
                "vacuum_cost",
            ],
        )

        # Add capital cost variables and constraints
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )
        blk.DCC_reactor = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Direct capital cost of reactor",
        )
        blk.DCC_mixer = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Direct capital cost of mixer",
        )
        blk.DCC_bead = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Direct capital cost of beads",
        )
        blk.DCC_membrane = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Direct capital cost of membrane",
        )
        blk.DCC_vacuum = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Direct capital cost of vacuum pump",
        )
        blk.eq_DCC_reactor = pyo.Constraint(
            expr=blk.DCC_reactor
            == pyo.units.convert(
                blk.unit_model.volume * reactor_cost,
                to_units=blk.config.flowsheet_costing_block.base_currency,
            )
        )
        blk.eq_DCC_mixer = pyo.Constraint(
            expr=blk.DCC_mixer
            == pyo.units.convert(
                blk.unit_model.energy_electric_mixer_vol
                * blk.unit_model.volume
                * mixer_cost,
                to_units=blk.config.flowsheet_costing_block.base_currency,
            )
        )
        blk.eq_DCC_bead = pyo.Constraint(
            expr=blk.DCC_bead
            == pyo.units.convert(
                blk.unit_model.volume * bead_bulk_density * bead_cost,
                to_units=blk.config.flowsheet_costing_block.base_currency,
            )
        )
        blk.eq_DCC_membrane = pyo.Constraint(
            expr=blk.DCC_membrane
            == pyo.units.convert(
                blk.unit_model.get_inlet_flow(t0)
                * membrane_sidestream_fraction
                * membrane_specific_size
                * membrane_cost,
                to_units=blk.config.flowsheet_costing_block.base_currency,
            )
        )
        blk.eq_DCC_vacuum = pyo.Constraint(
            expr=blk.DCC_vacuum
            == pyo.units.convert(
                blk.unit_model.properties_byproduct[t0].flow_mass_comp[
                    blk.unit_model._gas_comp
                ]
                * vacuum_cost,
                to_units=blk.config.flowsheet_costing_block.base_currency,
            )
        )

        expr = (
            blk.DCC_reactor
            + blk.DCC_mixer
            + blk.DCC_bead
            + blk.DCC_membrane
            + blk.DCC_vacuum
        )

        blk.costing_package.add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Add fixed operating cost variable and constraint
        blk.fixed_operating_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency
            / blk.config.flowsheet_costing_block.base_period,
            bounds=(0, None),
            doc="Fixed operating cost of unit",
        )
        blk.fixed_operating_cost_constraint = pyo.Constraint(
            expr=blk.fixed_operating_cost
            == pyo.units.convert(
                blk.DCC_bead * bead_replacement_factor,
                to_units=blk.config.flowsheet_costing_block.base_currency
                / blk.config.flowsheet_costing_block.base_period,
            )
        )

        # Register operating cost flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )
        blk.config.flowsheet_costing_block.cost_flow(blk.unit_model.heat[t0], "heat")
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.properties_byproduct[t0].flow_mass_comp[
                blk.unit_model._gas_comp
            ],
            blk.unit_model._gas_comp + "_product",
        )
