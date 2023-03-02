###############################################################################
# WaterTAP Copyright (c) 2021-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
###############################################################################

"""
This module contains a zero-order representation of an electrocoagulation unit.
"""

from copy import deepcopy

from pyomo.environ import (
    Var,
    Param,
    Constraint,
    Expression,
    log,
    Suffix,
    units as pyunits,
)
from pyomo.common.config import ConfigValue, In

from idaes.core import declare_process_block_class
from idaes.core.util.constants import Constants
from idaes.core.util.misc import StrEnum
import idaes.core.util.scaling as iscale
from watertap.core import build_sido, ZeroOrderBaseData


class ElectrodeMaterial(StrEnum):
    aluminum = "aluminum"
    iron = "iron"


class ReactorMaterial(StrEnum):
    pvc = "pvc"
    stainless_steel = "stainless_steel"


@declare_process_block_class("ElectrocoagulationZO")
class ElectrocoagulationZOData(ZeroOrderBaseData):
    """
    Zero-order model for an electrocoagulation unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    CONFIG.declare(
        "electrode_material",
        ConfigValue(
            default="aluminum",
            domain=In(ElectrodeMaterial),
            description="Electrode material",
        ),
    )

    CONFIG.declare(
        "reactor_material",
        ConfigValue(
            default="pvc",
            domain=In(ReactorMaterial),
            description="Reactor material",
        ),
    )

    def build(self):
        super().build()

        self._tech_type = "electrocoagulation"

        build_sido(self)

        if self.config.electrode_material == ElectrodeMaterial.aluminum:
            self.mw_electrode_material = Param(
                initialize=0.027,
                units=pyunits.kg / pyunits.mol,
                doc="Molecular weight of electrode material",
            )
            self.valence_electrode_material = Param(
                initialize=3,
                units=pyunits.dimensionless,
                doc="Number of valence electrons of electrode material",
            )
            self.density_electrode_material = Param(
                initialize=2710,
                units=pyunits.kg / pyunits.m**3,
                doc="Density of electrode material",
            )
            self.metal_dose_to_toc_ratio = Param(
                initialize=1,
                units=pyunits.kg / pyunits.kg,
                doc="Coagulant dose to inlet TOC ratio",
            )

        elif self.config.electrode_material == ElectrodeMaterial.iron:
            self.mw_electrode_material = Param(
                initialize=0.056,
                units=pyunits.kg / pyunits.mol,
                doc="Molecular weight of electrode material",
            )
            self.valence_electrode_material = Param(
                initialize=2,
                units=pyunits.dimensionless,
                doc="Number of valence electrons of electrode material",
            )
            self.density_electrode_material = Param(
                initialize=7860,
                units=pyunits.kg / pyunits.m**3,
                doc="Density of electrode material",
            )
            self.metal_dose_to_toc_ratio = Param(
                initialize=2,
                units=pyunits.kg / pyunits.kg,
                doc="Coagulant dose to inlet TOC ratio",
            )

        # Fixed parameters

        # self.current_density = Var(
        #     units=pyunits.mA / pyunits.cm**2,
        #     doc="Currenty density"
        # )

        # self.electrode_spacing = Var(
        #     units=pyunits.cm,
        #     doc="Distance between the electrodes",
        # )

        # self.solution_conductivity = Var(
        #     units=pyunits.S / pyunits.cm,
        #     doc="Electrical conductivity of solution",
        # )

        # self.power_density_k1 = Var(
        #     units=pyunits.dimensionless,
        #     doc="Constant k1 in power density equation",
        # )

        # self.power_density_k2 = Var(
        #     units=pyunits.dimensionless,
        #     doc="Constant k2 in power density equation",
        # )

        # self.z_valence = Var(
        #     units=pyunits.dimensionless,
        #     doc="Valence of dissolving metal ion for energy consumption equation"
        # )

        # self._fixed_perf_vars.append(self.current_density)
        # self._fixed_perf_vars.append(self.electrode_spacing)
        # self._fixed_perf_vars.append(self.solution_conductivity)
        # self._fixed_perf_vars.append(self.power_density_k1)
        # self._fixed_perf_vars.append(self.power_density_k2)
        # self._fixed_perf_vars.append(self.z_valence)

        # # Variable parameters and constraints

        # self.energy_consumption = Var(
        #     self.flowsheet().time,
        #     units=pyunits.kWh / pyunits.m**3 / pyunits.mM,
        #     doc="Energy required for treating a given volume",
        # )

        # @self.Constraint(self.flowsheet().time, doc="Energy requirement constraint")
        # def energy_consumption_constraint(b, t):
        #     ohmic_resistance = b.electrode_spacing / b.solution_conductivity
        #     return b.energy_consumption[t] == (
        #         b.current_density * ohmic_resistance
        #         + b.power_density_k1 * log(b.current_density)
        #         + b.power_density_k2
        #     ) * (b.z_valence * 96485 / (3600 * 10**6))

        # # Store contents for reporting output

        # self._perf_var_dict["Current density (mA/cm²)"] = self.current_density
        # self._perf_var_dict["Energy Consumption (kWh/m³/mM)"] = self.energy_consumption

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        self.current_per_reactor = Param(
            initialize=3000,
            units=pyunits.ampere,
            mutable=True,
            doc="Current required per reactor",
        )

        self.tds_to_cond_conversion = Param(
            initialize=5e3,
            mutable=True,
            units=(pyunits.mg * pyunits.m) / (pyunits.liter * pyunits.S),
            doc="Conersion factor for mg/L TDS to S/m",
        )

        self.electrode_width = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.m,
            doc="Electrode width",
        )

        self.electrode_height = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.m,
            doc="Electrode height",
        )

        self.electrode_thick = Var(
            initialize=0.001,
            bounds=(0, 0.1),
            units=pyunits.m,
            doc="Electrode thickness",
        )

        self.electrode_mass = Var(
            initialize=10,
            bounds=(0, None),
            units=pyunits.kg,
            doc="Electrode mass",
        )

        self.electrode_area_total = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.m**2,
            doc="Total electrode area",
        )

        self.electrode_area_per = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.m**2,
            doc="Electrode area",
        )

        self.electrode_volume_per = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.m**3,
            doc="Electrode volume",
        )

        self.electrode_gap = Var(
            initialize=0.005,
            bounds=(0.001, 0.2),
            units=pyunits.m,
            doc="Electrode gap",
        )

        self.electrolysis_time = Var(
            initialize=30,
            bounds=(2, 200),
            units=pyunits.minute,
            doc="Electrolysis time",
        )

        self.number_electrode_pairs = Var(
            initialize=5,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Number of electrode pairs",
        )

        self.number_cells = Var(
            initialize=5,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Number of cells",
        )

        self.applied_current = Var(
            initialize=1e4,
            bounds=(0, None),
            units=pyunits.ampere,
            doc="Applied current",
        )

        self.current_efficiency = Var(
            initialize=1,
            bounds=(0.9, 2.5),
            units=pyunits.kg / pyunits.kg,
            doc="Current efficiency",
        )

        self.cell_voltage = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.volt,
            doc="Cell voltage",
        )

        self.overpotential = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.volt,
            doc="Potential balance for voltage",  # ???
        )

        self.reactor_volume = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.m**3,
            doc="Reactor volume total (electrochemical + flotation + sedimentation)",
        )

        self.metal_loading = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.kg / pyunits.liter,
            doc="Metal loading",
        )

        self.ohmic_resistance = Var(
            initialize=1e-5,
            bounds=(0, None),
            units=pyunits.ohm,
            doc="Ohmic resistance of solution",
        )

        self.charge_loading_rate = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.coulomb / pyunits.liter,
            doc="Charge loading rate",
        )

        self.current_density = Var(
            initialize=1,
            bounds=(1, 2000),
            units=pyunits.ampere / pyunits.m**2,
            doc="Current density",
        )

        self._fixed_perf_vars.append(self.electrode_thick)
        self._fixed_perf_vars.append(self.current_density)
        self._fixed_perf_vars.append(self.electrolysis_time)
        self._fixed_perf_vars.append(self.metal_loading)
        self._fixed_perf_vars.append(self.number_electrode_pairs)
        self._fixed_perf_vars.append(self.electrode_gap)
        self._fixed_perf_vars.append(self.current_efficiency)
        self._fixed_perf_vars.append(self.overpotential)

        @self.Constraint(doc="Charge loading rate equation")
        def eq_charge_loading_rate(b):
            return b.charge_loading_rate == (
                b.applied_current
                * pyunits.convert(b.electrolysis_time, to_units=pyunits.second)
            ) / pyunits.convert(b.reactor_volume, to_units=pyunits.liter)

        @self.Constraint(doc="Metal loading rate equation")
        def eq_metal_loading_rate(b):
            return b.metal_loading == (
                b.current_efficiency * b.charge_loading_rate * b.mw_electrode_material
            ) / (b.valence_electrode_material * Constants.faraday_constant)

        @self.Constraint(doc="Total current required")
        def eq_applied_current(b):
            flow_in = pyunits.convert(
                b.properties_in[0].flow_vol, to_units=pyunits.liter / pyunits.second
            )
            return b.applied_current == (
                flow_in
                * b.metal_loading
                * b.valence_electrode_material
                * Constants.faraday_constant
            ) / (b.current_efficiency * b.mw_electrode_material)

        @self.Constraint(doc="Total electrode area required")
        def eq_electrode_area_total(b):
            return b.electrode_area_total == b.applied_current / b.current_density

        @self.Constraint(doc="Cell voltage")
        def eq_cell_voltage(b):
            return (
                b.cell_voltage
                == b.overpotential + b.applied_current * b.ohmic_resistance
            )

        @self.Constraint(doc="Area per electrode")
        def eq_electrode_area_per(b):
            return b.electrode_area_per == b.electrode_area_total / (
                b.number_electrode_pairs * 2
            )

        @self.Constraint(doc="Electrode width")
        def eq_electrode_width(b):
            return b.electrode_width == (2 * b.electrode_area_per) ** 0.5

        @self.Constraint(doc="Electrode height")
        def eq_electrode_height(b):
            return b.electrode_height == b.electrode_area_per / b.electrode_width

        @self.Constraint(doc="Electrode volume")
        def eq_electrode_volume_per(b):
            return (
                b.electrode_volume_per
                == b.electrode_width * b.electrode_height * b.electrode_thick
            )

        @self.Constraint(doc="Total reactor volume")
        def eq_reactor_volume(b):
            flow_vol = b.properties_in[0].flow_vol
            return (
                b.reactor_volume
                == pyunits.convert(
                    flow_vol * b.electrolysis_time,
                    to_units=pyunits.m**3,
                )
                / b.number_cells
            )

        @self.Expression(doc="Conductivity")
        def conductivity(b):
            tds = pyunits.convert(
                b.properties_in[0].conc_mass_comp["tds"],
                to_units=pyunits.mg / pyunits.L,
            )
            return tds / b.tds_to_cond_conversion

        @self.Constraint(doc="Ohmic resistance")
        def eq_ohmic_resistance(b):
            return b.ohmic_resistance == b.electrode_gap / (
                b.conductivity * b.electrode_area_per * b.number_cells
            )

        @self.Constraint(doc="Electrode mass")
        def eq_electrode_mass(b):
            return (
                b.electrode_mass
                == b.electrode_volume_per * b.density_electrode_material
            )

    @property
    def default_costing_method(self):
        return self.cost_electrocoagulation

    @staticmethod
    def cost_electrocoagulation(blk):
        costing = blk.config.flowsheet_costing_block
        base_currency = costing.base_currency
        ec = blk.unit_model
        flow_mgd = pyunits.convert(
            ec.properties_in[0].flow_vol, to_units=pyunits.Mgallons / pyunits.day
        )
        flow_m3_yr = pyunits.convert(
            ec.properties_in[0].flow_vol, to_units=pyunits.m**3 / pyunits.year
        )
        blk.annual_sludge_flow = pyunits.convert(
            sum(
                ec.properties_byproduct[0].flow_mass_comp[j] if j != "H2O" else 0
                for j in ec.properties_byproduct[0].params.component_list
            ),
            to_units=pyunits.kg / pyunits.year,
        )
        electrode_mat = ec.config.electrode_material

        # Add cost variable and constraint
        blk.capital_cost = Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        blk.fixed_operating_cost = Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency / pyunits.year,
            bounds=(0, None),
            doc="Fixed operating cost of unit operation",
        )
        # Get parameter dict from database
        blk.parameter_dict = (
            parameter_dict
        ) = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        (
            ec_reactor_cap_base,
            ec_reactor_cap_exp,
            ec_reactor_cap_material_coeff,
            ec_reactor_cap_safety_factor,
            ec_admin_lab_cap_base,
            ec_admin_lab_cap_exp,
            ec_power_supply_base,
            ec_admin_lab_op_base,
            ec_admin_lab_op_exp,
            sludge_handling_cost,
            ec_labor_maint_factor,
            current_per_reactor,
            number_redundant_reactors,
            electrode_material_cost,
        ) = blk.unit_model._get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "ec_reactor_cap_base",
                "ec_reactor_cap_exp",
                "ec_reactor_cap_material_coeff",
                "ec_reactor_cap_safety_factor",
                "ec_admin_lab_cap_base",
                "ec_admin_lab_cap_exp",
                "ec_power_supply_base",
                "ec_admin_lab_op_base",
                "ec_admin_lab_op_exp",
                "sludge_handling_cost",
                "ec_labor_maint_factor",
                "current_per_reactor",
                "number_redundant_reactors",
                "electrode_material_cost",
            ],
        )

        if electrode_mat == "aluminum":
            costing.defined_flows["aluminum"] = 2.23 * base_currency / pyunits.kg
            costing.register_flow_type("aluminum", 2.23 * base_currency / pyunits.kg)
        if electrode_mat == "iron":
            costing.defined_flows["iron"] = 3.41 * base_currency / pyunits.kg
            costing.register_flow_type("iron", 3.41 * base_currency / pyunits.kg)

        blk.number_chambers_system = Param(
            initialize=3,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Number total chambers for system - EC chamber > flotation chamber > sedimentation chamber. All made of same material.",
        )

        blk.number_EC_reactors = Var(
            initialize=3,
            units=pyunits.dimensionless,
            doc="Number EC cells and power supplies",
        )

        blk.capital_cost_reactor = Var(
            initialize=1e4,
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc="Cost of EC reactor",
        )

        blk.capital_cost_electrodes = Var(
            initialize=1e4,
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc="Cost of EC electrodes",
        )

        blk.capital_cost_power_supply = Var(
            initialize=1e6,
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc="Cost of EC power supply",
        )

        blk.capital_cost_admin_lab = Var(
            initialize=1e4,
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc="Cost of administration + lab + building, etc.",
        )

        blk.annual_labor_maintenance = Var(
            initialize=1e4,
            units=base_currency / pyunits.year,
            bounds=(0, None),
            doc="Annual labor + maintenance cost",
        )

        blk.annual_sludge_management = Var(
            initialize=1e4,
            units=base_currency / pyunits.year,
            bounds=(0, None),
            doc="Annual sludge management cost",
        )

        blk.annual_admin_lab = Var(
            initialize=1e4,
            units=base_currency / pyunits.year,
            bounds=(0, None),
            doc="Annual administration + lab cost",
        )

        blk.number_EC_reactors_constr = Constraint(
            expr=blk.number_EC_reactors
            == ec.applied_current / current_per_reactor + number_redundant_reactors
        )

        blk.capital_cost_reactor_constraint = Constraint(
            expr=blk.capital_cost_reactor
            == (
                (
                    ec_reactor_cap_base
                    * (
                        ec.reactor_volume
                        * blk.number_EC_reactors
                        * blk.number_chambers_system
                    )
                    ** ec_reactor_cap_exp
                )
                * ec_reactor_cap_material_coeff
            )
            * ec_reactor_cap_safety_factor
        )

        blk.capital_cost_electrodes_constraint = Constraint(
            expr=blk.capital_cost_electrodes
            == (
                ec.electrode_mass
                * ec.number_electrode_pairs
                * 2
                * blk.number_EC_reactors
            )
            * electrode_material_cost
        )

        blk.capital_cost_power_supply_constraint = Constraint(
            expr=blk.capital_cost_power_supply
            == ec_power_supply_base * blk.number_EC_reactors
        )

        blk.capital_cost_other_constraint = Constraint(
            expr=blk.capital_cost_admin_lab
            == ec_admin_lab_cap_base * flow_mgd**ec_admin_lab_cap_exp
        )

        blk.capital_cost_constraint = Constraint(
            expr=blk.capital_cost
            == blk.capital_cost_reactor
            + blk.capital_cost_electrodes
            + blk.capital_cost_power_supply
            + blk.capital_cost_admin_lab
        )

        blk.annual_labor_maintenance_constraint = Constraint(
            expr=blk.annual_labor_maintenance == flow_m3_yr * ec_labor_maint_factor
        )

        blk.annual_sludge_management_constraint = Constraint(
            expr=blk.annual_sludge_management
            == blk.annual_sludge_flow * sludge_handling_cost
        )

        blk.annual_admin_lab_constraint = Constraint(
            expr=blk.annual_admin_lab
            == ec_admin_lab_op_base * flow_mgd**ec_admin_lab_op_exp
        )

        blk.fixed_operating_cost_constraint = Constraint(
            expr=blk.fixed_operating_cost
            == blk.annual_labor_maintenance
            + blk.annual_sludge_management
            + blk.annual_admin_lab
        )

        blk.annual_electrode_replacement_mass_flow = Expression(
            expr=pyunits.convert(
                ec.metal_loading * flow_m3_yr, to_units=pyunits.kg / pyunits.year
            )
        )

        blk.electricity_flow = Expression(
            expr=pyunits.convert(
                ec.applied_current * ec.cell_voltage, to_units=pyunits.kW
            )
        )

        costing.cost_flow(
            blk.annual_electrode_replacement_mass_flow, ec.config.electrode_material
        )

        costing.cost_flow(blk.electricity_flow, "electricity")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        iscale.set_scaling_factor(self.electrode_width, 10)

        iscale.set_scaling_factor(self.electrode_height, 10)

        iscale.set_scaling_factor(self.electrode_thick, 1e3)

        iscale.set_scaling_factor(self.electrode_mass, 1)

        iscale.set_scaling_factor(self.electrode_area_total, 1e-2)

        iscale.set_scaling_factor(self.electrode_area_per, 10)

        iscale.set_scaling_factor(self.electrode_volume_per, 1)

        iscale.set_scaling_factor(self.electrode_gap, 10)

        iscale.set_scaling_factor(self.electrolysis_time, 0.1)

        iscale.set_scaling_factor(self.number_electrode_pairs, 0.1)

        iscale.set_scaling_factor(self.number_cells, 1)

        iscale.set_scaling_factor(self.applied_current, 1e-4)

        iscale.set_scaling_factor(self.current_efficiency, 1)

        iscale.set_scaling_factor(self.reactor_volume, 0.1)

        iscale.set_scaling_factor(self.metal_loading, 1e6)

        iscale.set_scaling_factor(self.ohmic_resistance, 1e5)

        iscale.set_scaling_factor(self.charge_loading_rate, 1e-2)

        iscale.set_scaling_factor(self.current_density, 1e-2)

        iscale.set_scaling_factor(self.cell_voltage, 0.1)

        iscale.set_scaling_factor(self.overpotential, 0.1)

        # transforming constraints
        sf = iscale.get_scaling_factor(self.metal_loading)
        iscale.constraint_scaling_transform(self.eq_metal_loading_rate, sf)

        sf = iscale.get_scaling_factor(self.charge_loading_rate)
        iscale.constraint_scaling_transform(self.eq_charge_loading_rate, sf)

        sf = iscale.get_scaling_factor(self.ohmic_resistance)
        iscale.constraint_scaling_transform(self.eq_ohmic_resistance, sf)
