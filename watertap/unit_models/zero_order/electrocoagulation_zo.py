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

from pyomo.environ import (
    Var,
    Param,
    PositiveReals,
    Constraint,
    Expression,
    units as pyunits,
)
from pyomo.common.config import ConfigValue, In

from idaes.core import declare_process_block_class
from idaes.core.util.constants import Constants
from idaes.core.util.misc import StrEnum
from watertap.core import build_sido, ZeroOrderBaseData


class ElectrodeMaterial(StrEnum):
    """
    Electrode material can be aluminum or iron. Default is aluminum. To run
    the model with iron electrodes, use :code:`ElectrocoagulationZO(electrode_material="iron")`;
    otherwise :code:`ElectrocoagulationZO()` will default to aluminum electrodes.
    """

    aluminum = "aluminum"
    iron = "iron"


class ReactorMaterial(StrEnum):
    """
    Reactor material can be PVC or stainless steel. Default is PVC. To run the
    model with stainless steel, use :code:`ElectrocoagulationZO(reactor_material="stainless_steel")`;
    otherwise, :code:`ElectrocoagulationZO()` will default to PVC.
    """

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
                within=PositiveReals,
                doc="Molecular weight of electrode material",
            )
            self.valence_electrode_material = Param(
                initialize=3,
                units=pyunits.dimensionless,
                within=PositiveReals,
                doc="Number of valence electrons of electrode material",
            )
            self.density_electrode_material = Param(
                initialize=2710,
                units=pyunits.kg / pyunits.m**3,
                within=PositiveReals,
                doc="Density of electrode material",
            )

        elif self.config.electrode_material == ElectrodeMaterial.iron:
            self.mw_electrode_material = Param(
                initialize=0.056,
                units=pyunits.kg / pyunits.mol,
                within=PositiveReals,
                doc="Molecular weight of electrode material",
            )
            self.valence_electrode_material = Param(
                initialize=2,
                units=pyunits.dimensionless,
                within=PositiveReals,
                doc="Number of valence electrons of electrode material",
            )
            self.density_electrode_material = Param(
                initialize=7860,
                units=pyunits.kg / pyunits.m**3,
                within=PositiveReals,
                doc="Density of electrode material",
            )

        self.current_per_reactor = Param(
            initialize=3000,
            units=pyunits.ampere,
            mutable=True,
            within=PositiveReals,
            doc="Current required per reactor",
        )

        self.tds_to_cond_conversion = Param(
            initialize=5e3,
            mutable=True,
            units=(pyunits.mg * pyunits.m) / (pyunits.liter * pyunits.S),
            within=PositiveReals,
            doc="Conversion factor for mg/L TDS to S/m",
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
            initialize=1,
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
            doc="Overpotential",
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

        self.power_required = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.watt,
            doc="Power required",
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

        @self.Constraint(doc="Metal loading equation")
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

        @self.Constraint(doc="Power required")
        def eq_power_required(b):
            return b.power_required == b.cell_voltage * b.applied_current

    @property
    def default_costing_method(self):
        return self.cost_electrocoagulation

    @staticmethod
    def cost_electrocoagulation(blk):
        """
        General method for costing electrocoagulation.
        """

        ec = blk.unit_model
        costing = blk.config.flowsheet_costing_block
        base_currency = costing.base_currency

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
        reactor_mat = ec.config.reactor_material

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
            ec_power_supply_base_slope,
            ec_admin_lab_op_base,
            ec_admin_lab_op_exp,
            sludge_handling_cost,
            ec_labor_maint_factor,
            current_per_reactor,
            number_redundant_reactors,
            electrode_material_cost,
            electrode_material_cost_coeff,
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
                "ec_power_supply_base_slope",
                "ec_admin_lab_op_base",
                "ec_admin_lab_op_exp",
                "sludge_handling_cost",
                "ec_labor_maint_factor",
                "current_per_reactor",
                "number_redundant_reactors",
                "electrode_material_cost",
                "electrode_material_cost_coeff",
            ],
        )

        costing_ec = costing.electrocoagulation

        if electrode_mat == "aluminum":
            # Reference for Al cost: Anuf et al., 2022 - https://doi.org/https://doi.org/10.1016/j.jwpe.2022.103074
            costing.register_flow_type("aluminum", 2.23 * base_currency / pyunits.kg)
            costing_ec.electrode_material_cost.fix(2.23)

        if electrode_mat == "iron":
            # Reference for Fe cost: Anuf et al., 2022 - https://doi.org/https://doi.org/10.1016/j.jwpe.2022.103074
            costing.register_flow_type("iron", 3.41 * base_currency / pyunits.kg)
            costing_ec.electrode_material_cost.fix(3.41)

        if reactor_mat == "stainless_steel":
            # default is for PVC, so only need to change if it is stainless steel
            # PVC coeff reference: Uludag-Demirer et al., 2020 - https://doi.org/10.3390/su12072697
            # Steel coeff reference: Smith, 2005 - https://doi.org/10.1205/cherd.br.0509
            costing_ec.ec_reactor_cap_material_coeff.fix(3.4)

        blk.number_chambers_system = Param(
            initialize=3,
            mutable=True,
            units=pyunits.dimensionless,
            within=PositiveReals,
            doc="Number total chambers for system - EC chamber > flotation chamber > sedimentation chamber. All made of same material.",
        )

        blk.number_EC_reactors = Var(
            initialize=3,
            units=pyunits.dimensionless,
            bounds=(0, None),
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
            * electrode_material_cost_coeff
        )

        blk.capital_cost_power_supply_constraint = Constraint(
            expr=blk.capital_cost_power_supply
            == (ec_power_supply_base_slope * ec.power_required) * blk.number_EC_reactors
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

        blk.electricity_flow = pyunits.convert(ec.power_required, to_units=pyunits.kW)

        costing.cost_flow(
            blk.annual_electrode_replacement_mass_flow, ec.config.electrode_material
        )

        costing.cost_flow(blk.electricity_flow, "electricity")
