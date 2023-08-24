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

# Some more information about this module
__author__ = "Chenyu Wang"

import pyomo.environ as pyo
from pyomo.environ import (
    value,
    units as pyunits,
)
from pyomo.network import Arc, SequentialDecomposition

from idaes.core import (
    FlowsheetBlock,
    UnitModelCostingBlock,
)
from idaes.models.unit_models import (
    CSTR,
    Feed,
    Mixer,
    Separator,
    PressureChanger,
    Product,
)
from idaes.models.unit_models.separator import SplittingType
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale
from idaes.core.util.tables import (
    create_stream_table_dataframe,
    stream_table_dataframe_to_string,
)
from watertap.unit_models.cstr_injection import CSTR_Injection
from watertap.property_models.anaerobic_digestion.modified_adm1_properties import (
    ModifiedADM1ParameterBlock,
)
from watertap.property_models.anaerobic_digestion.adm1_properties_vapor import (
    ADM1_vaporParameterBlock,
)
from watertap.property_models.anaerobic_digestion.modified_adm1_reactions import (
    ModifiedADM1ReactionParameterBlock,
)
from watertap.property_models.activated_sludge.modified_asm2d_properties import (
    ModifiedASM2dParameterBlock,
)
from watertap.property_models.activated_sludge.modified_asm2d_reactions import (
    ModifiedASM2dReactionParameterBlock,
)
from watertap.unit_models.translators.translator_adm1_asm2d import (
    Translator_ADM1_ASM2D,
)
from watertap.unit_models.translators.translator_asm2d_adm1 import Translator_ASM2d_ADM1
from watertap.unit_models.anaerobic_digestor import AD
from watertap.unit_models.electroNP_ZO import ElectroNPZO
from watertap.unit_models.dewatering import (
    DewateringUnit,
    ActivatedSludgeModelType as dewater_type,
)
from watertap.unit_models.thickener import (
    Thickener,
    ActivatedSludgeModelType as thickener_type,
)
from watertap.core.util.initialization import check_solve
from watertap.costing import WaterTAPCosting

# Set up logger
_log = idaeslog.getLogger(__name__)


def build_flowsheet():
    m = pyo.ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props_ASM2D = ModifiedASM2dParameterBlock()
    m.fs.rxn_props_ASM2D = ModifiedASM2dReactionParameterBlock(
        property_package=m.fs.props_ASM2D
    )

    m.fs.props_ADM1 = ModifiedADM1ParameterBlock()
    m.fs.props_vap_ADM1 = ADM1_vaporParameterBlock()
    m.fs.rxn_props_ADM1 = ModifiedADM1ReactionParameterBlock(
        property_package=m.fs.props_ADM1
    )

    m.fs.costing = WaterTAPCosting()

    # Feed water stream
    m.fs.feed = Feed(property_package=m.fs.props_ASM2D)

    # Activated sludge process
    # First reactor (anoxic) - standard CSTR
    m.fs.R1 = CSTR(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Second reactor (anoxic) - standard CSTR
    m.fs.R2 = CSTR(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Third reactor (aerobic) - CSTR with injection
    m.fs.R3 = CSTR_Injection(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Fourth reactor (aerobic) - CSTR with injection
    m.fs.R4 = CSTR_Injection(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    # Fifth reactor (aerobic) - CSTR with injection
    m.fs.R5 = CSTR_Injection(
        property_package=m.fs.props_ASM2D, reaction_package=m.fs.rxn_props_ASM2D
    )
    m.fs.SP5 = Separator(
        property_package=m.fs.props_ASM2D, outlet_list=["underflow", "overflow"]
    )
    # Clarifier
    # TODO: Replace with more detailed model when available
    m.fs.CL1 = Separator(
        property_package=m.fs.props_ASM2D,
        outlet_list=["underflow", "effluent"],
        split_basis=SplittingType.componentFlow,
    )
    # Sludge purge splitter
    m.fs.SP6 = Separator(
        property_package=m.fs.props_ASM2D,
        outlet_list=["recycle", "waste"],
        split_basis=SplittingType.totalFlow,
    )

    # Product Blocks
    m.fs.Treated = Product(property_package=m.fs.props_ASM2D)
    m.fs.Sludge = Product(property_package=m.fs.props_ASM2D)

    # Thickener
    m.fs.thickener = Thickener(
        property_package=m.fs.props_ASM2D,
        activated_sludge_model=thickener_type.modified_ASM2D,
    )

    # Translators
    m.fs.translator_asm2d_adm1 = Translator_ASM2d_ADM1(
        inlet_property_package=m.fs.props_ASM2D,
        outlet_property_package=m.fs.props_ADM1,
        inlet_reaction_package=m.fs.rxn_props_ASM2D,
        outlet_reaction_package=m.fs.rxn_props_ADM1,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    m.fs.translator_adm1_asm2d = Translator_ADM1_ASM2D(
        inlet_property_package=m.fs.props_ADM1,
        outlet_property_package=m.fs.props_ASM2D,
        reaction_package=m.fs.rxn_props_ADM1,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    # Anaerobic digestor
    m.fs.AD = AD(
        liquid_property_package=m.fs.props_ADM1,
        vapor_property_package=m.fs.props_vap_ADM1,
        reaction_package=m.fs.rxn_props_ADM1,
        has_heat_transfer=True,
        has_pressure_change=False,
    )
    m.fs.AD.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    # Dewatering Unit
    m.fs.dewater = DewateringUnit(
        property_package=m.fs.props_ASM2D,
        activated_sludge_model=dewater_type.modified_ASM2D,
    )

    # ElectroNP
    m.fs.electroNP = ElectroNPZO(property_package=m.fs.props_ASM2D)
    m.fs.electroNP.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    # Costing
    m.fs.costing.cost_process()
    m.fs.costing.add_annual_water_production(
        m.fs.electroNP.properties_treated[0].flow_vol
    )
    m.fs.costing.add_LCOW(m.fs.AD.inlet.flow_vol[0])

    # Link units
    m.fs.stream1 = Arc(source=m.fs.feed.outlet, destination=m.fs.R1.inlet)
    m.fs.stream3 = Arc(source=m.fs.R1.outlet, destination=m.fs.R2.inlet)
    m.fs.stream4 = Arc(source=m.fs.R2.outlet, destination=m.fs.R3.inlet)
    m.fs.stream5 = Arc(source=m.fs.R3.outlet, destination=m.fs.R4.inlet)
    m.fs.stream6 = Arc(source=m.fs.R4.outlet, destination=m.fs.R5.inlet)
    m.fs.stream7 = Arc(source=m.fs.R5.outlet, destination=m.fs.SP5.inlet)
    m.fs.stream8 = Arc(source=m.fs.SP5.overflow, destination=m.fs.CL1.inlet)
    m.fs.stream10 = Arc(source=m.fs.CL1.effluent, destination=m.fs.Treated.inlet)
    m.fs.stream11 = Arc(source=m.fs.CL1.underflow, destination=m.fs.SP6.inlet)

    m.fs.stream_SP_thickener = Arc(
        source=m.fs.SP6.waste, destination=m.fs.thickener.inlet
    )
    m.fs.stream_thickener_translator = Arc(
        source=m.fs.thickener.underflow, destination=m.fs.translator_asm2d_adm1.inlet
    )
    m.fs.stream_translator_AD = Arc(
        source=m.fs.translator_asm2d_adm1.outlet, destination=m.fs.AD.inlet
    )
    m.fs.stream_AD_translator = Arc(
        source=m.fs.AD.liquid_outlet, destination=m.fs.translator_adm1_asm2d.inlet
    )
    m.fs.stream_translator_dewater = Arc(
        source=m.fs.translator_adm1_asm2d.outlet, destination=m.fs.dewater.inlet
    )
    m.fs.stream_dewater_electroNP = Arc(
        source=m.fs.dewater.overflow, destination=m.fs.electroNP.inlet
    )

    pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    # Oxygen concentration in reactors 3 and 4 is governed by mass transfer
    # Add additional parameter and constraints
    m.fs.R3.KLa = pyo.Var(
        initialize=3.4,
        units=pyo.units.hour**-1,
        doc="Lumped mass transfer coefficient for oxygen",
    )
    m.fs.R4.KLa = pyo.Var(
        initialize=1.8,
        units=pyo.units.hour**-1,
        doc="Lumped mass transfer coefficient for oxygen",
    )
    m.fs.S_O_eq = pyo.Param(
        default=8e-3,
        units=pyo.units.kg / pyo.units.m**3,
        mutable=True,
        doc="Dissolved oxygen concentration at equilibrium",
    )

    @m.fs.R3.Constraint(m.fs.time, doc="Mass transfer constraint for R3")
    def mass_transfer_R3(self, t):
        return pyo.units.convert(
            m.fs.R3.injection[t, "Liq", "S_O2"], to_units=pyo.units.kg / pyo.units.hour
        ) == (
            m.fs.R3.KLa
            * m.fs.R3.volume[t]
            * (m.fs.S_O_eq - m.fs.R3.outlet.conc_mass_comp[t, "S_O2"])
        )

    @m.fs.R4.Constraint(m.fs.time, doc="Mass transfer constraint for R4")
    def mass_transfer_R4(self, t):
        return pyo.units.convert(
            m.fs.R4.injection[t, "Liq", "S_O2"], to_units=pyo.units.kg / pyo.units.hour
        ) == (
            m.fs.R4.KLa
            * m.fs.R4.volume[t]
            * (m.fs.S_O_eq - m.fs.R4.outlet.conc_mass_comp[t, "S_O2"])
        )

    # Feed inlet
    m.fs.feed.flow_vol.fix(20648 * pyo.units.m**3 / pyo.units.day)
    m.fs.feed.temperature.fix(298.15 * pyo.units.K)
    m.fs.feed.pressure.fix(1 * pyo.units.atm)
    m.fs.feed.conc_mass_comp[0, "S_O2"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "S_F"].fix(30 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "S_A"].fix(20 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "S_NH4"].fix(16 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "S_NO3"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "S_PO4"].fix(3.6 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "S_I"].fix(30 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "S_N2"].fix(15 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "X_I"].fix(25 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "X_S"].fix(125 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "X_H"].fix(30 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "X_PAO"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "X_PP"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "X_PHA"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.conc_mass_comp[0, "X_AUT"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.feed.properties[0].conc_mass_comp["S_IC"].fix(
        0.07899 * pyunits.kg / pyunits.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["S_K"].fix(
        1e-6 * pyo.units.g / pyo.units.m**3
    )
    m.fs.feed.properties[0].conc_mass_comp["S_Mg"].fix(
        1e-6 * pyo.units.g / pyo.units.m**3
    )

    # Reactor sizing
    m.fs.R1.volume.fix(1000 * pyo.units.m**3)
    m.fs.R2.volume.fix(1000 * pyo.units.m**3)
    m.fs.R3.volume.fix(1333 * pyo.units.m**3)
    m.fs.R4.volume.fix(1333 * pyo.units.m**3)
    m.fs.R5.volume.fix(1333 * pyo.units.m**3)

    # Injection rates to Reactions 3, 4 and 5
    for j in m.fs.props_ASM2D.component_list:
        if j != "S_O2":
            # All components except S_O have no injection
            m.fs.R3.injection[:, :, j].fix(0)
            m.fs.R4.injection[:, :, j].fix(0)
            m.fs.R5.injection[:, :, j].fix(0)
    # Then set injections rates for O2
    m.fs.R3.outlet.conc_mass_comp[:, "S_O2"].fix(2.72e-3)
    m.fs.R4.outlet.conc_mass_comp[:, "S_O2"].fix(2.43e-3)
    m.fs.R5.outlet.conc_mass_comp[:, "S_O2"].fix(4.49e-4)

    # Set fraction of outflow from reactor 5 that goes to recycle
    m.fs.SP5.split_fraction[:, "underflow"].fix(0.1)

    # Clarifier
    # TODO: Update once more detailed model available
    m.fs.CL1.split_fraction[0, "effluent", "H2O"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_A"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_F"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_I"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_N2"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_NH4"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_NO3"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_O2"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_PO4"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_IC"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_K"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "S_Mg"].fix(0.49986)
    m.fs.CL1.split_fraction[0, "effluent", "X_AUT"].fix(0.022117)
    m.fs.CL1.split_fraction[0, "effluent", "X_H"].fix(0.021922)
    m.fs.CL1.split_fraction[0, "effluent", "X_I"].fix(0.021715)
    m.fs.CL1.split_fraction[0, "effluent", "X_PAO"].fix(0.022)
    m.fs.CL1.split_fraction[0, "effluent", "X_PHA"].fix(0.02147)
    m.fs.CL1.split_fraction[0, "effluent", "X_PP"].fix(0.02144)
    m.fs.CL1.split_fraction[0, "effluent", "X_S"].fix(0.02221)

    # Sludge purge separator
    m.fs.SP6.split_fraction[:, "recycle"].fix(0.65)

    # AD
    m.fs.AD.volume_liquid.fix(3400)
    m.fs.AD.volume_vapor.fix(300)
    m.fs.AD.liquid_outlet.temperature.fix(298.15)

    # ElectroNP
    m.fs.electroNP.energy_electric_flow_mass.fix(0.044 * pyunits.kWh / pyunits.kg)
    m.fs.electroNP.magnesium_chloride_dosage.fix(0.388)

    # Costing
    m.fs.costing.electroNP.phosphorus_recovery_value = 0

    # Check degrees of freedom
    print(degrees_of_freedom(m))
    assert degrees_of_freedom(m) == 0

    def scale_variables(m):
        for var in m.fs.component_data_objects(pyo.Var, descend_into=True):
            if "flow_vol" in var.name:
                iscale.set_scaling_factor(var, 1e2)
            if "temperature" in var.name:
                iscale.set_scaling_factor(var, 1e-1)
            if "pressure" in var.name:
                iscale.set_scaling_factor(var, 1e-4)
            if "conc_mass_comp" in var.name:
                iscale.set_scaling_factor(var, 1e1)

    scale_variables(m)
    iscale.calculate_scaling_factors(m.fs)

    # Initialize flowsheet
    # Apply sequential decomposition - 1 iteration should suffice
    seq = SequentialDecomposition()
    seq.options.select_tear_method = "heuristic"
    seq.options.tear_method = "Wegstein"
    seq.options.iterLim = 0

    def function(unit):
        unit.initialize(outlvl=idaeslog.INFO, optarg={"bound_push": 1e-2})

    seq.run(m, function)
    m.fs.costing.initialize()

    results = solve(m, tee=True)

    return m, results


def solve(blk, solver=None, checkpoint=None, tee=False, fail_flag=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    check_solve(results, checkpoint=checkpoint, logger=_log, fail_flag=fail_flag)
    return results


def display_costing(m):
    print("\n----------Capital Cost----------")
    total_capital_cost = value(
        pyunits.convert(m.fs.costing.total_capital_cost, to_units=pyunits.USD_2018)
    )
    normalized_capex = total_capital_cost / value(
        pyunits.convert(m.fs.AD.inlet.flow_vol[0], to_units=pyunits.m**3 / pyunits.hr)
    )
    print(f"Total Capital Costs: {total_capital_cost:.3f} $")
    print(f"Normalized Capital Costs: {normalized_capex:.3f} $/m3/hr")
    print("Capital Cost Breakdown")
    for u in m.fs.costing._registered_unit_costing:
        print(
            u.name,
            " : {price:0.3f} $".format(
                price=value(pyunits.convert(u.capital_cost, to_units=pyunits.USD_2018))
            ),
        )
    print("\n----------Operation Cost----------")
    total_operating_cost = value(
        pyunits.convert(
            m.fs.costing.total_operating_cost, to_units=pyunits.USD_2018 / pyunits.year
        )
    )
    print(f"Total Operating Cost: {total_operating_cost:.3f} $/year")

    opex_fraction = value(
        pyunits.convert(
            m.fs.costing.total_operating_cost, to_units=pyunits.USD_2018 / pyunits.year
        )
        / pyunits.convert(
            m.fs.AD.inlet.flow_vol[0], to_units=pyunits.m**3 / pyunits.year
        )
        / m.fs.costing.LCOW
    )
    print(f"Operating cost fraction: {opex_fraction:.3f} $ opex / $ LCOW")

    print("Operating Cost Breakdown")
    for f in m.fs.costing.used_flows:
        print(
            f.title(),
            " :    {price:0.3f} $/m3 feed".format(
                price=value(
                    pyunits.convert(
                        m.fs.costing.aggregate_flow_costs[f],
                        to_units=pyunits.USD_2018 / pyunits.year,
                    )
                    / pyunits.convert(
                        m.fs.AD.inlet.flow_vol[0],
                        to_units=pyunits.m**3 / pyunits.year,
                    )
                )
            ),
        )

    print("\n----------Energy----------")

    electricity_intensity = value(
        pyunits.convert(
            m.fs.costing.aggregate_flow_electricity / m.fs.AD.inlet.flow_vol[0],
            to_units=pyunits.kWh / pyunits.m**3,
        )
    )
    print(f"Electricity Intensity: {electricity_intensity:.4f} kWh/m3")

    print("\n----------Levelized Cost----------")
    LCOW = value(
        pyunits.convert(m.fs.costing.LCOW, to_units=pyunits.USD_2018 / pyunits.m**3)
    )
    print(f"Levelized Cost of Water: {LCOW:.3f} $/m^3")


if __name__ == "__main__":
    m, results = build_flowsheet()

    stream_table = create_stream_table_dataframe(
        {
            "Feed": m.fs.feed.outlet,
            "ASP outlet": m.fs.R5.outlet,
            "thickener outlet": m.fs.thickener.underflow,
            "AD liquid inlet": m.fs.AD.inlet,
            "AD liquid outlet": m.fs.AD.liquid_outlet,
            "AD vapor outlet": m.fs.AD.vapor_outlet,
            "dewater outlet": m.fs.dewater.overflow,
            "ElectroNP treated": m.fs.electroNP.treated,
            "ElectroNP byproduct": m.fs.electroNP.byproduct,
        },
        time_point=0,
    )
    print(stream_table_dataframe_to_string(stream_table))
    display_costing(m)