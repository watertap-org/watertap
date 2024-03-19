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
import os
import idaes.logger as idaeslog
from pyomo.environ import (
    ConcreteModel,
    Block,
    Expression,
    TransformationFactory,
    units as pyunits,
)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.models.unit_models import (
    Product,
    Translator,
)
import idaes.core.util.scaling as iscale
from idaes.core import UnitModelCostingBlock

from watertap.core.util.initialization import assert_degrees_of_freedom, check_solve

from watertap.core.wt_database import Database
import watertap.core.zero_order_properties as prop_ZO
from watertap.unit_models.zero_order import (
    FeedZO,
    PumpElectricityZO,
    NanofiltrationZO,
)
from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
    DiffusivityCalculation,
)
from watertap.unit_models.gac import (
    GAC,
)
from watertap.costing.zero_order_costing import ZeroOrderCosting

from watertap.core.util.model_diagnostics.infeasible import *
from watertap.costing import WaterTAPCosting


# Set up logger
_log = idaeslog.getLogger(__name__)


def main():
    m = build(include_gac=False)

    set_operating_conditions(m)

    assert_degrees_of_freedom(m, 0)
    assert_units_consistent(m)

    initialize_system(m)

    results = solve(m, checkpoint="solve flowsheet after initializing system")

    display_results(m)

    add_costing(m)
    assert_degrees_of_freedom(m, 0)

    results = solve(m, checkpoint="solve flowsheet after costing")

    display_costing(m)
    return m, results


def build(include_gac=False):
    # flowsheet set up
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.prop = prop_ZO.WaterParameterBlock(solute_list=["dye", "tds"])

    # unit model
    m.fs.feed = FeedZO(property_package=m.fs.prop)

    # define block to integrate with dye_desalination_withRO
    dye_sep = m.fs.dye_separation = Block()

    dye_sep.P1 = PumpElectricityZO(
        property_package=m.fs.prop, database=m.db, process_subtype="default"
    )

    dye_sep.nanofiltration = NanofiltrationZO(
        property_package=m.fs.prop, database=m.db, process_subtype="rHGO_dye_rejection"
    )

    m.fs.permeate = Product(property_package=m.fs.prop)

    if include_gac == True:
        # GAC
        m.fs.prop_gac = MCASParameterBlock(
            material_flow_basis="mass",
            ignore_neutral_charge=True,
            solute_list=["tds", "dye"],
            mw_data={
                "H2O": 0.018,
                "tds": 0.05844,
                "dye": 0.696665,  # molecular weight of congo red dye
            },
            diffus_calculation=DiffusivityCalculation.none,
            diffusivity_data={("Liq", "tds"): 1e-09, ("Liq", "dye"): 2e-10},
        )
        m.fs.gac = GAC(
            property_package=m.fs.prop_gac,
            film_transfer_coefficient_type="calculated",
            surface_diffusion_coefficient_type="fixed",
            target_species={"dye"},
        )
        m.fs.adsorbed_dye = Product(property_package=m.fs.prop_gac)
        m.fs.treated = Product(property_package=m.fs.prop_gac)

        m.fs.tb_nf_gac = Translator(
            inlet_property_package=m.fs.prop, outlet_property_package=m.fs.prop_gac
        )

        @m.fs.tb_nf_gac.Constraint(["H2O", "dye", "tds"])
        def eq_flow_mass_comp(blk, j):
            if j == "dye":
                return (
                    blk.properties_in[0].flow_mass_comp["dye"]
                    == blk.properties_out[0].flow_mass_phase_comp["Liq", "dye"]
                )
            elif j == "tds":
                return (
                    blk.properties_in[0].flow_mass_comp["tds"]
                    == blk.properties_out[0].flow_mass_phase_comp["Liq", "tds"]
                )
            else:
                return (
                    blk.properties_in[0].flow_mass_comp["H2O"]
                    == blk.properties_out[0].flow_mass_phase_comp["Liq", "H2O"]
                )

    else:
        m.fs.dye_retentate = Product(property_package=m.fs.prop)

    # connections
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=dye_sep.P1.inlet)
    m.fs.s02 = Arc(source=dye_sep.P1.outlet, destination=dye_sep.nanofiltration.inlet)
    m.fs.s03 = Arc(
        source=dye_sep.nanofiltration.treated, destination=m.fs.permeate.inlet
    )
    if hasattr(m.fs, "gac"):
        m.fs.s04 = Arc(
            source=dye_sep.nanofiltration.byproduct, destination=m.fs.tb_nf_gac.inlet
        )
        m.fs.s05 = Arc(source=m.fs.tb_nf_gac.outlet, destination=m.fs.gac.inlet)
        # TODO: Recycle treated stream back to the feed via a mixer
        m.fs.s06 = Arc(source=m.fs.gac.outlet, destination=m.fs.treated.inlet)
        m.fs.s07 = Arc(source=m.fs.gac.adsorbed, destination=m.fs.adsorbed_dye.inlet)
    else:
        m.fs.s04 = Arc(
            source=dye_sep.nanofiltration.byproduct,
            destination=m.fs.dye_retentate.inlet,
        )

    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    iscale.calculate_scaling_factors(m)

    return m


def set_operating_conditions(m):
    dye_sep = m.fs.dye_separation
    # feed
    flow_vol = 280 / 3600 * pyunits.m**3 / pyunits.s
    conc_mass_dye = 0.2 * pyunits.kg / pyunits.m**3
    conc_mass_tds = 2 * pyunits.kg / pyunits.m**3

    m.fs.feed.flow_vol[0].fix(flow_vol)
    m.fs.feed.conc_mass_comp[0, "dye"].fix(conc_mass_dye)
    m.fs.feed.conc_mass_comp[0, "tds"].fix(conc_mass_tds)
    solve(m.fs.feed, checkpoint="solve feed block")

    # nanofiltration
    dye_sep.nanofiltration.load_parameters_from_database(use_default_removal=True)

    # pump
    dye_sep.P1.load_parameters_from_database(use_default_removal=True)
    dye_sep.P1.applied_pressure.fix(
        dye_sep.nanofiltration.applied_pressure.get_values()[0]
    )
    dye_sep.P1.eta_pump.fix(0.75)  # pump efficiency [-]
    dye_sep.P1.lift_height.unfix()

    if hasattr(m.fs, "gac"):
        m.fs.tb_nf_gac.properties_out[0].temperature.fix(298.15)
        m.fs.tb_nf_gac.properties_out[0].pressure.fix(101325)

        m.fs.gac.freund_k.fix(10)
        m.fs.gac.freund_ninv.fix(0.9)
        m.fs.gac.shape_correction_factor.fix()
        m.fs.gac.ds.fix(5e-13)
        # gac particle specifications
        m.fs.gac.particle_dens_app.fix(750)
        m.fs.gac.particle_dia.fix(0.001)
        # adsorber bed specifications
        m.fs.gac.ebct.fix(600)
        m.fs.gac.bed_voidage.fix(0.4)
        m.fs.gac.bed_length.fix(6)
        # design spec
        m.fs.gac.conc_ratio_replace.fix(0.50)
        # parameters
        m.fs.gac.a0.fix(3.68421)
        m.fs.gac.a1.fix(13.1579)
        m.fs.gac.b0.fix(0.784576)
        m.fs.gac.b1.fix(0.239663)
        m.fs.gac.b2.fix(0.484422)
        m.fs.gac.b3.fix(0.003206)
        m.fs.gac.b4.fix(0.134987)

        iscale.constraint_scaling_transform(m.fs.gac.eq_mass_adsorbed["dye"], 1e-2)
    else:
        pass

    return


def initialize_system(m):
    seq = SequentialDecomposition()
    seq.options.tear_set = []
    seq.options.iterLim = 1
    seq.run(m, lambda u: u.initialize())


def solve(blk, solver=None, checkpoint=None, tee=False, fail_flag=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    check_solve(results, checkpoint=checkpoint, logger=_log, fail_flag=fail_flag)
    return results


def add_costing(m):
    # initialize block
    dye_sep = m.fs.dye_separation

    source_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "dye_desalination_global_costing.yaml",
    )

    # zero order costing
    m.fs.zo_costing = ZeroOrderCosting(case_study_definition=source_file)

    costing_kwargs = {"flowsheet_costing_block": m.fs.zo_costing}

    # create costing blocks
    dye_sep.nanofiltration.costing = UnitModelCostingBlock(**costing_kwargs)
    dye_sep.P1.costing = UnitModelCostingBlock(**costing_kwargs)
    m.fs.zo_costing.pump_electricity.pump_cost["default"].fix(76)

    # aggregate unit level costs
    m.fs.zo_costing.cost_process()

    # create system level cost metrics
    m.fs.brine_recovery_revenue = Expression(
        expr=(
            m.fs.zo_costing.utilization_factor
            * (
                m.fs.zo_costing.brine_recovery_cost
                * pyunits.convert(
                    m.fs.permeate.properties[0].flow_vol,
                    to_units=pyunits.m**3 / m.fs.zo_costing.base_period,
                )
            )
        ),
        doc="Revenue from recovering saline brine/ NF permeate",
    )
    if hasattr(m.fs, "gac"):
        m.fs.dye_disposal_cost = Expression(
            expr=(
                m.fs.zo_costing.utilization_factor
                * m.fs.zo_costing.dye_disposal_cost
                * pyunits.convert(
                    m.fs.adsorbed_dye.properties[0].flow_vol,
                    to_units=pyunits.m**3 / m.fs.zo_costing.base_period,
                )
            ),
            doc="Cost of disposing of dye waste",
        )
    else:
        m.fs.dye_disposal_cost = Expression(
            expr=(
                m.fs.zo_costing.utilization_factor
                * m.fs.zo_costing.dye_disposal_cost
                * pyunits.convert(
                    m.fs.dye_retentate.properties[0].flow_vol,
                    to_units=pyunits.m**3 / m.fs.zo_costing.base_period,
                )
            ),
            doc="Cost of disposing of dye waste",
        )

    # combine results for system level costs - to be the same syntax as dye_desalination_withRO
    @m.fs.Expression(doc="Total capital cost")
    def total_capital_cost(b):
        return pyunits.convert(
            m.fs.zo_costing.total_capital_cost, to_units=pyunits.USD_2020
        )

    @m.fs.Expression(doc="Total operating cost")
    def total_operating_cost(b):
        return pyunits.convert(
            b.zo_costing.total_fixed_operating_cost,
            to_units=pyunits.USD_2020 / pyunits.year,
        ) + pyunits.convert(
            b.zo_costing.total_variable_operating_cost,
            to_units=pyunits.USD_2020 / pyunits.year,
        )

    @m.fs.Expression(doc="Total cost of brine recovery and dye disposal")
    def total_externalities(b):
        return pyunits.convert(
            m.fs.brine_recovery_revenue - m.fs.dye_disposal_cost,
            to_units=pyunits.USD_2020 / pyunits.year,
        )

    @m.fs.Expression(
        doc="Levelized cost of treatment with respect to volumetric feed flowrate"
    )
    def LCOT(b):
        return (
            b.total_capital_cost * b.zo_costing.capital_recovery_factor
            + b.total_operating_cost
            - b.total_externalities
        ) / (
            pyunits.convert(
                b.feed.properties[0].flow_vol,
                to_units=pyunits.m**3 / pyunits.year,
            )
            * b.zo_costing.utilization_factor
        )

    assert_units_consistent(m)
    m.fs.zo_costing.initialize()
    return


def display_results(m):
    unit_list = ["P1", "nanofiltration"]
    for u in unit_list:
        m.fs.dye_separation.component(u).report()

    if hasattr(m.fs, "gac"):
        m.fs.costing = WaterTAPCosting()
        m.fs.costing.base_currency = pyunits.USD_2021

        m.fs.gac.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={"contactor_type": "gravity"},
        )

        m.fs.costing.gac_gravity.regen_frac.fix(0.7)
        m.fs.costing.gac_gravity.num_contactors_op.fix(1)
        m.fs.costing.gac_gravity.num_contactors_redundant.fix(1)

        m.fs.costing.cost_process()

        m.fs.gac.costing.initialize()
        m.fs.costing.initialize()
        gac_capex = value(m.fs.gac.costing.capital_cost / 1e6)
        print(f"GAC capital expenses: {gac_capex}M$")
    else:
        pass


def display_costing(m):
    print("\nUnit Capital Costs")
    for u in m.fs.zo_costing._registered_unit_costing:
        print(
            u.name,
            " :   ",
            value(pyunits.convert(u.capital_cost, to_units=pyunits.USD_2020)),
            "$",
        )

    print("\nSystem Costs")
    total_capital_cost = value(
        pyunits.convert(m.fs.total_capital_cost, to_units=pyunits.MUSD_2020)
    )
    print(f"Total Capital Costs: {total_capital_cost:.4f} M$")

    total_operating_cost = value(
        pyunits.convert(
            m.fs.total_operating_cost, to_units=pyunits.MUSD_2020 / pyunits.year
        )
    )
    print(f"Total Operating Costs: {total_operating_cost:.4f} M$/year")

    total_externalities = value(
        pyunits.convert(
            m.fs.total_externalities, to_units=pyunits.MUSD_2020 / pyunits.year
        )
    )
    brr = value(
        pyunits.convert(
            m.fs.brine_recovery_revenue, to_units=pyunits.USD_2020 / pyunits.year
        )
    )
    print(f"Brine recovery revenue: {brr: .4f} M$/year")
    ddc = value(
        pyunits.convert(
            m.fs.dye_disposal_cost, to_units=pyunits.USD_2020 / pyunits.year
        )
    )
    print(f"Dye disposal cost: {ddc: .2f} $/year")
    print(f"Brine recovery revenue: {brr: .2f} $/year")
    print(f"Total Externalities: {total_externalities:.4f} M$/year")

    levelized_cost_treatment = value(
        pyunits.convert(m.fs.LCOT, to_units=pyunits.USD_2020 / pyunits.m**3)
    )
    print(
        f"Levelized Cost of Treatment (LCOT): {levelized_cost_treatment:.2f} $/m3 feed"
    )


if __name__ == "__main__":
    m, results = main()
