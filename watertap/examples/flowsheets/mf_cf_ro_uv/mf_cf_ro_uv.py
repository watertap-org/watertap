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
from pyomo.environ import (
    ConcreteModel,
    value,
    TransformationFactory,
    units as pyunits,
    Block,
    Constraint,
    assert_optimal_termination,
)
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock, UnitModelBlockData
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import (
    propagate_state,
    fix_state_vars,
    revert_state_vars,
)
from idaes.core.util import DiagnosticsToolbox

from idaes.core.util.exceptions import ConfigurationError
from idaes.models.unit_models.translator import Translator
from idaes.models.unit_models import Mixer, Separator, Product, Feed
from idaes.models.unit_models.mixer import MomentumMixingType
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core import UnitModelCostingBlock

from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.NaCl_T_dep_prop_pack import NaClParameterBlock

from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.unit_models.reverse_osmosis_1D import ReverseOsmosis1D
from watertap.unit_models.pressure_exchanger import PressureExchanger
from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice
from watertap.examples.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery import (
    ERDtype,
    erd_type_not_found,
)
from watertap.unit_models.uv_aop import Ultraviolet0D
from watertap.core.util.initialization import assert_degrees_of_freedom, check_solve
from watertap.core.wt_database import Database
from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
    MaterialFlowBasis,
)
from watertap.unit_models.zero_order import (
    MicroFiltrationZO,
    CartridgeFiltrationZO,
    UVZO,
    UVAOPZO,
    ChemicalAdditionZO,
    StaticMixerZO,
    StorageTankZO,
    MediaFiltrationZO,
    BackwashSolidsHandlingZO,
)
from watertap.costing.zero_order_costing import ZeroOrderCosting
from watertap.costing import WaterTAPCosting
from idaes.core.util.misc import StrEnum

__author__ = "Adam Atia"
# Set up logger
_log = idaeslog.getLogger(__name__)


class uvdimension(StrEnum):
    none = "none"
    zo = "zo"
    zero_d = "0d"
    one_d = "1d"


def main(
    ro_props="seawater",
    ro_dimension="0d",
    erd_config=ERDtype.no_ERD,
    uvdimension=uvdimension.zo,
    has_aop=False,
    diagnostics_active=False,
):
    m = build(
        ro_props=ro_props,
        ro_dimension=ro_dimension,
        erd_config=erd_config,
        uvdimension=uvdimension,
        has_aop=has_aop,
    )

    set_operating_conditions(m)

    dt = DiagnosticsToolbox(m)
    dt.report_structural_issues()

    if diagnostics_active:

        try:
            assert_degrees_of_freedom(m, 0)
            assert_units_consistent(m)
            initialize_system(m)
            if diagnostics_active:
                dt.report_numerical_issues()
            results = solve(m, checkpoint="solve flowsheet after initializing system")
            assert_optimal_termination(results)
            display_results(m)

            # TODO:
            # add_costing(m)
            # assert_degrees_of_freedom(m, 0)
            # m.fs.costing.initialize()

            # results = solve(m, checkpoint="solve flowsheet after costing")

            # display_results(m)
        except:
            return _, _, dt
    else:
        assert_degrees_of_freedom(m, 0)
        assert_units_consistent(m)
        initialize_system(m)

        results = solve(m, checkpoint="solve flowsheet after initializing system")
        assert_optimal_termination(results)
        display_results(m)

        # TODO:
        # add_costing(m)
        # assert_degrees_of_freedom(m, 0)
        # m.fs.costing.initialize()

        # results = solve(m, checkpoint="solve flowsheet after costing")

    return m, results, dt


def build(ro_props, ro_dimension, erd_config, uvdimension, has_aop):
    """
    ro_prop_model: choose between "NaCl" and "Seawater" prop models for RO
    """
    # flowsheet set up
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # Property Models: MCAS and either NaCl/Seawater

    # Choose between NaCl or Seawater prop models for RO
    m.fs.ro_props, m.fs.ro_ion = get_ro_props(ro_props)

    # Use MCAS for whole flowsheet, except RO
    m.fs.mcas_props = MCASParameterBlock(
        solute_list=[m.fs.ro_ion, "tss"],
        diffusivity_data={("Liq", m.fs.ro_ion): 1e-9, ("Liq", "tss"): 1e-9},
        mw_data={m.fs.ro_ion: None, "tss": None},
        material_flow_basis=MaterialFlowBasis.mass,
        ignore_neutral_charge=True,
    )

    # Add Database for initially parameterization of ZO models
    m.db = Database()

    # Feed block
    m.fs.feed = Feed(property_package=m.fs.mcas_props)

    # Microfiltration Pump
    m.fs.mf_pump = Pump(property_package=m.fs.mcas_props)

    # Microfiltration
    m.fs.mf = MicroFiltrationZO(property_package=m.fs.mcas_props, database=m.db)

    # Cartridge Filtration Pump
    m.fs.cf_pump = Pump(property_package=m.fs.mcas_props)

    # CF
    m.fs.cf = CartridgeFiltrationZO(property_package=m.fs.mcas_props, database=m.db)

    # Translate MCAS to RO property model
    m.fs.mcas_to_ro_translator = Translator(
        inlet_property_package=m.fs.mcas_props, outlet_property_package=m.fs.ro_props
    )

    @m.fs.mcas_to_ro_translator.Constraint(m.fs.ro_props.component_list)
    def eq_flow_mass_comp(blk, j):
        if j.lower() == "tss":
            Constraint.Skip
        else:
            return (
                blk.properties_in[0].flow_mass_phase_comp["Liq", j]
                == blk.properties_out[0].flow_mass_phase_comp["Liq", j]
            )

    # RO Train =====================================================================
    # High-pressure RO pump
    m.fs.hp_pump = Pump(property_package=m.fs.ro_props)

    # --- Reverse Osmosis Block ---
    m.fs.RO = get_ro_model(dimension=ro_dimension, ro_props=m.fs.ro_props)

    # --- ERD blocks ---
    if erd_config == ERDtype.pressure_exchanger:
        m.fs.feed_to_hp_and_erd_splitter = Separator(
            property_package=m.fs.ro_props, outlet_list=["hp_pump", "erd"]
        )

        m.fs.erd = PressureExchanger(property_package=m.fs.ro_props)
        m.fs.booster_pump = Pump(property_package=m.fs.ro_props)
        m.fs.hp_and_booster_mixer = Mixer(
            property_package=m.fs.ro_props,
            momentum_mixing_type=MomentumMixingType.equality,
            inlet_list=["hp_pump", "booster_pump"],
        )

    elif erd_config == ERDtype.pump_as_turbine:
        # add energy recovery turbine block
        m.fs.erd = EnergyRecoveryDevice(property_package=m.fs.ro_props)

    elif erd_config == ERDtype.no_ERD:
        pass
    else:
        erd_type_not_found(erd_config)
    # ===============================================================================

    # Translate RO to MCAS property model
    m.fs.ro_to_mcas_translator = Translator(
        inlet_property_package=m.fs.ro_props, outlet_property_package=m.fs.mcas_props
    )

    @m.fs.ro_to_mcas_translator.Constraint(m.fs.mcas_props.component_list)
    def eq_flow_mass_comp(blk, j):
        if j.lower() == "tss":
            return (
                m.fs.cf.treated.flow_mass_phase_comp[0, "Liq", j]
                == blk.properties_out[0].flow_mass_phase_comp["Liq", j]
            )
        else:
            return (
                blk.properties_in[0].flow_mass_phase_comp["Liq", j]
                == blk.properties_out[0].flow_mass_phase_comp["Liq", j]
            )

    # UV
    # TODO: Add UV as an option: None, UV, UV-AOP, UVZO, UVAOPZO?
    m.fs.uv = get_uv_model(
        m, uv_props=m.fs.mcas_props, dimension=uvdimension, has_aop=has_aop
    )

    # Product blocks for permeate and disposal
    m.fs.product_water = Product(property_package=m.fs.mcas_props)
    m.fs.waste_brine = Product(property_package=m.fs.ro_props)
    # TODO: consider adding mixer to collect MF, CF, and RO waste streams and send to drain (Product block)

    # connections
    m.fs.feed_to_mf_pump = Arc(source=m.fs.feed.outlet, destination=m.fs.mf_pump.inlet)
    m.fs.mf_pump_to_mf = Arc(source=m.fs.mf_pump.outlet, destination=m.fs.mf.inlet)
    m.fs.mf_to_cf_pump = Arc(source=m.fs.mf.treated, destination=m.fs.cf_pump.inlet)
    m.fs.cf_pump_to_cf = Arc(source=m.fs.cf_pump.outlet, destination=m.fs.cf.inlet)
    m.fs.cf_to_translator = Arc(
        source=m.fs.cf.treated, destination=m.fs.mcas_to_ro_translator.inlet
    )

    if erd_config == ERDtype.pressure_exchanger:
        raise NotImplementedError(
            f"While arc connections are ready, setting up the square problem (setting DOF=0) has not been completed for erd_config={erd_config}"
        )
        # TODO: finalize arc connections and square prob setup with this erd_config
        m.fs.translator_to_hp_erd_splitter = Arc(
            source=m.fs.mcas_to_ro_translator.outlet,
            destination=m.fs.feed_to_hp_and_erd_splitter.inlet,
        )
        m.fs.splitter_to_hp_pump = Arc(
            source=m.fs.feed_to_hp_and_erd_splitter.hp_pump,
            destination=m.fs.hp_pump.inlet,
        )
        m.fs.splitter_to_erd = Arc(
            source=m.fs.feed_to_hp_and_erd_splitter.erd,
            destination=m.fs.erd.low_pressure_inlet,
        )
        m.fs.erd_to_booster_pump = Arc(
            source=m.fs.erd.low_pressure_outlet, destination=m.fs.booster_pump.inlet
        )
        m.fs.booster_pump_to_mixer = Arc(
            source=m.fs.booster_pump.outlet,
            destination=m.fs.hp_and_booster_mixer.booster_pump,
        )
        m.fs.hp_pump_to_mixer = Arc(
            source=m.fs.hp_pump.outlet, destination=m.fs.hp_and_booster_mixer.hp_pump
        )
        m.fs.mixer_to_RO = Arc(
            source=m.fs.hp_and_booster_mixer.outlet, destination=m.fs.RO.inlet
        )
        m.fs.RO_brine_to_erd = Arc(
            source=m.fs.RO.retentate, destination=m.fs.erd.high_pressure_inlet
        )
        m.fs.erd_to_waste = Arc(
            source=m.fs.erd.high_pressure_outlet, destination=m.fs.waste_brine.inlet
        )
    elif erd_config == ERDtype.pump_as_turbine:
        raise NotImplementedError(
            f"While arc connections are ready, setting up the square problem (setting DOF=0) has not been completed for erd_config={erd_config}"
        )
        # TODO: finalize arc connections and square prob setup with this erd_config

        m.fs.translator_to_hp_pump = Arc(
            source=m.fs.mcas_to_ro_translator.outlet, destination=m.fs.hp_pump.inlet
        )
        m.fs.hp_pump_to_RO = Arc(source=m.fs.hp_pump.outlet, destination=m.fs.RO.inlet)
        m.fs.RO_brine_to_erd = Arc(source=m.fs.RO.retentate, destination=m.fs.erd.inlet)
        m.fs.erd_to_waste = Arc(
            source=m.fs.erd.outlet, destination=m.fs.waste_brine.inlet
        )
    elif erd_config == ERDtype.no_ERD:

        m.fs.translator_to_hp_pump = Arc(
            source=m.fs.mcas_to_ro_translator.outlet, destination=m.fs.hp_pump.inlet
        )
        m.fs.hp_pump_to_RO = Arc(source=m.fs.hp_pump.outlet, destination=m.fs.RO.inlet)
        m.fs.RO_brine_to_waste = Arc(
            source=m.fs.RO.retentate, destination=m.fs.waste_brine.inlet
        )
    else:
        # this case should be caught in the previous conditional
        erd_type_not_found(erd_config)

    m.fs.ro_permeate_to_translator = Arc(
        source=m.fs.RO.permeate, destination=m.fs.ro_to_mcas_translator.inlet
    )
    # TODO: UV conditionals to determine RO permeate connections
    if uvdimension == uvdimension.none:
        m.fs.translator_to_product = Arc(
            source=m.fs.ro_to_mcas_translator.outlet,
            destination=m.fs.product_water.inlet,
        )
    elif uvdimension == uvdimension.zo:
        m.fs.translator_to_uv = Arc(
            source=m.fs.ro_to_mcas_translator.outlet, destination=m.fs.uv.inlet
        )
        m.fs.uv_to_product = Arc(
            source=m.fs.uv.treated, destination=m.fs.product_water.inlet
        )
    elif uvdimension == uvdimension.zero_d:
        m.fs.translator_to_uv = Arc(
            source=m.fs.ro_to_mcas_translator.outlet, destination=m.fs.uv.inlet
        )
        m.fs.uv_to_product = Arc(
            source=m.fs.uv.outlet, destination=m.fs.product_water.inlet
        )
    else:
        raise NotImplementedError("This config isn't available.")
    # Apply connections
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling
    # set default property values
    m.fs.ro_props.set_default_scaling(
        "flow_mass_phase_comp", 1e-2, index=("Liq", "H2O")
    )
    m.fs.ro_props.set_default_scaling(
        "flow_mass_phase_comp", 1e-1, index=("Liq", m.fs.ro_ion)
    )
    # if m.fs.uv is not None:
    #     iscale.set_scaling_factor(m.fs.uv.control_volume.properties_in[0].flow_mass_phase_comp['Liq', 'tss'], 1e3)
    # set unit model values
    # iscale.set_scaling_factor(m.fs.hp_pump.control_volume.work, 1e-5)
    # iscale.set_scaling_factor(m.fs.RO.area, 1e-4)
    # calculate and propagate scaling factors
    iscale.calculate_scaling_factors(m)
    return m


def set_operating_conditions(m):
    # ---specifications---
    # feed
    flow_vol = 0.3092 * pyunits.m**3 / pyunits.s
    conc_mass_tds = 35 * pyunits.kg / pyunits.m**3
    conc_mass_tss = 0.03 * pyunits.kg / pyunits.m**3
    temperature = 298 * pyunits.K
    pressure = 1e5 * pyunits.Pa
    m.fs.feed.temperature[0].fix(temperature)
    m.fs.feed.pressure[0].fix(pressure)
    m.fs.feed.properties.calculate_state(
        var_args={
            ("conc_mass_phase_comp", ("Liq", m.fs.ro_ion)): conc_mass_tds,
            ("conc_mass_phase_comp", ("Liq", "tss")): conc_mass_tss,
            ("flow_vol_phase", "Liq"): flow_vol,
        },
        hold_state=True,
    )
    # TODO: add scaling at least for feed props before solve
    # solve(m.fs.feed, checkpoint="solve feed block")
    m.fs.feed.display()

    # ---pretreatment---
    # TODO: add option to eliminate underlying fixed energy calculations and shift to pump
    # mf pump
    m.fs.mf_pump.efficiency_pump.fix(0.8)
    m.fs.mf_pump.control_volume.properties_out[0].pressure.fix(2e5)

    # microfiltration
    m.db.get_unit_operation_parameters("microfiltration")
    m.fs.mf.load_parameters_from_database(use_default_removal=True)

    # cf pump
    m.fs.cf_pump.efficiency_pump.fix(0.8)
    m.fs.cf_pump.control_volume.properties_out[0].pressure.fix(2e5)

    # cartridge filtration
    m.db.get_unit_operation_parameters("cartridge_filtration")
    m.fs.cf.load_parameters_from_database(use_default_removal=True)

    m.fs.mcas_to_ro_translator.outlet.pressure[0].fix(pressure)
    m.fs.mcas_to_ro_translator.outlet.temperature[0].fix(temperature)

    # hp pump
    m.fs.hp_pump.efficiency_pump.fix(0.8)
    m.fs.hp_pump.control_volume.properties_out[0].pressure.fix(70e5)

    m.fs.RO.A_comp.fix(4.2e-12)  # membrane water permeability coefficient [m/s-Pa]
    m.fs.RO.B_comp.fix(3.5e-8)  # membrane salt permeability coefficient [m/s]
    m.fs.RO.feed_side.channel_height.fix(1e-3)  # channel height in membrane stage [m]
    m.fs.RO.feed_side.spacer_porosity.fix(0.97)  # spacer porosity in membrane stage [-]
    m.fs.RO.permeate.pressure[0].fix(101325)  # atmospheric pressure [Pa]
    m.fs.RO.width.fix(1000)  # stage width [m]
    m.fs.RO.area.fix(flow_vol * 4.5e4 * pyunits.s / pyunits.m)

    m.fs.ro_to_mcas_translator.outlet.pressure[0].fix(pressure)
    m.fs.ro_to_mcas_translator.outlet.temperature[0].fix(temperature)

    if hasattr(m.fs.uv, "_tech_type"):
        m.fs.uv.load_parameters_from_database(use_default_removal=True)
    elif isinstance(m.fs.uv, UnitModelBlockData):
        m.fs.uv.uv_intensity.fix(1 * pyunits.mW / pyunits.cm**2)
        m.fs.uv.exposure_time.fix(500 * pyunits.s)
        m.fs.uv.inactivation_rate["Liq", "tss"].fix(2.3 * pyunits.cm**2 / pyunits.J)
        # m.fs.uv.outlet.pressure[0].fix(101325)
        m.fs.uv.outlet.temperature[0].fix(temperature)

        m.fs.uv.electrical_efficiency_phase_comp[0, "Liq", "tss"].fix(
            0.1 * pyunits.kWh / pyunits.m**3
        )
        m.fs.uv.lamp_efficiency.fix(0.8)
        if m.fs.uv.config.has_aop:
            m.fs.uv.second_order_reaction_rate_constant["Liq", "tss"].fix(
                3.3e8 * pyunits.M**-1 * pyunits.s**-1
            )
            # TODO: hydrogen_peroxide_conc should be generalized to oxidant_conc or something of the like
            m.fs.uv.hydrogen_peroxide_conc.fix(5.05e-13 * pyunits.M)


def initialize_system(m):

    m.fs.feed.initialize()
    propagate_state(m.fs.feed_to_mf_pump)
    m.fs.mf_pump.initialize()
    propagate_state(m.fs.mf_pump_to_mf)
    m.fs.mf.initialize()

    propagate_state(m.fs.mf_to_cf_pump)
    m.fs.cf_pump.initialize()

    propagate_state(m.fs.cf_pump_to_cf)
    m.fs.cf.initialize()

    propagate_state(m.fs.cf_to_translator)
    m.fs.mcas_to_ro_translator.initialize()

    propagate_state(m.fs.translator_to_hp_pump)
    m.fs.hp_pump.initialize()

    propagate_state(m.fs.hp_pump_to_RO)

    try:
        m.fs.RO.initialize()
    except:
        pass
    propagate_state(m.fs.ro_permeate_to_translator)
    m.fs.ro_to_mcas_translator.initialize()

    if hasattr(m.fs, "translator_to_uv"):
        propagate_state(m.fs.translator_to_uv)
        m.fs.uv.initialize(outlvl=idaeslog.DEBUG)
        propagate_state(m.fs.uv_to_product)
    else:
        propagate_state(m.fs.translator_to_product)


def solve(blk, solver=None, checkpoint=None, tee=False, fail_flag=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    check_solve(results, checkpoint=checkpoint, logger=_log, fail_flag=fail_flag)
    return results


def get_uv_model(m, uv_props, dimension, has_aop):
    if dimension == "zo" and not has_aop:
        return UVZO(
            property_package=uv_props,
            database=m.db,
        )
    if dimension == "zo" and has_aop:
        return UVAOPZO(
            property_package=uv_props,
            database=m.db,
        )
    elif dimension == "0d" and not has_aop:
        return Ultraviolet0D(property_package=uv_props, target_species=["tss"])
    elif dimension == "0d" and has_aop:
        return Ultraviolet0D(
            property_package=uv_props, target_species=["tss"], has_aop=True
        )
    elif dimension == "none":
        return None
    else:
        if dimension not in uvdimension:
            raise ConfigurationError(
                f"Either 'none', 'zo', or '0d' should be provided for uvdimension instead of {dimension}"
            )
        elif has_aop is not bool:
            raise ConfigurationError(f"has_aop should be True or False, not {has_aop}")
        else:
            raise ConfigurationError("There's something wrong...but what?")


def get_ro_props(ro_props):
    if ro_props.lower() == "nacl":
        return NaClParameterBlock(), "NaCl"
    elif ro_props.lower() == "seawater":
        return SeawaterParameterBlock(), "TDS"
    else:
        raise ConfigurationError(
            f"Either 'nacl' or 'seawater' should be provided for ro_props. Instead, {ro_props} was provided."
        )


def get_ro_model(dimension, ro_props):
    if dimension == "0d":
        return ReverseOsmosis0D(
            property_package=ro_props,
            has_pressure_change=True,
            pressure_change_type=PressureChangeType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            has_full_reporting=True,
        )
    elif dimension == "1d":
        return ReverseOsmosis1D(
            property_package=ro_props,
            has_pressure_change=True,
            pressure_change_type=PressureChangeType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
        )
    else:
        raise ConfigurationError(
            f"Either '0d' or '1d' should be provided for dimension instead of {dimension}"
        )


def display_results(m):
    m.fs.report()
    # call report() on every unit in the flowsheet:
    for block in m.fs.component_objects(Block, descend_into=True):
        if isinstance(block, UnitModelBlockData):
            block.report()


def add_costing(m):
    # TODO: add costing
    # process costing and add system level metrics
    # m.fs.costing.cost_process()
    # m.fs.costing.add_annual_water_production(m.fs.product.properties[0].flow_vol)
    # m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol)
    # m.fs.costing.add_specific_energy_consumption(m.fs.product.properties[0].flow_vol)
    # m.fs.costing.add_specific_electrical_carbon_intensity(
    #     m.fs.product.properties[0].flow_vol
    # )
    pass


if __name__ == "__main__":
    m, results, diagnostics = main(
        ro_props="seawater",
        ro_dimension="1d",
        erd_config=ERDtype.no_ERD,
        uvdimension=uvdimension.zero_d,
        has_aop=True,
        diagnostics_active=True,
    )
