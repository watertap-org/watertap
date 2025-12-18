from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    value,
    Var,
    NonNegativeReals,
    assert_optimal_termination,
    Objective,
    units as pyunits,
)

from pyomo.util.calc_var_value import calculate_variable_from_constraint

from watertap.core.solvers import get_solver

import watertap.flowsheets.ccro.utils.utils as cc_utils
from watertap.flowsheets.ccro.ccro_flowsheet_functions import scaling as ccro_scaling


def set_operating_conditions(m, cc_configuration=None, **kwargs):

    if m.fs.operation_mode == "filtration":

        # Feed block operating conditions
        m.fs.raw_feed.properties[0].pressure.fix(101325 * pyunits.Pa)
        m.fs.raw_feed.properties[0].pressure_osm_phase["Liq"]
        m.fs.raw_feed.properties[0].flow_vol_phase["Liq"]
        m.fs.raw_feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"]
        m.fs.raw_feed.properties[0].temperature.fix(cc_configuration["temperature"])
        m.fs.raw_feed.properties[0].flow_vol_phase["Liq"].fix(
            cc_configuration["raw_feed_flowrate"]
        )
        m.fs.raw_feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].fix(
            cc_configuration["raw_feed_conc"]
        )

        solver = get_solver()
        solver.solve(m.fs.raw_feed)
        # Dead Volume operating conditions
        # Fix volume
        if cc_configuration["dead_volume"] == "base_on_dead_volume_to_area_ratio":
            dead_volume = (
                cc_configuration["dead_volume_to_area_ratio"]
                * cc_configuration["membrane_area"]
            )
            print(f"Dead Volume set to {dead_volume}")
        else:
            dead_volume = cc_configuration["dead_volume"] = cc_configuration[
                "dead_volume"
            ]
        if m.fs.ro_model_with_hold_up:
            dead_volume = dead_volume * cc_configuration["pipe_to_module_ratio"]
        m.fs.dead_volume.volume.fix(dead_volume)
        m.fs.dead_volume.delta_state.volume.fix(dead_volume)

        m.fs.dead_volume.accumulation_time.fix(cc_configuration["accumulation_time"])

        # Fixing the flow rate of the dead volume delta state to the feed to calculate the mass fraction and density
        # I found fixing mass fraction and density is easiest way to get initial state
        # we will also use these as connection points between current and future state.

        m.fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"].fix(
            m.fs.raw_feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"].value
        )
        m.fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"].fix(
            m.fs.raw_feed.properties[0].dens_mass_phase["Liq"].value
        )
        if m.fs.ro_model_with_hold_up:
            m.fs.RO.feed_side.volume.fix(
                cc_configuration["dead_volume_to_area_ratio"]
                * cc_configuration["membrane_area"]
            )
            m.fs.RO.feed_side.accumulation_time.fix(
                cc_configuration["accumulation_time"]
            )
            idx = m.fs.RO.difference_elements
            for i in idx:
                m.fs.RO.feed_side.delta_state.node_mass_frac_phase_comp[
                    0, i, "Liq", "NaCl"
                ].fix(
                    m.fs.raw_feed.properties[0]
                    .mass_frac_phase_comp["Liq", "NaCl"]
                    .value
                )
                m.fs.RO.feed_side.delta_state.node_dens_mass_phase[0, i, "Liq"].fix(
                    m.fs.raw_feed.properties[0].dens_mass_phase["Liq"].value
                )
            solver.solve(m.fs.RO.feed_side.delta_state)
        # Pump 1 operating conditions
        m.fs.P1.efficiency_pump.fix(cc_configuration["p1_eff"])

        cc_utils.set_pump_operating_pressure(
            unit=m.fs.P1,
            configuration_options=cc_configuration,
            osmotic_pressure=m.fs.raw_feed.properties[0].pressure_osm_phase["Liq"],
        )
        cc_utils.set_pump_operating_pressure(
            unit=m.fs.P2,
            configuration_options=cc_configuration,
            osmotic_pressure=m.fs.raw_feed.properties[0].pressure_osm_phase["Liq"],
        )
        # Pump 2 operating conditions
        m.fs.P2.efficiency_pump.fix(cc_configuration["p2_eff"])

        # Set RO configuration parameters
        m.fs.RO.A_comp.fix(cc_configuration["A_comp"])
        m.fs.RO.B_comp.fix(cc_configuration["B_comp"])
        m.fs.RO.area.fix(cc_configuration["membrane_area"])
        m.fs.RO.length.unfix()
        m.fs.RO.width.unfix()
        m.fs.RO.feed_side.velocity[0, 0].fix(0.15)
        m.fs.RO.feed_side.channel_height.fix(cc_configuration["channel_height"])
        m.fs.RO.feed_side.spacer_porosity.fix(cc_configuration["spacer_porosity"])

        m.fs.RO.permeate.pressure[0].fix(101325 * pyunits.Pa)

        m.fs.RO.feed_side.friction_factor_darcy.setub(200)
        # m.fs.RO.flux_mass_phase_comp.setub(1)
        # m.fs.RO.feed_side.cp_modulus.setub(50)
        # m.fs.RO.feed_side.cp_modulus.setlb(0.1)
        m.fs.RO.deltaP.setlb(None)

        ccro_scaling.scale_filtration_system(m)

    elif m.fs.operation_mode == "flushing":

        # Pump 1 operating conditions
        m.fs.P1.efficiency_pump.fix(cc_configuration["p1_eff"])

        cc_utils.set_pump_operating_pressure(
            unit=m.fs.P1,
            configuration_options=cc_configuration,
            osmotic_pressure=m.fs.raw_feed.properties[0].pressure_osm_phase["Liq"],
        )

        # Pump 2 operating conditions - Add only for costing function. No work is done by this pump
        m.fs.P2.efficiency_pump.fix(cc_configuration["p2_eff"])

        # Concentration of the flushing water is the raw feed concentration
        if cc_configuration["flushing_conc"] == "raw_feed_conc":
            m.fs.flushing.flushing_feed_concentration.fix(
                cc_configuration["raw_feed_conc"]
            )
        else:
            m.fs.flushing.flushing_feed_concentration.fix(
                cc_configuration["flushing_conc"]
            )

        # Calculating values for the dead volume and delta state
        # Feed block operating conditions
        m.fs.raw_feed.properties[0].pressure.fix(101325 * pyunits.Pa)
        m.fs.raw_feed.properties[0].temperature.fix(cc_configuration["temperature"])

        m.fs.raw_feed.properties[0].pressure_osm_phase["Liq"]
        m.fs.raw_feed.properties[0].flow_vol_phase["Liq"]
        m.fs.raw_feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"]

        # Reassign the raw feed flowrate and concentration
        if cc_configuration["flushing_flowrate"] == "raw_feed_flowrate":
            m.fs.raw_feed.properties[0].flow_vol_phase["Liq"].fix(
                cc_configuration["raw_feed_flowrate"]
            )
        else:
            m.fs.raw_feed.properties[0].flow_vol_phase["Liq"].fix(
                cc_configuration["flushing_flowrate"]
            )
        if cc_configuration["flushing_conc"] == "raw_feed_conc":
            m.fs.raw_feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].fix(
                cc_configuration["raw_feed_conc"]
            )
        else:
            m.fs.raw_feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"].fix(
                cc_configuration["flushing_conc"]
            )
        solver = get_solver()
        solver.solve(m.fs.raw_feed)

        # Dead Volume operating conditions
        # Fix volume
        if cc_configuration["dead_volume"] == "base_on_dead_volume_to_area_ratio":
            dead_volume = (
                cc_configuration["dead_volume_to_area_ratio"]
                * cc_configuration["membrane_area"]
            )
        else:
            dead_volume = cc_configuration["dead_volume"] = cc_configuration[
                "dead_volume"
            ]
        if m.fs.ro_model_with_hold_up:
            dead_volume = dead_volume * (1 + cc_configuration["pipe_to_module_ratio"])
        m.fs.dead_volume.volume.fix(dead_volume)
        m.fs.dead_volume.delta_state.volume.fix(dead_volume)
        m.fs.dead_volume.accumulation_time.fix(cc_configuration["accumulation_time"])
        m.fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"].fix()
        m.fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"].fix()
        m.fs.flushing.flushing_efficiency.fix(cc_configuration["flushing_efficiency"])

        # Constraints
        @m.fs.Constraint()
        def dead_volume_residence_time_constraint(m):
            return m.flushing.mean_residence_time == (
                pyunits.convert(
                    m.dead_volume.volume[0, "Liq"]
                    / m.raw_feed.properties[0].flow_vol_phase["Liq"],
                    to_units=pyunits.s,
                )
            )

        # Calculate pre-flushing/ dead volume delta state concentration. Concentration before
        # flushing should be the delta state concentration
        @m.fs.Constraint()
        def pre_flushing_conc_constraint(m):
            return (
                m.flushing.pre_flushing_concentration
                == m.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"]
                * m.dead_volume.delta_state.dens_mass_phase[0, "Liq"]
            )

        # Concentration after flushing should be the dead volume properties out concentration
        @m.fs.Constraint()
        def post_flushing_conc_constraint(m):
            return (
                m.dead_volume.dead_volume.properties_out[0].conc_mass_phase_comp[
                    "Liq", "NaCl"
                ]
                == m.flushing.post_flushing_concentration
            )

        ccro_scaling.scale_flushing_system(m)
    return m


def unfix_dof(m, unfix_dead_volume_state=True, cc_configuration=None, **kwargs):
    """
    Unfix linking variables in MP model
    """

    if unfix_dead_volume_state:
        m.fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"].unfix()
        m.fs.dead_volume.delta_state.dens_mass_phase[0, "Liq"].unfix()

    if m.fs.operation_mode == "filtration":

        m.fs.P1.control_volume.properties_out[0].pressure.unfix()
        m.fs.P2.control_volume.properties_out[0].pressure.unfix()

        m.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"].fix(
            cc_configuration["recycle_flowrate"]
        )
        if m.fs.ro_model_with_hold_up:
            idx = m.fs.RO.difference_elements
            for i in idx:
                m.fs.RO.feed_side.delta_state.node_mass_frac_phase_comp[
                    0, i, "Liq", "NaCl"
                ].unfix()
                m.fs.RO.feed_side.delta_state.node_dens_mass_phase[0, i, "Liq"].unfix()
    elif m.fs.operation_mode == "flushing":
        m.fs.dead_volume.accumulation_time.unfix()
        m.fs.flushing.flushing_efficiency.fix(cc_configuration["flushing_efficiency"])
        m.fs.flushing.flushing_time.unfix()
        m.fs.flushing.mean_residence_time.unfix()
        m.fs.dead_volume.dead_volume.mass_frac_phase_comp[0, "Liq", "NaCl"].unfix()
        m.fs.flushing.pre_flushing_concentration.unfix()
        m.fs.flushing.post_flushing_concentration.unfix()
