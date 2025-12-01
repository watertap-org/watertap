from pyomo.environ import (
    value,
    units as pyunits,
)
from watertap.flowsheets.ccro.multiperiod.model_state_tool import ModelState
from idaes.core import UnitModelCostingBlock


def copy_time_period_links(m_old, m_new, vars_to_copy=None):
    """
    Copy linking variables between time periods.
    Args:
        m_old: model of the old time period
        m_new: model of the new time period
        vars_to_copy: list of dictionaries with old and new model variable names to copy for example:
            {"old_model_var": 'fs.dead_volume.dead_volume.properties_out[0].mass_frac_phase_comp["Liq", "NaCl"]',
             "new_model_var": 'fs.dead_volume.delta_state.mass_frac_phase_comp[0, "Liq", "NaCl"]'}
    """

    if vars_to_copy is not None:
        for var in vars_to_copy:
            assert (
                m_new.find_component(var["new_model_var"]) is not None
            ), f"{var} not found in new model"
            assert (
                m_new.find_component(var["old_model_var"]) is not None
            ), f"{var} not found in new model"
            m_new.find_component(var["new_model_var"]).fix(
                m_old.find_component(var["old_model_var"]).value
            )


def copy_state(old_model, new_model):
    model_state = ModelState()
    model_state.get_model_state(old_model)
    model_state.set_model_state(new_model)


def register_costed_unit(
    mp, unit, costing_method_arguments={}, register_electricity_flow_only=False
):
    if register_electricity_flow_only:
        lb = unit.work_mechanical[0.0].lb
        # set lower bound to 0 to avoid negative defined flow warning when lb is not >= 0
        unit.work_mechanical.setlb(0)
        mp.costing.cost_flow(
            pyunits.convert(unit.work_mechanical[0.0], to_units=pyunits.kW),
            "electricity",
        )
        # set lower bound back to its original value that was assigned to lb
        unit.work_mechanical.setlb(lb)
    else:
        unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=mp.costing,
            costing_method_arguments=costing_method_arguments,
        )


def set_pump_operating_pressure(unit, osmotic_pressure, configuration_options):
    """
    Set pump operating pressure based on configuration option.
    Args:
        unit: pump unit model
        osm_pressure: osmotic pressure at the pump inlet
        configuration_options: CCROConfiguration object
    """
    if configuration_options["p1_pressure_start"] == "osmotic_pressure":
        set_pressure = (
            value(osmotic_pressure)
            * value(configuration_options["osmotic_overpressure"])
            + 2e5
        )  # Pa
        print(
            f"Setting P1 outlet pressure based on osmotic pressure of: {value(osmotic_pressure)/1e5} bar, pressure is set to {(set_pressure)/1e5} bar"
        )
        unit.control_volume.properties_out[0].pressure.fix(set_pressure)
    else:
        unit.control_volume.properties_out[0].pressure.fix(
            configuration_options["p1_pressure_start"]
        )
