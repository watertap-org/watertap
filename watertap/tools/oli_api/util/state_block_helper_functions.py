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
#
#################################################################################

__author__ = "Paul Vecchiarelli, Adam Atia"

from pyomo.environ import (
    units as pyunits,
    value,
    ConcreteModel,
    assert_optimal_termination,
)

from idaes.core import FlowsheetBlock, MaterialFlowBasis

from idaes.core.util.scaling import calculate_scaling_factors

from idaes.core.solvers import get_solver

from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock

from watertap.tools.oli_api.util.watertap_to_oli_helper_functions import (
    get_oli_names,
    oli_reverse_lookup,
    names_db,
)


def create_state_block(source_water):
    """
    Creates a state block using the Multi Component Aqueous Solution property model

    :param source_water: dictionary containing state variables and units

    :returns m: ConcreteModel containing MCAS state block
    """

    solver = get_solver()

    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    props = create_property_model_input(
        source_water["components"], property_model_type="mcas"
    )

    m.fs.properties = MCASParameterBlock(
        **props, material_flow_basis=MaterialFlowBasis.mass
    )

    stream = m.fs.stream = m.fs.properties.build_state_block([0])

    for comp, conc in stream[0].flow_mass_phase_comp.items():

        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1 / conc.value, index=comp
        )

    for comp, conc in stream[0].conc_mass_phase_comp.items():

        m.fs.properties.set_default_scaling(
            "conc_mass_phase_comp", 1 / conc.value, index=comp
        )

    m.fs.properties.set_default_scaling(
        "flow_vol_phase", 1 / stream[0].flow_vol_phase["Liq"].value, index=("Liq")
    )

    calculate_scaling_factors(m)

    # setting state variables
    convert_to_state_block_units(stream[0].temperature, source_water, "temperature")

    convert_to_state_block_units(stream[0].pressure, source_water, "pressure")

    stream[0].conc_mass_phase_comp

    convert_conc = lambda conc: pyunits.convert_value(
        conc,
        from_units=source_water["units"]["components"],
        to_units=stream[0].conc_mass_phase_comp._units,
    )

    conc_basis = "conc_mass_phase_comp"

    phase = "Liq"

    var_args = {}

    for comp in source_water["components"]:

        var_args[("conc_mass_phase_comp", ("Liq", comp.oli_name))] = convert_conc(
            source_water["components"][comp]
        )

    var_args.update({("flow_vol_phase", "Liq"): 1e-3})

    stream.calculate_state(var_args=var_args, hold_state=True)

    stream.initialize()

    result = solver.solve(m)

    assert_optimal_termination(result)

    return m


def convert_to_state_block_units(state_variable, source: dict, key):

    converted_value = pyunits.convert_value(
        source[key], from_units=source["units"][key], to_units=state_variable._units
    )

    state_variable.fix(converted_value)


def extract_state_vars(state_block, conc_var, units):

    update_conc = lambda conc, var, key: pyunits.convert_value(
        conc, var._units, units[key]
    )

    components = {}

    for comp, conc in conc_var.items():

        if comp[1] != "H2O":

            oli_name = oli_reverse_lookup(comp[1], names_db)

            components[oli_name] = update_conc(conc.value, conc_var, "components")

    update_t_p = lambda var, key: pyunits.convert_value(
        var.value, var._units, units[key]
    )

    temp = update_t_p(state_block.temperature, "temperature")

    pressure = update_t_p(state_block.pressure, "pressure")

    state_vars = {
        "components": components,
        "temperature": temp,
        "pressure": pressure,
    }

    return state_vars


def create_property_model_input(components, property_model_type: str = ""):
    """
    Builds property package inputs.

    :param solutes: dictionary containing solute concentrations in mg/L.
    """

    if property_model_type == "mcas":

        try:

            property_model_input = {
                "solute_list": [component.oli_name for component in components],
                "mw_data": {
                    component.oli_name: component.molar_mass * 1e-3
                    for component in components
                },
                "charge": {
                    component.oli_name: component.charge
                    for component in components
                    if component.charge != 0
                },
            }

            return property_model_input

        except AttributeError:

            raise RuntimeError(
                " Run get_oli_names to fetch required attributes for create_property_model_input."
            )
