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

import yaml

from pyomo.environ import units as pyunits

__author__ = "Alexander V. Dudchenko"


def get_source_water_data(file_location, use_watertap_convention=True):
    # simple function to load feed water compostion from yaml file
    with open(file_location, "r") as ymlfile:
        data_dict = yaml.safe_load(ymlfile)
    # print(data_dict)
    if use_watertap_convention:
        data_dict = key_convert(data_dict)
    default_dict = {}
    default_dict["solute_list"] = get_solute_dict(data_dict)
    default_dict["diffusivity_data"] = gen_diffusivity_dict(data_dict)
    default_dict["mw_data"] = gen_mw_dict(data_dict)
    default_dict["stokes_radius_data"] = gen_stoke_dict(data_dict)
    default_dict["charge"] = gen_charge_dict(data_dict)
    tracked_solids = data_dict.get("tracked_solids")
    mass_comp_dict, pH = get_mole_comp(data_dict)
    # print(mass_comp_dict)
    return default_dict, mass_comp_dict, pH, tracked_solids


def key_convert(data_dict):
    new_dict = {}
    new_dict["solvent_list"] = {}
    new_dict["solute_list"] = {}
    for solvent in data_dict["solvent_list"].keys():
        new_dict["solvent_list"][solvent] = data_dict["solvent_list"][solvent]

    for solute in data_dict["solute_list"].keys():
        sol_charge = int(data_dict["solute_list"][solute]["elemental charge"])
        print(sol_charge)
        if sol_charge == 1:
            charge = "+"
        if sol_charge > 1:
            charge = "+{}".format(sol_charge)
        if sol_charge == -1:
            charge = "-"
        if sol_charge < -1:
            charge = "-{}".format(sol_charge)
        if sol_charge == 0:
            charge = ""
        new_dict["solute_list"]["{}_{}".format(solute, charge)] = data_dict[
            "solute_list"
        ][solute]
    return new_dict


def get_solute_dict(data_dict):
    solute_list = list(data_dict["solute_list"].keys())
    return solute_list


def gen_diffusivity_dict(data_dict):
    diff_dict = {}
    for solute in data_dict["solute_list"].keys():
        diff_dict[("Liq", solute)] = float(
            data_dict["solute_list"][solute]["diffusivity"]
        )
    return diff_dict


def gen_mw_dict(data_dict):
    mw_dict = {}
    for solute in data_dict["solvent_list"].keys():
        # print(data_dict["solvent_list"][solute])
        mw_dict[solute] = float(
            data_dict["solvent_list"][solute]["molecular_weight (kg/mol)"]
        )
    for solute in data_dict["solute_list"].keys():
        mw_dict[solute] = float(
            data_dict["solute_list"][solute]["molecular_weight (kg/mol)"]
        )
    return mw_dict


def gen_stoke_dict(data_dict):
    stokes_dict = {}
    for solute in data_dict["solute_list"].keys():
        stokes_dict[solute] = float(
            data_dict["solute_list"][solute]["stokes_radius (m)"]
        )
    return stokes_dict


def gen_charge_dict(data_dict):
    charge_dict = {}
    for solute in data_dict["solute_list"].keys():
        charge_dict[solute] = float(
            data_dict["solute_list"][solute]["elemental charge"]
        )
    return charge_dict


def get_mole_comp(data_dict):
    mass_loading_dict = {}
    for solute in data_dict["solute_list"].keys():
        if "concetration (kg/kg)" in data_dict["solute_list"][solute]:
            mass_loading_dict[solute] = float(
                data_dict["solute_list"][solute]["concetration (kg/kg)"]
            )
        elif "concetration (mg/kg)" in data_dict["solute_list"][solute]:
            mass_loading_dict[solute] = (
                float(data_dict["solute_list"][solute]["concetration (mg/kg)"]) / 1e6
            )
    print(data_dict)
    return mass_loading_dict, float(data_dict["pH"])


def set_feed(blk, feed_mass_frac, flow_basis=1):
    mass_flow_in = flow_basis * pyunits.kg / pyunits.s
    H2O_mass_frac = 1 - sum(x for x in feed_mass_frac.values())
    H2O_mass_comp_flow = H2O_mass_frac * pyunits.kg / pyunits.kg * mass_flow_in
    # blk.feed.properties[0].flow_vol_phase["Liq"].fix(mass_flow_in / 1000)
    blk.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(mass_flow_in)
    for ion, x in feed_mass_frac.items():
        mass_comp_flow = x * pyunits.kg / pyunits.kg * mass_flow_in
        blk.feed.properties[0].flow_mass_phase_comp["Liq", ion].fix(mass_comp_flow)
    if ("Liq", "Cl") in blk.feed.properties[0].flow_mass_phase_comp.keys():
        blk.feed.properties[0].assert_electroneutrality(
            defined_state=True,
            adjust_by_ion="Cl",
            get_property="flow_mass_phase_comp",
        )
    blk.feed.properties[0].temperature.fix(298.15)
    blk.feed.properties[0].pressure.fix(101325)
