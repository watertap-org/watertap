###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

# TODO: Load and spot-check base cases corresponding to 2, 3, and 4 stage systems
# for 35g/L, 70g/L, and 125 g/L feed concentrations as base-line coverage for LSRRO flowsheet

# TODO: Maybe we want to randomly generate which of the various test
#       scenarios we run, so we get coverage in expectation.

# TODO: If needed for coverage, *build* and set operating conditions for,
#       but do not *solve*, the cross of all the options.

import pytest

from watertap.examples.flowsheets.lsrro.lsrro import ACase, BCase, ABTradeoff, run_lsrro_case


_input_headers = {
    "cin (kg/m3)": float,
    "recovery (-)": float,
    "num_stages": int,
    "A_case": ACase,
    "B_case": BCase,
    "AB_Tradeoff": ABTradeoff,
    "AB_gamma_factor": float,
}

_results_headers = {
    "final brine concentration": "fs.disposal.properties[0].conc_mass_phase_comp[Liq, NaCl]",
    "final perm (ppm)": "fs.final_permeate_concentration",
    "Membrane area": "fs.total_membrane_area",
    "SEC": "fs.costing.specific_energy_consumption",
    "LCOW": "fs.costing.LCOW",
    "LCOW_feed": "fs.costing.LCOW_feed",
    "primary_pump_capex":
        "fs.costing.primary_pump_capex_lcow"
    ,
    "booster_pump_capex":
        "fs.costing.booster_pump_capex_lcow"
    ,
    "erd_capex": 
        "fs.costing.erd_capex_lcow"
    ,
    "membrane_capex":
        "fs.costing.membrane_capex_lcow"
    ,
    "indirect_capex":
        "fs.costing.indirect_capex_lcow"
    ,
    "electricity":
        "fs.costing.electricity_lcow"
    ,
    "membrane_replacement": 
        "fs.costing.membrane_replacement_lcow"
    ,
    "chem_lab_main": 
        "fs.costing.chemical_labor_maintenance_lcow",
    # TODO: somewhat number of stages dependent
    #*(f"A_stage {i}" for i in range(1,6)),
    #*(f"B_stage {i}" for i in range(1,6)),
    "pumping_energy_agg_costs":"fs.costing.pumping_energy_aggregate_lcow",
    "membrane_agg_costs":"fs.costing.membrane_aggregate_lcow",
}
