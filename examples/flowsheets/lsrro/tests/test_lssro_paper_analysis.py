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

# TODO: If needed for coverage, *build* and set operating conditions for,
#       but do not *solve*, the cross of all the options.

import glob
import math
import pytest
import os

import pandas as pd

from pyomo.environ import check_optimal_termination, value

from watertap.examples.flowsheets.lsrro.lsrro import (
    ACase,
    BCase,
    ABTradeoff,
    run_lsrro_case,
)

from .gha_divider import get_test_cases_subset

_this_file_path = os.path.dirname(os.path.abspath(__file__))

_input_headers = {
    "cin (kg/m3)": (float, "Cin"),
    "recovery (-)": (float, "water_recovery"),
    "num_stages": (int, "number_of_stages"),
    "A_case": (ACase, "A_case"),
    "B_case": (BCase, "B_case"),
    "AB_Tradeoff": (ABTradeoff, "AB_tradeoff"),
    "AB_gamma_factor": (float, "AB_gamma_factor"),
}

_results_headers = {
    "final brine concentration": "fs.disposal.properties[0].conc_mass_phase_comp[Liq, NaCl]",
    "final perm (ppm)": "fs.final_permeate_concentration",
    "Membrane area": "fs.total_membrane_area",
    "SEC": "fs.costing.specific_energy_consumption",
    "LCOW": "fs.costing.LCOW",
    "LCOW_feed": "fs.costing.LCOW_feed",
    "primary_pump_capex": "fs.costing.primary_pump_capex_lcow",
    "booster_pump_capex": "fs.costing.booster_pump_capex_lcow",
    "erd_capex": "fs.costing.erd_capex_lcow",
    "membrane_capex": "fs.costing.membrane_capex_lcow",
    "indirect_capex": "fs.costing.indirect_capex_lcow",
    "electricity": "fs.costing.electricity_lcow",
    "membrane_replacement": "fs.costing.membrane_replacement_lcow",
    "chem_lab_main": "fs.costing.chemical_labor_maintenance_lcow",
    "pumping_energy_agg_costs": "fs.costing.pumping_energy_aggregate_lcow",
    "membrane_agg_costs": "fs.costing.membrane_aggregate_lcow",
}

_csv_files = sorted(
    glob.glob(os.path.join(_this_file_path, "paper_analysis_baselines", "*.csv"))
)

_dfs = {os.path.basename(csv_file): pd.read_csv(csv_file) for csv_file in _csv_files}

_test_cases = [
    (csv_file, idx) for csv_file, df in _dfs.items() for idx in range(len(df))
]

# comment out this line if you want to run the entire baseline
_test_cases = get_test_cases_subset(_test_cases)


@pytest.mark.parametrize("csv_file, row_index", _test_cases)
@pytest.mark.component
def test_against_paper_analysis(csv_file, row_index):

    row = _dfs[csv_file].iloc[row_index]
    input_arguments = {
        argument: converter(row[property_name])
        for property_name, (converter, argument) in _input_headers.items()
    }
    number_of_stages = input_arguments["number_of_stages"]
    model, results = run_lsrro_case(
        **input_arguments,
        has_NaCl_solubility_limit=True,
        permeate_quality_limit=1000e-6,
        has_calculated_concentration_polarization=True,
        has_calculated_ro_pressure_drop=True,
        A_value=5 / 3.6e11,
        B_max=None,
        number_of_RO_finite_elements=10
    )

    if check_optimal_termination(results):
        for property_name, flowsheet_attribute in _results_headers.items():
            assert value(model.find_component(flowsheet_attribute)) == pytest.approx(
                float(row[property_name]),
                rel=1e-3,
            )
    else:
        for property_name in _results_headers:
            assert math.isnan(float(row[property_name]))
