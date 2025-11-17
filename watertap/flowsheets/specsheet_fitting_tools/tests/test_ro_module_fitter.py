import pytest

from watertap.flowsheets.specsheet_fitting_tools.ro_module_fitter import (
    fit_bw30_4040_to_spec_sheet,
    fit_sw30_4040_to_spec_sheet,
    fit_ro_module_to_spec_sheet,
)
import os

_test_file_path = os.path.dirname(os.path.abspath(__file__))

import yaml


def load_test_yaml(file_name, remove_file=False):
    with open(os.path.join(_test_file_path, file_name), "r") as f:
        data = yaml.safe_load(f)
    if remove_file:
        os.remove(file_name)
    return data


def validate_dict(dict1, dict2, rel_tol=1e-2):
    for key in dict1:
        assert pytest.approx(dict1[key]["value"], rel=rel_tol) == dict2[key]["value"]
        assert dict1[key]["units"] == dict2[key]["units"]


def validate_result(ro_module, result):
    expected_data = load_test_yaml(
        f"{_test_file_path}/expected_ro_results/{ro_module}.yaml"
    )
    validate_dict(result, expected_data)
    load_result = load_test_yaml(
        f"{_test_file_path}/{ro_module}.yaml", remove_file=True
    )
    validate_dict(load_result, expected_data)


def test_fit_bw30_4040_to_spec_sheet():
    result = fit_bw30_4040_to_spec_sheet(save_location=_test_file_path)
    validate_result("BW30 PRO-4040", result)


def test_fit_sw30_4040_to_spec_sheet():
    result = fit_sw30_4040_to_spec_sheet(save_location=_test_file_path)
    validate_result("SW30-4040", result)


def test_fit_ro_module_to_spec_sheet():
    result = fit_ro_module_to_spec_sheet(save_location=_test_file_path)
    validate_result("RO_module", result)
