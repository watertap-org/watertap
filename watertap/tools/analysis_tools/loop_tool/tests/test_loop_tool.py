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
import pytest
import os
import numpy as np

from watertap.tools.analysis_tools.loop_tool.tests import ro_setup
from watertap.tools.analysis_tools.loop_tool.loop_tool import loopTool, get_working_dir
import yaml
import h5py

_this_file_path = os.path.dirname(os.path.abspath(__file__))


@pytest.fixture()
def loop_test_options_setup():
    cwd = get_working_dir()
    lp = loopTool(
        _this_file_path + "/test_all_options.yaml",
        build_function=ro_setup.ro_build,
        initialize_function=ro_setup.ro_init,
        optimize_function=ro_setup.ro_solve,
        saving_dir=_this_file_path,
        save_name="ro_with_erd",
        execute_simulations=False,
        number_of_subprocesses=1,
    )
    lp.build_run_dict()
    """ used to generate test file"""
    # with open("test_expected_option_directory.yaml", "w") as file:
    #     documents = yaml.dump(lp.sweep_directory, file)
    with open(_this_file_path + "/test_expected_option_directory.yaml", "r") as infile:
        expected_run_dict = yaml.safe_load(infile)

    return lp, expected_run_dict


@pytest.fixture()
def loop_sweep_setup():
    cwd = get_working_dir()
    lp = loopTool(
        _this_file_path + "/test_sweep.yaml",
        build_function=ro_setup.ro_build,
        initialize_function=ro_setup.ro_init,
        optimize_function=ro_setup.ro_solve,
        saving_dir=_this_file_path,
        save_name="ro_with_erd",
        execute_simulations=False,
        number_of_subprocesses=1,
    )
    lp.build_run_dict()
    """ used to generate test file"""
    # with open("test_expected_sweep_directory.yaml", "w") as file:
    #     documents = yaml.dump(lp.sweep_directory, file)
    with open(_this_file_path + "/test_expected_sweep_directory.yaml", "r") as infile:
        expected_run_dict = yaml.safe_load(infile)

    return lp, expected_run_dict


@pytest.fixture()
def loop_diff_setup():
    # Test without parallel implementation as its broken in
    # water tap paramtersweep tool
    cdw = get_working_dir()
    lp = loopTool(
        _this_file_path + "/test_diff.yaml",
        build_function=ro_setup.ro_build,
        initialize_function=ro_setup.ro_init,
        optimize_function=ro_setup.ro_solve,
        saving_dir=_this_file_path,
        save_name="ro_with_erd",
        execute_simulations=False,
        number_of_subprocesses=1,
    )
    lp.build_run_dict()
    """ used to generate test file"""
    # with open("test_expected_diff_directory.yaml", "w") as file:
    #     documents = yaml.dump(lp.sweep_directory, file)
    with open(_this_file_path + "/test_expected_diff_directory.yaml", "r") as infile:
        expected_run_dict = yaml.safe_load(infile)

    return lp, expected_run_dict


@pytest.mark.component
def test_options_setups(loop_test_options_setup):
    lp, expected_run_dict = loop_test_options_setup
    lp.build_run_dict()

    def test_diff_dict(dicta, dictb):
        for key in dicta:
            if key != "dir":
                if isinstance(dicta[key], dict):
                    test_diff_dict(dicta[key], dictb[key])

                elif dicta[key] != dictb[key]:
                    # print(dicta[key], dictb[key])
                    return False
                # else:
                #   break

        return True

    assert test_diff_dict(lp.sweep_directory, expected_run_dict)


@pytest.mark.component
def test_sweep_setup(loop_sweep_setup):
    lp, expected_run_dict = loop_sweep_setup
    lp.build_run_dict()

    def test_diff_dict(dicta, dictb):
        for key in dicta:
            if key != "dir":
                if isinstance(dicta[key], dict):
                    test_diff_dict(dicta[key], dictb[key])

                elif dicta[key] != dictb[key]:
                    # print(dicta[key], dictb[key])
                    return False
                # else:
                #   break

        return True

    assert test_diff_dict(lp.sweep_directory, expected_run_dict)


@pytest.mark.component
def test_diff_setup(loop_diff_setup):
    lp, expected_run_dict = loop_diff_setup
    lp.build_run_dict()

    def test_diff_dict(dicta, dictb):
        for key in dicta:
            if key != "dir":
                if isinstance(dicta[key], dict):
                    test_diff_dict(dicta[key], dictb[key])

                elif dicta[key] != dictb[key]:
                    # print(dicta[key], dictb[key])
                    return False
                # else:
                #   break

        return True

    assert test_diff_dict(lp.sweep_directory, expected_run_dict)


@pytest.mark.component
def test_sweep_run(loop_sweep_setup):
    lp, test_file = loop_sweep_setup
    lp.build_run_dict()
    # remove any existing file before test
    if os.path.isfile(lp.h5_file_location_default + "_analysisType_ro_analysis.h5"):
        os.remove(lp.h5_file_location_default + "_analysisType_ro_analysis.h5")

    lp.run_simulations()

    h5file = h5py.File(
        lp.h5_file_location_default + "_analysisType_ro_analysis.h5", "r"
    )
    data = h5file[
        "ro_analysis/erd_type/pressure_exchanger/membrane_cost/outputs/fs.costing.LCOW/value"
    ]

    true_vals = [0.37203417, 0.39167574, 0.41117995]
    d = data[()]
    # print(true_vals, d)
    for i, tv in enumerate(true_vals):
        assert d[i] == pytest.approx(tv, rel=1e-2)
    data = h5file[
        "ro_analysis/erd_type/pump_as_turbine/membrane_cost/outputs/fs.costing.LCOW/value"
    ]

    true_vals = [0.50886109, 0.52850266, 0.54814424]
    d = data[()]
    # print(true_vals, d)
    for i, tv in enumerate(true_vals):
        assert d[i] == pytest.approx(tv, rel=1e-2)
    data = h5file[
        "ro_analysis/erd_type/pressure_exchanger/membrane_group/outputs/fs.costing.LCOW/value"
    ]

    true_vals = [
        0.3810009713006634,
        0.3916757385992817,
        0.3985075912766517,
        0.4111799488092862,
    ]
    d = data[()]
    # print(true_vals, d)
    for i, tv in enumerate(true_vals):
        assert d[i] == pytest.approx(tv, rel=1e-2)
    data = h5file[
        "ro_analysis/erd_type/pump_as_turbine/membrane_group/outputs/fs.costing.LCOW/value"
    ]

    true_vals = [
        0.5178278972844322,
        0.5285026645830471,
        0.5353345156541605,
        0.5481442364124981,
    ]
    d = data[()]
    # print(true_vals, d)
    for i, tv in enumerate(true_vals):
        assert d[i] == pytest.approx(tv, rel=1e-2)
    h5file.close()
    """test that backup works, will not run actual simulation ,create a back up file, and
    load data from it into sim file. The lp.back_file_name should not be None
    """


@pytest.mark.component
def test_sweep_backup(loop_sweep_setup):
    """test that backup works, will not run actual simulation ,create a back up file, and
    load data from it into sim file. The lp.back_file_name should not be None
    """
    lp, test_file = loop_sweep_setup
    lp.build_run_dict()
    lp.run_simulations()

    assert lp.h5_backup_location != None

    h5file = h5py.File(
        lp.h5_file_location_default + "_analysisType_ro_analysis.h5", "r"
    )
    data = h5file[
        "ro_analysis/erd_type/pressure_exchanger/membrane_cost/outputs/fs.costing.LCOW/value"
    ]

    true_vals = [0.37203417, 0.39167574, 0.41117995]
    d = data[()]
    for i, tv in enumerate(true_vals):
        assert d[i] == pytest.approx(tv, rel=1e-2)
    data = h5file[
        "ro_analysis/erd_type/pump_as_turbine/membrane_cost/outputs/fs.costing.LCOW/value"
    ]
    true_vals = [0.50886109, 0.52850266, 0.54814424]
    d = data[()]

    for i, tv in enumerate(true_vals):
        assert d[i] == pytest.approx(tv, rel=1e-2)
    h5file.close()

    os.remove(lp.h5_file_location_default + "_analysisType_ro_analysis.h5")
    # try:
    os.remove(lp.h5_backup_location)


@pytest.mark.component
def test_diff_run(loop_diff_setup):
    lp, test_file = loop_diff_setup
    lp.build_run_dict()
    lp.run_simulations()
    h5file = h5py.File(
        lp.h5_file_location_default + "_analysisType_ro_diff_analysis.h5", "r"
    )

    data = h5file["ro_diff_analysis/membrane_cost/outputs/fs.costing.LCOW/value"][()]

    # for i, tv in enumerate(true_vals):
    assert len(data) == 4

    data_a = h5file[
        "ro_diff_analysis/membrane_group/sweep_params/fs.costing.reverse_osmosis.factor_membrane_replacement/value"
    ][()]
    data_b = h5file[
        "ro_diff_analysis/membrane_group/sweep_params/fs.costing.reverse_osmosis.membrane_cost/value"
    ][()]
    # for i, tv in enumerate(true_vals):
    assert len(data_a) == 4
    assert len(data_b) == 4
    h5file.close()
    # try:
    os.remove(lp.h5_file_location_default + "_analysisType_ro_diff_analysis.h5")
    # except OSError:
    #     pass
