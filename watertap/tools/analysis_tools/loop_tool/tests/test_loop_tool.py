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


@pytest.fixture()
def loop_sweep_setup():
    # Test without parallel implementation as its broken in
    # water tap paramtersweep tool
    lp = loopTool(
        "test_sweep.yaml",
        build_function=ro_setup.ro_build,
        initialize_function=ro_setup.ro_init,
        optimize_function=ro_setup.ro_solve,
        saving_dir=get_working_dir(),
        save_name="ro_with_erd",
        execute_simulations=False,
        number_of_subprocesses=False,
    )
    lp.build_run_dict()
    """ used to generate test file"""
    # with open("test_expected_sweep_directory.yaml", "w") as file:
    #     documents = yaml.dump(lp.sweep_directory, file)
    with open("test_expected_sweep_directory.yaml", "r") as infile:
        expected_run_dict = yaml.safe_load(infile)

    return lp, expected_run_dict


@pytest.fixture()
def loop_diff_setup():
    # Test without parallel implementation as its broken in
    # water tap paramtersweep tool
    cdw = get_working_dir()
    lp = loopTool(
        "test_diff.yaml",
        build_function=ro_setup.ro_build,
        initialize_function=ro_setup.ro_init,
        optimize_function=ro_setup.ro_solve,
        saving_dir=get_working_dir(),
        save_name="ro_with_erd",
        execute_simulations=False,
        number_of_subprocesses=False,
    )
    lp.build_run_dict()
    """ used to generate test file"""
    # with open("test_expected_diff_directory.yaml", "w") as file:
    #     documents = yaml.dump(lp.sweep_directory, file)
    with open("test_expected_diff_directory.yaml", "r") as infile:
        expected_run_dict = yaml.safe_load(infile)

    return lp, expected_run_dict


# def test_sweep_setup(loop_sweep_setup):
#     lp, expected_run_dict = loop_sweep_setup
#     lp.build_run_dict()

#     def test_diff_dict(dicta, dictb):
#         for key in dicta:
#             if key != "dir":
#                 if isinstance(dicta[key], dict):
#                     test_diff_dict(dicta[key], dictb[key])

#                 elif dicta[key] != dictb[key]:
#                     # print(dicta[key], dictb[key])
#                     return False
#                 # else:
#                 #   break

#         return True

#     assert test_diff_dict(lp.sweep_directory, expected_run_dict)


# def test_diff_setup(loop_diff_setup):
#     lp, expected_run_dict = loop_diff_setup
#     lp.build_run_dict()

#     def test_diff_dict(dicta, dictb):
#         for key in dicta:
#             if key != "dir":
#                 if isinstance(dicta[key], dict):
#                     test_diff_dict(dicta[key], dictb[key])

#                 elif dicta[key] != dictb[key]:
#                     # print(dicta[key], dictb[key])
#                     return False
#                 # else:
#                 #   break

#         return True

#     assert test_diff_dict(lp.sweep_directory, expected_run_dict)


# def test_sweep_run(loop_sweep_setup):
#     lp, test_file = loop_sweep_setup
#     lp.build_run_dict()
#     try:
#         os.remove(lp.h5_file_location_default + "_analysisType_ro_analysis.h5")
#     except OSError:
#         pass
#     lp.run_simulations()

#     h5file = h5py.File(
#         lp.h5_file_location_default + "_analysisType_ro_analysis.h5", "r"
#     )
#     data = h5file[
#         "ro_analysis/erd_type/pressure_exchanger/membrane_cost/outputs/fs.costing.LCOW/value"
#     ]

#     true_vals = [0.79001381, 0.79909567, 0.80820627]
#     d = data[()]
#     for i, tv in enumerate(true_vals):
#         assert d[i] == pytest.approx(tv, rel=1e-2)
#     data = h5file[
#         "ro_analysis/erd_type/pump_as_turbine/membrane_cost/outputs/fs.costing.LCOW/value"
#     ]

#     true_vals = [0.46957795, 0.50886109, 0.54814424]
#     d = data[()]
#     for i, tv in enumerate(true_vals):
#         assert d[i] == pytest.approx(tv, rel=1e-2)
#     h5file.close()
#     """test that backup works, will not run actual simulation ,create a back up file, and
#     load data from it into sim file. The lp.back_file_name should not be None
#     """


# def test_sweep_backup(loop_sweep_setup):
#     """test that backup works, will not run actual simulation ,create a back up file, and
#     load data from it into sim file. The lp.back_file_name should not be None
#     """
#     lp, test_file = loop_sweep_setup
#     lp.build_run_dict()
#     lp.run_simulations()
#     assert lp.h5_backup_location != None

#     h5file = h5py.File(
#         lp.h5_file_location_default + "_analysisType_ro_analysis.h5", "r"
#     )
#     data = h5file[
#         "ro_analysis/erd_type/pressure_exchanger/membrane_cost/outputs/fs.costing.LCOW/value"
#     ]

#     true_vals = [0.79001381, 0.79909567, 0.80820627]
#     d = data[()]
#     for i, tv in enumerate(true_vals):
#         assert d[i] == pytest.approx(tv, rel=1e-2)
#     data = h5file[
#         "ro_analysis/erd_type/pump_as_turbine/membrane_cost/outputs/fs.costing.LCOW/value"
#     ]
#     true_vals = [0.46957795, 0.50886109, 0.54814424]
#     d = data[()]
#     for i, tv in enumerate(true_vals):
#         assert d[i] == pytest.approx(tv, rel=1e-2)
#     h5file.close()
#     try:
#         os.remove(lp.h5_file_location_default + "_analysisType_ro_analysis.h5")
#     except OSError:
#         pass
#     try:
#         os.remove(lp.h5_backup_location)
#     except OSError:
#         pass


# Differentail parmater sweep tool does not work with current
# frame work, this will be incldued once it is.
def test_diff_setup(loop_diff_setup):
    lp, test_file = loop_diff_setup
    lp.build_run_dict()
    lp.run_simulations()
    h5file = h5py.File(
        lp.h5_file_location_default + "_analysisType_ro_analysis.h5", "r"
    )
    data = h5file[
        "ro_analysis/erd_type/pressure_exchanger/membrane_cost/outputs/fs.costing.LCOW/value"
    ]
    print(data[()])
