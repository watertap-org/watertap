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

from watertap.tools.analysis_tools.loop_tool import (
    data_merging_tool,
)
import h5py

__author__ = "Alexander V. Dudchenko (SLAC)"

_this_file_path = os.path.dirname(os.path.abspath(__file__))


def create_test_h5_file():
    f_name = _this_file_path + "/test_h5_file.h5"
    test_group = "test_vals"
    h5file = h5py.File(f_name, "w")
    h5file[test_group + "/solve_successful/solve_successful"] = np.ones(10)
    h5file.close()
    return f_name, test_group
    # test creating new file, no back exists, so we should run sweep


@pytest.mark.component
def test_new_file_run():
    """test create new file, and ensure we run a sweep"""
    h5_file, test_group = create_test_h5_file()
    run_sweep = data_merging_tool.merge_data_into_file(
        h5_file,
        None,
        test_group,
        expected_solved_values=10,
        min_solve_values=10,
        force_rerun=None,
    )
    try:
        os.remove(h5_file)
    except OSError:
        pass
    assert run_sweep == True


@pytest.mark.component
def test_runs_from_backup():
    h5_file, test_group = create_test_h5_file()
    backup_name = data_merging_tool.create_backup_file(h5_file, None, "test_vals")
    run_sweep = data_merging_tool.merge_data_into_file(
        h5_file,
        backup_name,
        test_group,
        expected_solved_values=10,
        min_solve_values=None,
        force_rerun=None,
    )
    # solutions shouild already be present in solved file, so we re not
    # reruning
    assert run_sweep == False

    # we expect more solutions then present in file, so
    # run_sweep shold be true
    run_sweep = data_merging_tool.merge_data_into_file(
        h5_file,
        backup_name,
        test_group,
        expected_solved_values=20,
        min_solve_values=None,
        force_rerun=None,
    )
    assert run_sweep == True

    # we expect min 9 solutions, so should not re-run as we have 10
    run_sweep = data_merging_tool.merge_data_into_file(
        h5_file,
        backup_name,
        test_group,
        expected_solved_values=20,
        min_solve_values=9,
        force_rerun=None,
    )
    assert run_sweep == False

    # We do not want to re-run, so should not run_sweep
    run_sweep = data_merging_tool.merge_data_into_file(
        h5_file,
        backup_name,
        test_group,
        expected_solved_values=20,
        min_solve_values=10,
        force_rerun=False,
    )
    assert run_sweep == False

    # We are forcing a re-run even thouhg there as many values as expected
    run_sweep = data_merging_tool.merge_data_into_file(
        h5_file,
        backup_name,
        test_group,
        expected_solved_values=10,
        min_solve_values=10,
        force_rerun=True,
    )
    assert run_sweep == False
    try:
        os.remove(h5_file)
        os.remove(backup_name)
    except OSError:
        pass
