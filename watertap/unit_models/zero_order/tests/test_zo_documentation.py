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
"""
Tests that documentation exists for all zero-order unit models
"""
import pytest
import os

rst_path = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..",
    "..",
    "..",
    "..",
    "docs",
    "technical_reference",
    "unit_models",
    "zero_order_unit_models",
)
rst_list = []
exclude_rst_files = ["index.rst"]
for file in os.listdir(rst_path):
    filename = os.fsdecode(file)
    if filename.endswith(".rst") and filename not in exclude_rst_files:
        rst_list.append(filename[:-4])

py_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
py_list = []
exclude_py_files = ["__init__.py"]
for file in os.listdir(py_path):
    filename = os.fsdecode(file)
    if filename.endswith(".py") and filename not in exclude_py_files:
        py_list.append(filename[:-3])


@pytest.mark.unit
def test_rst_path():
    assert os.path.exists(rst_path) == True

    assert rst_path == os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..",
        "..",
        "..",
        "..",
        "docs",
        "technical_reference",
        "unit_models",
        "zero_order_unit_models",
    )

    for file in os.listdir(rst_path):
        filename = os.fsdecode(file)
        if filename.endswith(".rst") and filename not in exclude_rst_files:
            file_path = os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                "..",
                "..",
                "..",
                "..",
                "docs",
                "technical_reference",
                "unit_models",
                "zero_order_unit_models",
                f"{file}",
            )
            assert os.path.isfile(file_path) == True


@pytest.mark.unit
def test_py_path():
    assert os.path.exists(py_path) == True

    assert py_path == os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

    for file in os.listdir(py_path):
        filename = os.fsdecode(file)
        if filename.endswith(".py") and filename not in exclude_py_files:
            file_path = os.path.join(
                os.path.dirname(os.path.abspath(__file__)), "..", f"{file}"
            )
            assert os.path.isfile(file_path) == True


@pytest.mark.integration
def test_lists_match():
    # Run pytest -vv to see how the lists are different
    assert sorted(py_list) == sorted(rst_list)

    assert len(py_list) == len(rst_list)
