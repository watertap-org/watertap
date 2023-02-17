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
Tests that documentation exists for all zero-order unit models.
To create new documentation, update the unit classification excel sheet and run the automation script.
"""
import pytest
import os

from pyomo.common.fileutils import this_file_dir
from pathlib import Path

this_file_dir = this_file_dir()

doc_path = os.path.join(
    this_file_dir,
    "..",
    "..",
    "..",
    "..",
    "docs",
)

doc_path = Path(doc_path).resolve()

rst_path = os.path.join(
    doc_path,
    "technical_reference",
    "unit_models",
    "zero_order_unit_models",
)

py_path = os.path.join(this_file_dir, "..")


def generate_model_list(path, extension, exclude):
    file_list = []
    for file in os.listdir(path):
        filename = os.fsdecode(file)
        if filename.endswith(extension) and filename not in exclude:
            file_path = os.path.join(
                path,
                file,
            )
            assert os.path.isfile(file_path)
            file_list.append(filename[: -len(extension)])
    return file_list


@pytest.mark.integration
@pytest.mark.skipif(
    not os.path.exists(doc_path / "index.rst"), reason="No RST files found"
)
def test_lists_match():

    rst_list = generate_model_list(rst_path, ".rst", ("index.rst",))
    py_list = generate_model_list(py_path, ".py", ("__init__.py",))

    # remove the AOPMixin
    py_list.remove("aop_addition_zo")

    # Run pytest -vv to see how the lists are different
    assert sorted(py_list) == sorted(rst_list)
