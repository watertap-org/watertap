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

doc_path = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..",
    "..",
    "..",
    "..",
    "docs",
)

rst_path = os.path.join(
    doc_path,
    "technical_reference",
    "unit_models",
    "zero_order_unit_models",
)
py_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def generate_model_list(path, extension, exclude):
    file_list = []
    for file in os.listdir(path):
        filename = os.fsdecode(file)
        if filename.endswith(extension) and filename not in exclude:
            file_path = os.path.join(
                path,
                f"{file}",
            )
            assert os.path.isfile(file_path)
            file_list.append(filename[: -len(extension)])
    return file_list


@pytest.mark.integration
@pytest.mark.skipif(not os.path.exists(doc_path), reason="No docs directory")
def test_lists_match():

    rst_list = generate_model_list(rst_path, ".rst", ("index.rst",))
    py_list = generate_model_list(py_path, ".py", ("__init__.py",))

    # Run pytest -vv to see how the lists are different
    assert sorted(py_list) == sorted(rst_list)
