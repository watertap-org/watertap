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
Tests for meta_api module
"""
import pytest
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.metab import (
    metab_api,
)

WA = metab_api.WorkflowActions


@pytest.mark.unit
def test_flowsheet_interface():
    r = metab_api.flowsheet_interface()
    assert r is not None


@pytest.mark.unit
def test_build():
    fsi = metab_api.flowsheet_interface()
    fsi.run_action(WA.build)


@pytest.mark.component
def test_solve():
    fsi = metab_api.flowsheet_interface()
    fsi.run_action(WA.solve)


@pytest.mark.unit
def test_io(tmp_path):
    fsi = metab_api.flowsheet_interface()
    fsi.run_action(WA.build)
    filename = "test_io.json"
    with (tmp_path / filename).open("w") as f:
        fsi.save(f)
    with (tmp_path / filename).open("r") as f:
        fsi.load(f)


@pytest.mark.unit
def test_action_set_was_run():
    WA = metab_api.WorkflowActions
    interface = metab_api.FlowsheetInterface(None)
    # pretend we ran all the steps
    action_list = (WA.build, WA.solve, WA.results)
    for action in action_list:
        interface._action_set_was_run(action)
    # They should all be marked as run
    for action in action_list:
        assert interface._action_was_run(action)
    # Now saying we re-ran the build step should reset run status on solve and results
    interface._action_set_was_run(WA.build)
    for action in action_list:
        if action == WA.build:
            assert interface._action_was_run(action)
        else:
            assert not interface._action_was_run(action)
