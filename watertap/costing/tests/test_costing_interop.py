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

import pyomo.environ as pyo
from idaes.core import UnitModelCostingBlock

from watertap.costing import WaterTAPCosting, ZeroOrderCosting
import watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.amo_1575_hrcs.hrcs as hrcs

hrcs_original_costing = hrcs.add_costing


def simple_main():
    m = hrcs.build()

    hrcs.set_operating_conditions(m)
    hrcs.assert_degrees_of_freedom(m, 0)
    hrcs.assert_units_consistent(m)

    hrcs.initialize_system(m)

    results = hrcs.solve(m, checkpoint="solve flowsheet after initializing system")

    hrcs.add_costing(m)
    m.fs.costing.initialize()

    hrcs.assert_degrees_of_freedom(m, 0)
    hrcs.assert_units_consistent(m)
    results = hrcs.solve(m, checkpoint="solve flowsheet with costing")

    return m, results


def add_wt_costing(m):
    m.fs.costing = WaterTAPCosting()

    m.fs.costing.register_flow_type(
        "ferric_chloride", 0.6 * pyo.units.USD_2020 / pyo.units.kg
    )

    m.fs.mixer.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    m.fs.HRCS.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.clarifier.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.sep.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    m.fs.costing.cost_process()

    m.fs.costing.add_electricity_intensity(m.fs.feed.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.feed.properties[0].flow_vol)

    m.fs.costing.add_annual_water_production(m.fs.product.properties[0].flow_vol)

    m.fs.costing.add_annual_water_production(
        m.fs.feed.properties[0].flow_vol, "annual_water_inlet"
    )

    m.fs.costing.LCOW_with_revenue = pyo.Expression(
        expr=(
            (
                m.fs.costing.total_annualized_cost
                + m.fs.costing.aggregate_flow_costs["ferric_chloride"]
            )
            / m.fs.costing.annual_water_production
        ),
        doc="Levelized Cost of Water With Revenue",
    )

    m.fs.costing.add_LCOW(m.fs.feed.properties[0].flow_vol, "LCOT")

    m.fs.costing.LCOT_with_revenue = pyo.Expression(
        expr=(
            (
                m.fs.costing.total_annualized_cost
                + m.fs.costing.aggregate_flow_costs["ferric_chloride"]
            )
            / m.fs.costing.annual_water_inlet
        ),
        doc="Levelized Cost of Treatment With Revenue",
    )


def test_hrcs_case_1575_wtcosting():
    hrcs.add_costing = add_wt_costing

    try:
        m, results = simple_main()
        pyo.assert_optimal_termination(results)

        # check costing -- baseline is 0.02003276
        assert pyo.value(m.fs.costing.LCOW) == pytest.approx(
            0.02419106, rel=1e-3
        )  # in $/m**3
    except:
        raise
    finally:
        hrcs.add_costing = hrcs_original_costing


def add_zo_costing(m):
    source_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..",
        "..",
        "examples",
        "flowsheets",
        "case_studies",
        "wastewater_resource_recovery",
        "amo_1575_hrcs",
        "hrcs_case_1575.yaml",
    )

    m.fs.costing = ZeroOrderCosting(case_study_definition=source_file)

    m.fs.mixer.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    m.fs.HRCS.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.clarifier.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.sep.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    m.fs.costing.cost_process()

    m.fs.costing.add_electricity_intensity(m.fs.feed.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.feed.properties[0].flow_vol)

    m.fs.costing.add_annual_water_production(m.fs.product.properties[0].flow_vol)

    m.fs.costing.add_annual_water_production(
        m.fs.feed.properties[0].flow_vol, "annual_water_inlet"
    )

    m.fs.costing.LCOW_with_revenue = pyo.Expression(
        expr=(
            (
                m.fs.costing.total_annualized_cost
                + m.fs.costing.aggregate_flow_costs["ferric_chloride"]
            )
            / m.fs.costing.annual_water_production
        ),
        doc="Levelized Cost of Water With Revenue",
    )

    m.fs.costing.add_LCOW(m.fs.feed.properties[0].flow_vol, "LCOT")

    m.fs.costing.LCOT_with_revenue = pyo.Expression(
        expr=(
            (
                m.fs.costing.total_annualized_cost
                + m.fs.costing.aggregate_flow_costs["ferric_chloride"]
            )
            / m.fs.costing.annual_water_inlet
        ),
        doc="Levelized Cost of Treatment With Revenue",
    )


def test_hrcs_case_1575_zocosting():
    hrcs.add_costing = add_zo_costing

    try:
        m, results = simple_main()
        pyo.assert_optimal_termination(results)

        # check costing -- baseline is 0.02003276
        assert pyo.value(m.fs.costing.LCOW) == pytest.approx(
            0.02172723, rel=1e-3
        )  # in $/m**3
    except:
        raise
    finally:
        hrcs.add_costing = hrcs_original_costing
