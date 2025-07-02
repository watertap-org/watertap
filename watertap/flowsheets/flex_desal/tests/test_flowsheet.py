#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
Tests for full Water Resource Recovery Facility with ASM2d and ADM1
(WRRF; a.k.a., wastewater treatment plant) flowsheet example with ASM1 and ADM1.
The flowsheet follows the same formulation as benchmark simulation model no.2 (BSM2)
but comprises different specifications for default values than BSM2.
"""

# Some more information about this module
__author__ = "Chenyu Wang"

import pandas as pd
import pyomo.environ as pyo
import pytest
import os

from watertap.flowsheets.flex_desal import flowsheet as fs
from watertap.flowsheets.flex_desal import utils
from watertap.flowsheets.flex_desal.params import FlexDesalParams
from watertap.flowsheets.flex_desal.price_taker_model import (
    PriceTakerModel,
)  # Use this until IDAES is updated
from watertap.core.solvers import get_solver

solver = get_solver()


@pytest.mark.requires_idaes_solver
class TestPriceTakerWorkflow:
    @pytest.fixture(scope="class")
    def system_frame(self):
        filepath = os.path.join("..", "sbce_pricesignal.csv")
        price_data = pd.read_csv(filepath)
        price_data["Energy Rate"] = (
            price_data["electric_energy_0_2022-07-05_2022-07-14_0"]
            + price_data["electric_energy_1_2022-07-05_2022-07-14_0"]
            + price_data["electric_energy_2_2022-07-05_2022-07-14_0"]
            + price_data["electric_energy_3_2022-07-05_2022-07-14_0"]
        )
        price_data["Fixed Demand Rate"] = price_data[
            "electric_demand_maximum_2022-07-05_2022-07-14_0"
        ]
        price_data["Var Demand Rate"] = price_data[
            "electric_demand_peak-summer_2022-07-05_2022-07-14_0"
        ]
        price_data["Emissions Intensity"] = 0
        price_data["Customer Cost"] = price_data[
            "electric_customer_0_2022-07-05_2022-07-14_0"
        ]

        m = PriceTakerModel()

        m.params = FlexDesalParams(
            start_date="2022-07-05 00:00:00",
            end_date="2022-07-15 00:00:00",
            annual_production_AF=3125,  # acrft/yr
        )
        m.params.intake.nominal_flowrate = 1063.5  # m3/hr
        m.params.ro.update(
            {
                "startup_delay": 8,  # hours
                "minimum_downtime": 4,  # hours
                "nominal_flowrate": 337.670,  # m3/hr
                "surrogate_type": "quadratic_surrogate",
                "surrogate_a": 11.509,
                "surrogate_b": -10.269,
                "surrogate_c": 5.627,
                "surrogate_d": 0,
                "minimum_recovery": 0.4,
                "nominal_recovery": 0.465,
                "maximum_recovery": 0.52,
                "allow_variable_recovery": True,
            }
        )

        # Append LMP data to the model
        m.append_lmp_data(lmp_data=price_data["Energy Rate"])

        m.build_multiperiod_model(
            flowsheet_func=fs.build_desal_flowsheet,
            flowsheet_options={"params": m.params},
        )

        # Update the time-varying parameters other than the LMP, such as
        # demand costs and emissions intensity. LMP value is updated by default
        m.update_operation_params(
            {
                "fixed_demand_rate": price_data["Fixed Demand Rate"],
                "variable_demand_rate": price_data["Var Demand Rate"],
                "emissions_intensity": price_data["Emissions Intensity"],
                "customer_cost": price_data["Customer Cost"],
            }
        )

        return m, price_data

    @pytest.mark.unit
    def test_price_data_structure(self, system_frame):
        m, price_data = system_frame

        assert "Energy Rate" in price_data.columns
        assert "Fixed Demand Rate" in price_data.columns
        assert "Var Demand Rate" in price_data.columns
        assert "Emissions Intensity" in price_data.columns
        assert "Customer Cost" in price_data.columns
        assert "electric_energy_0_2022-07-05_2022-07-14_0" in price_data.columns
        assert "electric_energy_1_2022-07-05_2022-07-14_0" in price_data.columns
        assert "electric_energy_2_2022-07-05_2022-07-14_0" in price_data.columns
        assert "electric_energy_3_2022-07-05_2022-07-14_0" in price_data.columns
        assert "electric_demand_maximum_2022-07-05_2022-07-14_0" in price_data.columns
        assert (
            "electric_demand_peak-summer_2022-07-05_2022-07-14_0" in price_data.columns
        )
        assert "electric_customer_0_2022-07-05_2022-07-14_0" in price_data.columns

        assert hasattr(m.params.intake, "nominal_flowrate")
        assert hasattr(m.params.ro, "startup_delay")
        assert hasattr(m.params.ro, "minimum_downtime")
        assert hasattr(m.params.ro, "nominal_flowrate")
        assert hasattr(m.params.ro, "surrogate_type")
        assert hasattr(m.params.ro, "surrogate_a")
        assert hasattr(m.params.ro, "surrogate_b")
        assert hasattr(m.params.ro, "surrogate_c")
        assert hasattr(m.params.ro, "surrogate_d")
        assert hasattr(m.params.ro, "minimum_recovery")
        assert hasattr(m.params.ro, "nominal_recovery")
        assert hasattr(m.params.ro, "maximum_recovery")
        assert hasattr(m.params.ro, "allow_variable_recovery")

    @pytest.mark.unit
    def test_add_constraints(self, system_frame):
        m, price_data = system_frame

        # Add demand cost and fixed cost calculation constraints
        fs.add_demand_and_fixed_costs(m)

        assert isinstance(m.fixed_demand_cost, pyo.Var)
        assert isinstance(m.variable_demand_cost, pyo.Var)
        assert isinstance(m.fixed_monthly_cost, pyo.Var)
        assert isinstance(m.calculate_fixed_demand_cost, pyo.Constraint)
        assert isinstance(m.calculate_variable_demand_cost, pyo.Constraint)
        assert isinstance(m.calculate_fixed_monthly_cost, pyo.Constraint)

        # Add the startup delay constraints
        fs.add_delayed_startup_constraints(m)
        assert isinstance(m.posttreatment_unit_commitment, pyo.Constraint)
        assert isinstance(m.brine_pump_unit_commitment, pyo.Constraint)

    @pytest.mark.unit
    def test_add_capacity_limits(self, system_frame):
        m, price_data = system_frame

        for skid in range(1, m.params.ro.num_ro_skids + 1):
            m.add_capacity_limits(
                op_block_name=f"reverse_osmosis.ro_skid[{skid}]",  # Name of the operation model block
                commodity="feed_flowrate",
                # Name of the commodity on the operation model that capacity constraints will be applied to (flow rate, power, etc.)
                capacity=m.params.ro.nominal_flowrate,  # Maximum capacity on the commodity
                op_range_lb=1,
                # Ratio of the capacity at minimum stable operation to the maximum capacity. Must be between [0, 1]
            )

    @pytest.mark.unit
    def test_constrain_water_production(self, system_frame):
        m, price_data = system_frame

        m.total_water_production = pyo.Expression(
            expr=m.params.timestep_hours
            * sum(m.period[:, :].posttreatment.product_flowrate)
        )
        m.total_energy_cost = pyo.Expression(expr=sum(m.period[:, :].energy_cost))
        m.total_demand_cost = pyo.Expression(
            expr=m.fixed_demand_cost + m.variable_demand_cost
        )
        m.total_customer_cost = pyo.Expression(
            expr=sum(m.period[:, :].customer_cost) * m.params.num_months
        )
        m.total_electricity_cost = pyo.Expression(
            expr=m.total_energy_cost + m.total_demand_cost + m.total_customer_cost
        )

        # Feed flow to the intake does not vary with time
        m.fix_operation_var("intake.feed_flowrate", m.params.intake.nominal_flowrate)
        # Pretreatment is either active (1) or inactive (0) for the entire run
        m.fix_operation_var("pretreatment.op_mode", 1)

        fs.constrain_water_production(m)

        # If water recovery is static, it must be fixed
        if not m.params.ro.allow_variable_recovery:
            utils.fix_recovery(m, recovery=m.params.ro.nominal_recovery)

    @pytest.mark.unit
    def test_update_recovery_bounds(self, system_frame):
        m, price_data = system_frame
        utils.update_recovery_bounds(m, lb=0.4, ub=0.5)

    @pytest.mark.unit
    def test_get_baseline_model(self, system_frame):
        m, price_data = system_frame
        utils.get_baseline_model(m)

    @pytest.mark.component
    @pytest.mark.xfail
    # This test will fail if the user does not have a Gurobi license
    def test_gurobi_solve(self, system_frame):
        m, price_data = system_frame

        solver = pyo.SolverFactory("gurobi")
        solver.options["MIPGap"] = 0.03
        solver.solve(m)

    @pytest.mark.component
    @pytest.mark.xfail
    # This test will fail if the user does not have a Gurobi license
    def test_gurobi_util_solve(self, system_frame):
        m, price_data = system_frame

        solver = utils.get_gurobi_solver_model(m)
        solver.solve(m)

    # @pytest.mark.component
    # def test_scip_solve(self, system_frame):
    #     # Timed out after 2s
    #     m, price_data = system_frame
    #
    #     solver = pyo.SolverFactory("scip")
    #     solver.solve(m)

    # @pytest.mark.component
    # # Took 6 hours to solve locally
    # def test_baron_solve(self, system_frame):
    #     m, price_data = system_frame
    #
    #     solver = pyo.SolverFactory("gams")
    #     solver.solve(
    #         m, solver="baron", add_options=[f"options optcr={0.03};"]
    #     )

    # @pytest.mark.component
    # # cplex unsuitable for model type (minlp)
    # def test_cplex_solve(self, system_frame):
    #     m, price_data = system_frame
    #
    #     solver = pyo.SolverFactory("gams")
    #     solver.solve(
    #         m, solver="cplex", add_options=[f"options optcr={0.03};"]
    #     )

    # @pytest.mark.component
    # def test_results(self, system_frame):
    #     m = system_frame
    #
    #     assert value(m.fs.Treated.properties[0].flow_vol) == pytest.approx(
    #         0.24219, rel=1e-3
    #     )
    #     assert value(m.fs.Treated.properties[0].conc_mass_comp["S_A"]) == pytest.approx(
    #         2.7392e-06, abs=1e-6
    #     )
    #     assert value(m.fs.Treated.properties[0].conc_mass_comp["S_F"]) == pytest.approx(
    #         0.00026922, rel=1e-3
    #     )
    #     assert value(m.fs.Treated.properties[0].conc_mass_comp["S_I"]) == pytest.approx(
    #         0.057450006, rel=1e-3
    #     )
    #     assert value(
    #         m.fs.Treated.properties[0].conc_mass_comp["S_N2"]
    #     ) == pytest.approx(0.0516457, rel=1e-3)
    #     assert value(
    #         m.fs.Treated.properties[0].conc_mass_comp["S_NH4"]
    #     ) == pytest.approx(0.000209, rel=1e-3)
    #     assert value(
    #         m.fs.Treated.properties[0].conc_mass_comp["S_NO3"]
    #     ) == pytest.approx(0.00452157, rel=1e-3)
    #     assert value(
    #         m.fs.Treated.properties[0].conc_mass_comp["S_O2"]
    #     ) == pytest.approx(0.0014803, rel=1e-3)
    #     assert value(
    #         m.fs.Treated.properties[0].conc_mass_comp["S_PO4"]
    #     ) == pytest.approx(0.755433, rel=1e-3)
    #     assert value(m.fs.Treated.properties[0].conc_mass_comp["S_K"]) == pytest.approx(
    #         0.3667916, rel=1e-3
    #     )
    #     assert value(
    #         m.fs.Treated.properties[0].conc_mass_comp["S_Mg"]
    #     ) == pytest.approx(0.0182828, rel=1e-3)
    #     assert value(
    #         m.fs.Treated.properties[0].conc_mass_comp["S_IC"]
    #     ) == pytest.approx(0.1497356, rel=1e-3)
    #     assert value(
    #         m.fs.Treated.properties[0].conc_mass_comp["X_AUT"]
    #     ) == pytest.approx(0.0004246397, rel=1e-3)
    #     assert value(m.fs.Treated.properties[0].conc_mass_comp["X_H"]) == pytest.approx(
    #         0.01328946, rel=1e-3
    #     )
    #     assert value(m.fs.Treated.properties[0].conc_mass_comp["X_I"]) == pytest.approx(
    #         0.0120139, rel=1e-3
    #     )
    #     assert value(
    #         m.fs.Treated.properties[0].conc_mass_comp["X_PAO"]
    #     ) == pytest.approx(0.01305288, rel=1e-3)
    #     assert value(
    #         m.fs.Treated.properties[0].conc_mass_comp["X_PHA"]
    #     ) == pytest.approx(7.7306e-06, abs=1e-6)
    #     assert value(
    #         m.fs.Treated.properties[0].conc_mass_comp["X_PP"]
    #     ) == pytest.approx(0.0043593, rel=1e-3)
    #     assert value(m.fs.Treated.properties[0].conc_mass_comp["X_S"]) == pytest.approx(
    #         0.00021958, rel=1e-3
    #     )
    #
    #     # Check electricity consumption for each aerobic reactor
    #     assert value(m.fs.R5.electricity_consumption[0]) == pytest.approx(
    #         48.0849, rel=1e-3
    #     )
    #     assert value(m.fs.R6.electricity_consumption[0]) == pytest.approx(
    #         48.0849, rel=1e-3
    #     )
    #     assert value(m.fs.R7.electricity_consumption[0]) == pytest.approx(
    #         48.0849, rel=1e-3
    #     )
    #
