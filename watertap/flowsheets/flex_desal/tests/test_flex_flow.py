#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
import os
import importlib.metadata

import pyomo.environ as pyo
import pytest
import pandas as pd

from idaes.apps.grid_integration import PriceTakerModel

from watertap.flowsheets.flex_desal import wrd_ro_flowsheet as fs
from watertap.flowsheets.flex_desal import utils
from watertap.flowsheets.flex_desal.params import FlexDesalParams
from idaes.core.solvers import get_solver


@pytest.mark.requires_idaes_solver
class TestPriceTakerWorkflow:
    @pytest.fixture(scope="class")
    @classmethod
    def system_frame(cls):
        price_data_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..",
            "wrd_pricesignal_summer_week_DR.csv",
        )
        price_data = pd.read_csv(price_data_path)
        price_data["Energy Rate"] = (
            price_data["electric_energy_on_peak"]
            + price_data["electric_energy_mid_peak"]
            + price_data["electric_energy_off_peak"]
            + price_data["electric_energy_super_off_peak"]
        )
        price_data["Fixed Demand Rate"] = price_data["electric_demand_fixed"]
        price_data["Var Demand Rate"] = price_data["electric_demand_peak"]
        price_data["Customer Cost"] = price_data["electric_customer_fixed_charge"]
        price_data["Demand_Response_Price"] = price_data[
            "electric_demand_response_price"
        ]

        price_data["Emissions Intensity"] = 0
        peak_hours = price_data["Var Demand Rate"].to_numpy() != 0  # For plotting

        # Load PV data
        pv_kW = price_data["solar_output_kW"]
        pv_capacity = max(pv_kW)
        pv_capacity_factors = pv_kW / pv_capacity

        m = PriceTakerModel()
        m.params = FlexDesalParams(
            start_date="2022-07-05 00:00:00",
            end_date="2022-07-15 00:00:00",
            annual_production_AF=4000,  # Very low water production to ensure feasiblity. The additional constraints force longer periods of downtime
            timestep_hours=1,
            include_onsite_solar=True,
            onsite_capacity=pv_capacity,
            nonworking_hours=list(range(0, 8))
            + list(
                range(18, 24)
            ),  # 6pm-8am are nonworking hours (assuming time index starts at 0 for 12am-1am)
            CAPEX_yr=6498300,  # For WRD, this assumes a 30 yr lifetime
            include_demand_response=True,
            max_daily_shutdowns=1,
        )

        m.params.intake.update(
            {
                "energy_intensity": 0,
                "nominal_flowrate": 2500,  # m3/hr
                "feed_cost": 0.16,
                "chemical_cost": 0.0332,
            }
        )

        m.params.wrd_uf.update(
            {
                "minimum_downtime": 2,
                "startup_delay": 2,
                "minimum_flowrate": 344,  # m3/hr
                "nominal_flowrate": 900,
                "maximum_flowrate": 989,
                "surrogate_type": "quadratic_energy_intensity",
                "surrogate_a": 2.71e-1,
                "surrogate_b": -3.32e-4,
                "surrogate_c": 2.39e-7,
                "nominal_recovery": 0.96,
                "num_uf_pumps": 3,
            }
        )

        m.params.wrd_ro.update(
            {
                "startup_delay": 2,  # hours
                "minimum_downtime": 2,
                "minimum_flowrate": 520,  # m3/hr
                "nominal_flowrate": 602,
                "maximum_flowrate": 635,
                "allow_variable_recovery": True,
                "surrogate_type": "PySMO_polyfit",
                "surrogate_file": os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
                    "ro_SEC_poly_fit_order_2.json",
                ),
                "minimum_recovery": 0.88,
                "nominal_recovery": 0.925,
                "maximum_recovery": 0.925,
                "num_ro_skids": 4,
                "replacement_types": ["membranes", "motors"],
                "replacement_costs": [
                    500 * 4 * (72 + 30 + 15),
                    125000,
                ],  # $ per replacement
                "replacement_lifetimes": [5, 20],  # years
                "replacement_max_flex_penalty": [
                    0.1,
                    0.1,
                ],  # Reduction in lifetime if shutdowns occur every day
            }
        )

        m.params.posttreatment.update(
            {
                "energy_intensity": 0.101,  # kWh/m3
                "chemical_cost": 0.0310,  # $/m3
            }
        )

        m.params.brinedischarge.update({"brine_cost": 0.43, "energy_intensity": 0})

        m.append_lmp_data(lmp_data=price_data["Energy Rate"])

        m.build_multiperiod_model(
            flowsheet_func=fs.build_wrd_desal_flowsheet,
            flowsheet_options={"params": m.params},
        )

        m.update_operation_params(
            {
                "fixed_demand_rate": price_data["Fixed Demand Rate"],
                "variable_demand_rate": price_data["Var Demand Rate"],
                "emissions_intensity": price_data["Emissions Intensity"],
                "customer_cost": price_data["Customer Cost"],
                "demand_response_price": price_data["Demand_Response_Price"],
            }
        )

        m.update_operation_params(
            {"power_generation.capacity_factor": pv_capacity_factors}
        )

        return m, price_data, peak_hours

    @pytest.mark.unit
    def test_build(self, system_frame):
        m, price_data, peak_hours = system_frame

        # Check Demand Response
        assert "Demand_Response_Price" in price_data.columns

        # Check params added
        assert hasattr(m.params.wrd_ro, "surrogate_file")
        assert hasattr(m.params.wrd_ro, "replacement_types")
        assert hasattr(m.params.wrd_ro, "replacement_costs")
        assert hasattr(m.params.wrd_ro, "replacement_lifetimes")
        assert hasattr(m.params.wrd_ro, "replacement_max_flex_penalty")
        assert hasattr(m.params.wrd_uf, "surrogate_a")
        assert hasattr(m.params.wrd_uf, "surrogate_b")
        assert hasattr(m.params.wrd_uf, "surrogate_c")
        assert hasattr(m.params.wrd_uf, "num_uf_pumps")
        assert hasattr(m.params.posttreatment, "chemical_cost")
        assert hasattr(m.params.brinedischarge, "brine_cost")
        assert hasattr(m.params.intake, "chemical_cost")
        assert hasattr(m.params.intake, "feed_cost")

        for blk in m.period.values():
            # Check PV is added
            assert hasattr(blk, "power_generation")
            # Check Demand Response price added
            assert hasattr(blk, "demand_response_price")

    @pytest.mark.unit
    def test_add_constraints(self, system_frame):
        m, price_data, peak_hours = system_frame

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

        # Ensure consistent ending and starting states of the plant
        fs.begin_and_end_constraint(m)
        assert isinstance(m.match_train_1_at_start_and_end, pyo.Constraint)

        # Limit the number of shutdowns per day
        fs.add_maximum_shutdowns(m)
        assert hasattr(m, "max_shutdowns_per_24h_window")

        # Add the slow shutdown constraint
        fs.add_delayed_shutdown_constraints(m)
        assert hasattr(m, "posttreatment_unit_commitment_shutdown")

        # Limit hours when start-ups and shutdowns can occur to reflect when operators are on-site
        fs.add_working_hours_constraint(m)

        # Limit number of trains that can turn on and off
        fs.restrict_flexible_trains(m, num_flexible_trains=2)
        for p in m.period:
            for skid in [1, 2]:
                ro_skid = m.period[p].reverse_osmosis.ro_skid[skid]
                assert ro_skid.startup.fixed
                assert ro_skid.startup() == 0
                assert ro_skid.shutdown.fixed
                assert ro_skid.shutdown() == 0

    @pytest.mark.unit
    def test_add_expressions(self, system_frame):
        m, price_data, peak_hours = system_frame

        fs.add_useful_expressions(m)
        fs.add_flow_costs(m)  #  Flow costs = Feed, Brine, and Chemicals

        m.total_water_production = pyo.Expression(
            expr=m.params.timestep_hours
            * sum(m.period[:, :].posttreatment.product_flowrate)
        )
        m.total_energy_cost = pyo.Expression(expr=sum(m.period[:, :].energy_cost))

        # Demand costs are automatically normalized by number of months.
        # If the time horizon is a week, the cost is multiplied by 7/31
        m.total_demand_cost = pyo.Expression(
            expr=m.fixed_demand_cost + m.variable_demand_cost
        )
        m.total_customer_cost = pyo.Expression(
            expr=sum(m.period[:, :].customer_cost) * m.params.num_months
        )

        m.total_op_cost = pyo.Expression(
            expr=m.total_energy_cost
            + m.total_demand_cost
            + m.total_customer_cost
            - m.total_demand_response_revenue
            + m.total_feed_cost
            + m.total_brine_cost
            + m.total_chemical_cost
        )

        # add CAPEX as a fixed cost to calculate LCOW
        m.fixed_cost = pyo.Expression(expr=m.params.CAPEX_yr * m.params.num_months / 12)
        m.total_cost = pyo.Expression(expr=m.total_op_cost + m.fixed_cost)

        m.LCOW = pyo.Expression(expr=m.total_cost / m.total_water_production)  # $/m3

        assert isinstance(m.total_water_production, pyo.Expression)
        assert isinstance(m.total_energy_cost, pyo.Expression)
        assert isinstance(m.total_demand_cost, pyo.Expression)
        assert isinstance(m.total_customer_cost, pyo.Expression)
        assert isinstance(m.total_op_cost, pyo.Expression)
        assert isinstance(m.fixed_cost, pyo.Expression)
        assert isinstance(m.total_cost, pyo.Expression)
        assert isinstance(m.LCOW, pyo.Expression)

    @pytest.mark.unit
    def test_fixing_operations(self, system_frame):
        m, price_data, peak_hours = system_frame
        fs.fix_operations_for_first_four_days(m, peak_hours=peak_hours)

        utils.wrd_fix_ro_recovery(
            m,
            ro_recovery=m.params.wrd_ro.nominal_recovery,
        )

        utils.wrd_fix_uf_recovery(
            m,
            uf_recovery=m.params.wrd_uf.nominal_recovery,
        )

        # Set the water production target over the entire time period
        fs.constrain_water_production(m)

    @pytest.mark.unit
    def test_flow_changes_penalty(self, system_frame):
        m, price_data, peak_hours = system_frame

        fs.add_flow_changes_penalty(m)
        assert isinstance(m.flow_changes_penalty, pyo.Expression)

        m.obj = pyo.Objective(
            expr=1e-4
            * (
                m.total_energy_cost
                + m.total_demand_cost
                + m.total_customer_cost
                - m.total_demand_response_revenue
                + m.total_feed_cost
                + m.total_brine_cost
                + m.total_chemical_cost
                + m.flow_changes_penalty
            ),
            sense=pyo.minimize,
        )

    @pytest.mark.component
    # This test uses ipopt to solve the model, but it will give a nonsense solution
    # with binary variable between 0 and 1. However, the solution allows the post_solve test
    # to be performed
    def test_ipopt_solve(self, system_frame):
        m, price_data, peak_hours = system_frame

        solver = get_solver(
            options={
                "bound_push": 1e-8,
                "bound_frac": 1e-8,
                "tol": 1e-6,
                "max_iter": 5000,
            }
        )
        results = solver.solve(m, tee=True)

        # This is line that is causing the failure on Linux. But I don't know why.
        # results = solver.solve(m)

        pyo.assert_optimal_termination(results)

    @pytest.mark.component
    @pytest.mark.xfail
    # This test will fail if the user does not have a Gurobi license
    def test_gurobi_solve(self, system_frame):
        m, price_data, peak_hours = system_frame

        solver = pyo.SolverFactory("gurobi_direct_minlp")
        solver.options["MIPGap"] = 0.03
        results = solver.solve(m)

        pyo.assert_optimal_termination(results)

    @pytest.mark.unit
    @pytest.mark.xfail
    # This test will fail if the model is not solved aready
    def test_post_solve_calculations(self, system_frame):
        m, price_data, peak_hours = system_frame

        fs.calculate_replacement_costs(m)
        assert hasattr(m, "degree_of_flex")
        assert isinstance(m.total_replacement_cost, pyo.Expression)

        # These baseline values are not realistic- just for testing
        fs.calculate_flexibility_metrics(
            m,
            baseline_power=600,
            baseline_electricity_cost=20843,  # $/kWh
            baseline_replacement_cost=992,
        )
        assert isinstance(m.maximum_power, pyo.Var)
        assert isinstance(m.energy_capacity, pyo.Var)
        assert isinstance(m.power_capacity, pyo.Var)
        assert isinstance(m.LVOF, pyo.Var)

        design_var_values = m.get_design_var_values()
        filtered_design_var_values = {
            k: v
            for k, v in design_var_values.items()
            if "flow_change" not in k
            and "flow_changed" not in k
            and "reduction" not in k
        }

        assert filtered_design_var_values
