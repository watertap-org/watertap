import os
import glob

import numpy as np
import pandas as pd

from matplotlib.pyplot import subplots

from pyomo.common.config import ConfigDict, ConfigValue, In, Path, NonNegativeFloat
from pyomo.environ import (
    ConcreteModel,
    value,
    TransformationFactory,
    units as pyunits,
)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent

from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.generic_models.unit_models import Product

from idaes.core import FlowsheetBlock
from idaes.generic_models.costing import UnitModelCostingBlock
import idaes.logger as idaeslogger

from watertap.core.wt_database import Database
import watertap.core.zero_order_properties as prop_ZO
from watertap.core.zero_order_costing import ZeroOrderCosting
from watertap.unit_models.zero_order import (
    FeedZO,
)

idaeslogger.getLogger("idaes.init").setLevel("CRITICAL")
_this_dir = os.path.dirname(os.path.abspath(__file__))


def _guess_and_load_csv_file(unit):
    unit_name = unit._tech_type
    files = list(glob.glob(os.path.join(_this_dir, "data", "*.csv")))

    winners = []
    for fn in files:
        if unit_name in fn:
            winners.append(fn)

    if len(winners) == 0:
        raise RuntimeError(f"Could not file a file for {unit_name}")
    elif len(winners) >= 2:
        raise RuntimeError(
            f"Could not disambiguate {unit_name}, found {len(winners)} candidates: {winners}"
        )

    df = _load_csv_file(winners[0])

    unit = df["unit"]

    assert (unit == unit[0]).all()

    if unit_name not in unit[0]:
        raise RuntimeError(f"Opened {winners[0]}, but not certain unit matches")

    return df


def _load_csv_file(file_name):
    print(f"loading {file_name}")
    return pd.read_csv(file_name)


def _build_flowsheet(unit_model_class, process_subtype, water_source):
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.prop = prop_ZO.WaterParameterBlock(
        default={
            "database": m.db,
            "water_source": water_source,
        }
    )

    # unit model
    m.fs.feed = FeedZO(default={"property_package": m.fs.prop})
    m.fs.unit = unit_model_class(
        default={
            "property_package": m.fs.prop,
            "database": m.db,
            "process_subtype": process_subtype,
        }
    )
    m.fs.product = Product(default={"property_package": m.fs.prop})

    # arcs
    m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.unit.inlet)
    m.fs.s02 = Arc(
        source=m.fs.unit.outlet if hasattr(m.fs.unit, "outlet") else m.fs.unit.treated,
        destination=m.fs.product.inlet,
    )
    TransformationFactory("network.expand_arcs").apply_to(m)

    # costing
    m.fs.costing = ZeroOrderCosting(
        default={
            "case_study_definition": os.path.join(_this_dir, "wt3_test_tea_data.yaml")
        }
    )
    m.fs.unit.costing = UnitModelCostingBlock(
        default={"flowsheet_costing_block": m.fs.costing}
    )
    m.fs.costing.cost_process()
    m.fs.costing.add_electricity_intensity(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol)

    # fix concentration for feed
    m.fs.feed.load_feed_data_from_database(overwrite=True)
    m.fs.feed.flow_vol[0].unfix()

    # load database parameters
    m.fs.unit.load_parameters_from_database(use_default_removal=True)

    assert_units_consistent(m)

    return m


def _initialize_flowsheet(m):
    m.fs.feed.initialize()
    m.fs.unit.initialize()
    m.fs.product.initialize()

    m.fs.costing.initialize()


def _run_flow_in_only(m, df):

    s = get_solver()
    watertap_costing_attributes = {
        "total_capital_cost": [],
        "total_fixed_operating_cost": [],
        "total_operating_cost": [],
        "LCOW": [],
    }
    print("Starting analsys")
    for _, row in df.iterrows():
        if _ % 10 == 0:
            print(f"At {row['flow_in']:.6f} m^3/s")
        m.fs.feed.flow_vol[0].fix(row["flow_in"] * pyunits.m**3 / pyunits.s)
        _initialize_flowsheet(m)
        s.solve(m)
        for att, vals in watertap_costing_attributes.items():
            if att == "LCOW":
                vals.append(value(m.fs.costing.component(att)) * 1e06)
            else:
                vals.append(value(m.fs.costing.component(att)))

    for att, vals in watertap_costing_attributes.items():
        df["WT_" + att] = vals


_WT3_stone = {
    "tci": "WT_total_capital_cost",
    "fixed_op_cost": "WT_total_fixed_operating_cost",
    "annual_op_cost": "WT_total_operating_cost",
    "lcow": "WT_LCOW",
}

_long_name_to_WT3_name = {
    "Total Capital Investment (M\$)": "tci",
    "Unit LCOW (\$/m^3)": "lcow",
    "Annual Fixed Operating Cost (M\$/yr)": "fixed_op_cost",
    "Annual Total Operating Cost (M\$/yr)": "annual_op_cost",
}


def _worst_relative_difference(a, b):
    return (np.abs(a - b) / np.maximum(np.maximum(np.abs(a), np.abs(b)), 1e-06)).max()


def ZeroOrderModel(data):
    if isinstance(data, type):
        return data
    else:
        assert isinstance(data, str)
        import watertap.unit_models.zero_order as zo

        return getattr(zo, data)


class ZeroOrderUnitChecker:

    CONFIG = ConfigDict()

    CONFIG.declare(
        "zero_order_model",
        ConfigValue(
            domain=ZeroOrderModel,
            description="Zero Order Model to validate",
        ),
    ).declare_as_argument()

    CONFIG.declare(
        "process_subtype",
        ConfigValue(
            domain=str,
            default=None,
            description="Process subtype, if needed",
        ),
    ).declare_as_argument()

    CONFIG.declare(
        "water_source",
        ConfigValue(
            domain=str,
            default=None,
            description="Database water_source, if needed",
        ),
    ).declare_as_argument()

    CONFIG.declare(
        "csv_file",
        ConfigValue(
            domain=Path(),
            default=None,
            description="CSV file with validation data",
        ),
    ).declare_as_argument()

    def __init__(self, **options):
        self.config = self.CONFIG()
        self.config.set_value(options)

        self.model = _build_flowsheet(
            self.config.zero_order_model,
            self.config.process_subtype,
            self.config.water_source,
        )

        if self.config.csv_file is None:
            self.comparison_dataframe = _guess_and_load_csv_file(self.model.fs.unit)
        else:
            self.comparison_dataframe = _load_csv_file(self.config.csv_file)

        df = self.comparison_dataframe

        expected_df_columns = (
            "unit",
            "flow_in",
            "flow_out",
            "fci",
            "tci",
            "electricity_intensity",
            "electricity_cost",
            "chem_cost",
            "other_cost",
            "fixed_op_cost",
            "annual_op_cost",
            "lcow",
        )

        if set(df.columns) - set(expected_df_columns):
            raise RuntimeError(
                f"Unexpected columns in csv file: {set(df.columns) - set(expected_df_columns)}"
            )

        if degrees_of_freedom(self.model) != 1:
            raise RuntimeError(
                f"Unexpected number of degrees of freedom {degrees_of_freedom(self.model)}"
            )

        self.worst_difference = None

    def check_unit(self):
        df = self.comparison_dataframe
        msg = f"Checking {self.config.zero_order_model.__name__}"
        if self.config.process_subtype is not None:
            msg += f" with subtype {self.config.process_subtype}"
        print(msg)
        _run_flow_in_only(self.model, df)

        self.worst_difference = max(
            _worst_relative_difference(np.array(df[wt3_key]), np.array(df[wt_k]))
            for wt3_key, wt_k in _WT3_stone.items()
        )
        print(f"Worst relative difference: {self.worst_difference*100:.4f}%")
        return self.worst_difference

    def get_differnces_figure(self):
        assert self.comparison_dataframe is not None

        fig, ax = subplots(1, len(_long_name_to_WT3_name), sharex=True, figsize=(18, 6))

        xaxis_label = "Flow In (m^3/s)"

        df = self.comparison_dataframe

        for idx, (ln, wt3n) in enumerate(_long_name_to_WT3_name.items()):
            ax[idx].plot(df["flow_in"], df[wt3n], label="WT3", color="blue")
            ax[idx].plot(
                df["flow_in"], df[_WT3_stone[wt3n]], label="WT", color="orange"
            )
            ax[idx].set_xlabel(xaxis_label)
            ax[idx].set_ylabel(ln)

        name = f"Unit Model {self.config.zero_order_model.__name__}"
        if self.config.process_subtype:
            name += " ,subtype {self.config.process_subtype}"
        if self.worst_difference:
            name += f"\nWorst Relative Difference: {self.worst_difference*100:.4f}%"
        fig.suptitle(name)
        # only use the last legend
        fig.legend(*ax[idx].get_legend_handles_labels())
        fig.tight_layout()

        return fig


def check_unit(**kwargs):
    from matplotlib.pyplot import show

    checker = ZeroOrderUnitChecker(**kwargs)
    checker.check_unit()
    fig = checker.get_differnces_figure()
    show()

    return checker


if __name__ == "__main__":
    import sys
    import argparse

    parser = argparse.ArgumentParser("unit_validator")
    config = ZeroOrderUnitChecker.CONFIG
    config.initialize_argparse(parser)
    if len(sys.argv) == 1:
        parser.print_help()
        exit()
    args = parser.parse_args(sys.argv[1:])
    config.import_argparse(args)

    checker = check_unit(**config)
