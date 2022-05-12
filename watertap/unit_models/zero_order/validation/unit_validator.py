import os
import glob
import pandas as pd

from pyomo.environ import (
    ConcreteModel,
    value,
    TransformationFactory,
    units as pyunits,
)
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent

from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.generic_models.unit_models import Product

from idaes.core import FlowsheetBlock
from idaes.generic_models.costing import UnitModelCostingBlock

from watertap.core.wt_database import Database
import watertap.core.zero_order_properties as prop_ZO
from watertap.core.zero_order_costing import ZeroOrderCosting
from watertap.unit_models.zero_order import (
    FeedZO,
)

_this_dir = os.path.dirname(os.path.abspath(__file__))

def _guess_and_load_csv_file(unit):
    unit_name = unit._tech_type
    files = list(glob.glob(os.path.join(_this_dir, 'data', '*.csv')))

    winners = []
    for fn in files:
        if unit_name in fn:
            winners.append(fn)

    if len(winners) == 0:
        raise RuntimeError(f"Could not file a file for {unit_name}")
    elif len(winners) >= 2:
        raise RuntimeError(f"Could not disambiguate {unit_name}, found {len(winners)} candidates: {winners}")

    df = _load_csv_file(winners[0])

    unit = df['unit']

    assert (unit == unit[0]).all()

    if unit_name not in unit[0]:
        raise RuntimeError(f"Opened {winners[0]}, but not certain unit matches")

    return df

def _load_csv_file(file_name):
    print(f"loading {file_name}")
    return pd.read_csv(file_name)

def _build_flowsheet(unit_model_class, process_subtype):
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.prop = prop_ZO.WaterParameterBlock(
        default={"solute_list": ["tds", "tss", "toc"]}
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
    m.fs.s02 = Arc(source=m.fs.unit.outlet if hasattr(m.fs.unit, "outlet") else m.fs.unit.treated,
                    destination=m.fs.product.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # costing
    m.fs.costing = ZeroOrderCosting()
    m.fs.unit.costing = UnitModelCostingBlock(default={"flowsheet_costing_block":m.fs.costing})
    m.fs.costing.cost_process()
    m.fs.costing.add_electricity_intensity(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol)

    # fix concentration for feed
    conc_mass_tds = 0.63 * pyunits.kg / pyunits.m**3
    conc_mass_tss = 0.006525 * pyunits.kg / pyunits.m**3
    conc_mass_toc = 0.004 * pyunits.kg / pyunits.m**3
    m.fs.feed.conc_mass_comp[0, "tds"].fix(conc_mass_tds)
    m.fs.feed.conc_mass_comp[0, "tss"].fix(conc_mass_tss)
    m.fs.feed.conc_mass_comp[0, "toc"].fix(conc_mass_toc)

    # load database parameters
    m.fs.unit.load_parameters_from_database(use_default_removal=True)

    return m


def check_unit(unit_model_class, process_subtype=None, csv_file_name=None):

    m = _build_flowsheet(unit_model_class, process_subtype)

    if csv_file_name is None:
        df = _guess_and_load_csv_file(m.fs.unit)
    else:
        df = _load_csv_file(csv_file_name)

    expected_df_columns = ('unit', 'flow_in', 'flow_out', 'fci', 'tci', 'electricity_intensity',
            'electricity_cost', 'chem_cost', 'other_cost', 'fixed_op_cost', 'annual_op_cost', 'lcow')

    if (set(df.columns) - set(expected_df_columns)):
        raise RuntimeError(f"Unexpected columns in csv file: {set(df.columns) - set(expected_df_columns)}")

    if degrees_of_freedom(m) != 1:
        raise RuntimeError(f"Unexpected number of degrees of freedom {degrees_of_freedom(m)}")

    return m, df
