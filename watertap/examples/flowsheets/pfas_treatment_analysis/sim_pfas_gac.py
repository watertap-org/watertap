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

import pytest
import pyomo.environ as pyo
from pyomo.network import Port
from pyomo.util.check_units import assert_units_consistent

import idaes.core.util.scaling as iscale
from idaes.core import (
    FlowsheetBlock,
    EnergyBalanceType,
    MaterialBalanceType,
    MomentumBalanceType,
)
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    badly_scaled_var_generator,
)
from idaes.core import UnitModelCostingBlock

from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
)
from idaes.models.unit_models import Feed, Product
from watertap.unit_models.gac import (
    GAC,
    FilmTransferCoefficientType,
    SurfaceDiffusionCoefficientType,
)
from watertap.costing import WaterTAPCosting

__author__ = "Hunter Barber"


def main():
    # set up solver
    solver = get_solver()

    # run simulation
    m = build()

    # run simulation
    results = solver.solve(m, tee=False)
    assert results.solver.termination_condition == pyo.TerminationCondition.optimal

    m.fs.feed.display()


def build():
    # build models
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(
        solute_list=["PFOS"], mw_data={"H2O": 0.018, "PFOS": 0.50013}
    )
    m.fs.feed = Feed(property_package=m.fs.properties)

    # scaling property variables
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e7, index=("Liq", "PFOS")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e-3, index=("Liq", "H2O")
    )

    # touch properties
    m.fs.feed.properties[0].conc_mass_phase_comp
    m.fs.feed.properties[0].flow_vol_phase["Liq"]

    # propagate scaling factors
    iscale.calculate_scaling_factors(m)

    # feed specifications
    m.fs.feed.properties[0].pressure.fix(101325)  # feed pressure [Pa]
    m.fs.feed.properties[0].temperature.fix(273.15 + 25)  # feed temperature [K]
    # properties (cannot be fixed for initialization routines, must calculate the state variables)
    m.fs.feed.properties.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 170 / 60 / 60,  # 170 m**3/hr or 1 MGD
            ("conc_mass_phase_comp", ("Liq", "PFOS")): 1585.1e-9,  # 1585.1 ng/L
        },
        hold_state=True,  # fixes the calculated component mass flow rates
    )

    return m


if __name__ == "__main__":
    main()
