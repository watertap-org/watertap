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

"""Zero order nanofiltration model based on specifying solvent flux and solute rejection"""

from pyomo.environ import ConcreteModel, Constraint
from idaes.core import FlowsheetBlock
from watertap.unit_models.nanofiltration_ZO import NanofiltrationZO
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    constraint_scaling_transform,
)
from watertap.examples.flowsheets.full_treatment_train.model_components import (
    property_models,
)
from watertap.examples.flowsheets.full_treatment_train.util import (
    solve_block,
    check_dof,
)


def build_ZONF(m, base="ion"):
    """
    Builds a ZONF model based on specified rejection and solvent flux.
    Requires prop_ion property package.
    """

    if base not in ["ion"]:
        raise ValueError(
            "Unexpected property base {base} for build_ZONF" "".format(base=base)
        )
    prop = property_models.get_prop(m, base=base)

    m.fs.NF = NanofiltrationZO(property_package=prop, has_pressure_change=False)

    # specify
    m.fs.NF.flux_vol_solvent.fix(1.67e-6)
    m.fs.NF.area.fix(500)
    m.fs.NF.properties_permeate[0].pressure.fix(101325)

    m.fs.NF.rejection_phase_comp[0, "Liq", "Na"].fix(0.01)
    m.fs.NF.rejection_phase_comp[0, "Liq", "Ca"].fix(0.79)
    m.fs.NF.rejection_phase_comp[0, "Liq", "Mg"].fix(0.94)
    m.fs.NF.rejection_phase_comp[0, "Liq", "SO4"].fix(0.87)
    m.fs.NF.rejection_phase_comp[
        0, "Liq", "Cl"
    ] = 0.15  # guess, but electroneutrality enforced below
    charge_comp = {"Na": 1, "Ca": 2, "Mg": 2, "SO4": -2, "Cl": -1}
    m.fs.NF.eq_electroneutrality = Constraint(
        expr=0
        == sum(
            charge_comp[j]
            * m.fs.NF.feed_side.properties_out[0].conc_mol_phase_comp["Liq", j]
            for j in charge_comp
        )
    )
    constraint_scaling_transform(m.fs.NF.eq_electroneutrality, 1)


def solve_ZONF(base="ion"):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    property_models.build_prop(m, base=base)
    build_ZONF(m, base=base)
    property_models.specify_feed(m.fs.NF.feed_side.properties_in[0], base="ion")

    check_dof(m)
    calculate_scaling_factors(m)
    solve_block(m)

    m.fs.NF.inlet.display()
    m.fs.NF.permeate.display()
    m.fs.NF.retentate.display()

    return m


if __name__ == "__main__":
    solve_ZONF(base="ion")
