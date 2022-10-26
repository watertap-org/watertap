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

from pyomo.environ import (
    ConcreteModel,
    value,
    Var,
    Constraint,
    Component,
    Boolean,
    Expression,
    Objective,
    Param,
    RangeSet,
    Block,
    Suffix,
    TransformationFactory,
    NonNegativeReals,
    units as pyunits,
    assert_optimal_termination,
)
import pyomo
from pyomo.network import Arc
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    ControlVolume0DBlock,
)


from pyomo.common.collections import ComponentSet, ComponentMap

from idaes.core.util import get_solver
from idaes.core.util.exceptions import InitializationError
from watertap.core.util.infeasible import *
from idaes.core.util.model_statistics import degrees_of_freedom, report_statistics
from idaes.core.util.initialization import solve_indexed_blocks, propagate_state
from idaes.generic_models.unit_models import (
    Mixer,
    Separator,
    Product,
    Feed,
    StateJunction,
)
from idaes.generic_models.unit_models.mixer import MomentumMixingType, MixingType
from idaes.generic_models.costing import UnitModelCostingBlock

import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
import copy

from watertap.unit_models.nanofiltration_DSPMDE_0D import (
    NanofiltrationDSPMDE0D,
    MassTransferCoefficient,
    ConcentrationPolarizationType,
)

from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice

from watertap.unit_models.pressure_exchanger import PressureExchanger
from pyomo.util.check_units import assert_units_consistent
from watertap.core.util.initialization import (
    assert_degrees_of_freedom,
    assert_no_degrees_of_freedom,
)
from watertap.costing import WaterTAPCosting
from pyomo.environ import *

# importing analysis specific code and data
import os
import sys
import math

from watertap.property_models.ion_DSPMDE_prop_pack import (
    DSPMDEParameterBlock,
    ActivityCoefficientModel,
    DensityCalculation,
    DSPMDEStateBlock,
)
import numpy as np

# fr
def main():
    m = build()
    feed_mass_frac = {
        "Ca_2+": 4.0034374454637006e-04,
        "HCO3_-": 0.00022696833343821863,
        "SO4_2-": 0.00020497140244420624,
        "Cl_-": 0.0004559124032433401,
        "Na_+": 0.00043333830389924205,
    }
    # feed_mass_frac = {
    #     "Cl_-": 0.0004559124032433401,
    #     "Na_+": 0.00043333830389924205,
    # }
    set_NF_feed(m.fs, feed_mass_frac)
    init(m)
    print(m.fs.nfUnit.area.value)
    solver = get_solver()
    solver.solve(m, tee=True)

    for k in np.linspace(2, 20, 10):
        print("solving for {} bar".format(k))
        m.fs.feed.properties[0].pressure.fix(k * 1e5)
        # init(m)
        res = solver.solve(m, tee=False)
        assert_optimal_termination(res)
        print(
            "current solution: area {} m^2 and water flux {} LMH".format(
                m.fs.nfUnit.area.value,
                value(m.fs.nfUnit.flux_vol_water_avg[0] * 3600 * 1000),
            )
        )


def build():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    default = {
        "solute_list": [
            "Ca_2+",
            "SO4_2-",
            "HCO3_-",
            "Na_+",
            "Cl_-",
        ],
        "diffusivity_data": {
            ("Liq", "Ca_2+"): 9.2e-10,
            ("Liq", "SO4_2-"): 1.06e-9,
            ("Liq", "HCO3_-"): 1.19e-9,
            ("Liq", "Na_+"): 1.33e-9,
            ("Liq", "Cl_-"): 2.03e-9,
        },
        "mw_data": {
            "H2O": 18e-3,
            "Ca_2+": 40e-3,
            "HCO3_-": 61.0168e-3,
            "SO4_2-": 96e-3,
            "Na_+": 23e-3,
            "Cl_-": 35e-3,
        },
        "stokes_radius_data": {
            "Ca_2+": 0.309e-9,
            "HCO3_-": 2.06e-10,
            "SO4_2-": 0.230e-9,
            "Cl_-": 0.121e-9,
            "Na_+": 0.184e-9,
        },
        "charge": {
            "Ca_2+": 2,
            "HCO3_-": -1,
            "SO4_2-": -2,
            "Na_+": 1,
            "Cl_-": -1,
        },
        "activity_coefficient_model": ActivityCoefficientModel.ideal,
        "density_calculation": DensityCalculation.constant,
    }
    # default = {
    #     "solute_list": [
    #         "Na_+",
    #         "Cl_-",
    #     ],
    #     "diffusivity_data": {
    #         ("Liq", "Na_+"): 1.33e-9,
    #         ("Liq", "Cl_-"): 2.03e-9,
    #     },
    #     "mw_data": {
    #         "H2O": 18e-3,
    #         "Na_+": 23e-3,
    #         "Cl_-": 35e-3,
    #     },
    #     "stokes_radius_data": {
    #         "Cl_-": 0.121e-9,
    #         "Na_+": 0.184e-9,
    #     },
    #     "charge": {
    #         "Na_+": 1,
    #         "Cl_-": -1,
    #     },
    #     "activity_coefficient_model": ActivityCoefficientModel.ideal,
    #     "density_calculation": DensityCalculation.constant,
    # }
    m.fs.properties = DSPMDEParameterBlock(default=default)

    m.fs.feed = Feed(default={"property_package": m.fs.properties})

    m.fs.nfUnit = NanofiltrationDSPMDE0D(default={"property_package": m.fs.properties})

    m.fs.feed_to_nf = Arc(source=m.fs.feed.outlet, destination=m.fs.nfUnit.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)
    return m


def init(m):
    solver = get_solver()
    m.fs.feed.initialize(optarg=solver.options)
    m.fs.feed.properties[0].pressure.fix(1.5 * 1e5)
    propagate_state(m.fs.feed_to_nf)
    # m.fs.nfUnit.inlet.pressure[0].fix(2 * 1e5)
    m.fs.nfUnit.recovery_vol_phase.fix(0.1)
    m.fs.nfUnit.area.unfix()
    m.fs.nfUnit.area = 100

    # m.fs.nfUnit.area.setub(1000)
    # m.fs.nfUnit.area.setlb(1)
    # m.fs.nfUnit.width.setub(10000)
    # m.fs.nfUnit.width.setlb(0.1)
    # m.fs.nfUnit.length.setub(10000)
    # m.fs.nfUnit.length.setlb(1)

    m.fs.nfUnit.spacer_porosity.fix(0.85)
    m.fs.nfUnit.channel_height.fix(1e-3)
    m.fs.nfUnit.velocity[0, 0].fix(0.25)
    m.fs.nfUnit.spacer_mixing_efficiency.fix()
    m.fs.nfUnit.spacer_mixing_length.fix()

    m.fs.nfUnit.radius_pore.fix(0.5e-9)
    # m.fs.nfUnit.membrane_thickness_effective.fix(1.33e-6)
    m.fs.nfUnit.membrane_thickness_effective.fix(8.598945196055952e-07)
    m.fs.nfUnit.membrane_charge_density.fix(-680)
    m.fs.nfUnit.dielectric_constant_pore.fix(41.3)
    m.fs.nfUnit.mixed_permeate[0].pressure.fix(101325)
    m.fs.nfUnit.initialize(
        optarg=solver.options,
        automate_rescale=False,
        initialize_guess={
            "solvent_recovery": 0.01,
            "solute_recovery": 0.01,
            "cp_modulus": 1,
        },
    )
    # m.fs.nfUnit.display()


def set_NF_feed(blk, feed_mass_frac):
    mass_flow_in = 1 * pyunits.kg / pyunits.s
    # blk.feed.display()
    for ion, x in feed_mass_frac.items():
        mol_comp_flow = (
            x
            * pyunits.kg
            / pyunits.kg
            * mass_flow_in
            / blk.feed.properties[0].mw_comp[ion]
        )
        blk.feed.properties[0].flow_mol_phase_comp["Liq", ion].fix(mol_comp_flow)
    H2O_mass_frac = 1 - sum(x for x in feed_mass_frac.values())
    H2O_mol_comp_flow = (
        H2O_mass_frac
        * pyunits.kg
        / pyunits.kg
        * mass_flow_in
        / blk.feed.properties[0].mw_comp["H2O"]
    )
    blk.feed.properties[0].flow_mol_phase_comp["Liq", "H2O"].fix(H2O_mol_comp_flow)
    set_NF_feed_scaling(blk)
    blk.feed.properties[0].assert_electroneutrality(
        defined_state=True,
        adjust_by_ion="Cl_-",
        get_property="flow_mol_phase_comp",
    )
    # for index in blk.feed.properties[0].flow_mol_phase_comp:
    #     print(value(blk.feed.properties[0].flow_mol_phase_comp[index]))

    blk.feed.properties[0].temperature.fix(298.15)
    # blk.feed.properties[0].pressure.fix(101325)
    iscale.calculate_scaling_factors(blk)


def calc_scale(value):
    return -1 * math.floor(math.log(value, 10))


def set_NF_feed_scaling(blk):
    _add = 0
    for index in blk.feed.properties[0].flow_mol_phase_comp:
        scale = calc_scale(blk.feed.properties[0].flow_mol_phase_comp[index].value)
        # if "SO4_2-" in index[1]:
        #    _add = 2
        print(f"{index} flow_mol_phase_comp scaling factor = {10**(scale+_add)}")
        blk.properties.set_default_scaling(
            "flow_mol_phase_comp", 10 ** (scale + 0), index=index
        )


if __name__ == "__main__":
    main()
