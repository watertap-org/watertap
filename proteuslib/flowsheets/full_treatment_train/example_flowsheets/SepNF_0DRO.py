###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################

"""Flowsheets with a separator NF and 0D RO examples"""

from pyomo.environ import ConcreteModel, Constraint, TransformationFactory
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import Separator, Mixer
from idaes.generic_models.unit_models.separator import SplittingType, EnergySplittingType
from idaes.core.util.initialization import propagate_state, fix_state_vars, revert_state_vars
from proteuslib.unit_models.pump_isothermal import Pump
from proteuslib.flowsheets.full_treatment_train.example_models import unit_separator, unit_0DRO, property_models
from proteuslib.flowsheets.full_treatment_train.example_flowsheets import translator_block
from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof
from idaes.core.util.scaling import calculate_scaling_factors, set_scaling_factor


def build_flowsheet_simple_example(m):
    """
    Build a flowsheet with NF and RO including bypass. Both the NF and RO are
    modeled as separators with specified split fractions. The NF uses prop_salt and
    the RO uses prop_TDS.
    """
    # build flowsheet
    property_models.build_prop(m, base='salt')
    property_models.build_prop(m, base='TDS')
    unit_separator.build_SepNF(m, base='salt')
    unit_0DRO.build_RO(m, level='simple')
    m.fs.splitter = Separator(default={
        "property_package": m.fs.prop_salt,
        "outlet_list": ['pretreatment', 'bypass'],
        "split_basis": SplittingType.totalFlow,
        "energy_split_basis": EnergySplittingType.equal_temperature})
    m.fs.mixer = Mixer(default={
        "property_package": m.fs.prop_salt,
        "inlet_list": ['pretreatment', 'bypass']})
    m.fs.pump = Pump(default={'property_package': m.fs.prop_TDS})

    translator_block.build_tb_salt_to_TDS(m)

    # connect models
    m.fs.s01 = Arc(source=m.fs.splitter.bypass, destination=m.fs.mixer.bypass)
    m.fs.s02 = Arc(source=m.fs.splitter.pretreatment, destination=m.fs.NF.inlet)
    m.fs.s03 = Arc(source=m.fs.NF.permeate, destination=m.fs.mixer.pretreatment)
    m.fs.s04 = Arc(source=m.fs.mixer.outlet, destination=m.fs.tb_salt_to_TDS.inlet)
    m.fs.s05 = Arc(source=m.fs.tb_salt_to_TDS.outlet, destination=m.fs.pump.inlet)
    m.fs.s06 = Arc(source=m.fs.pump.outlet, destination=m.fs.RO.inlet)

    # specify (NF and RO models are already specified, DOF: mixer 0, splitter 1, pump 2)
    # splitter
    m.fs.splitter.split_fraction[0, 'bypass'].fix(0.25)
    # pump
    m.fs.pump.efficiency_pump.fix(0.80)
    m.fs.pump.control_volume.properties_out[0].pressure.fix(50e5)

    # scaling (NF and RO models and translator are already scaled)
    calculate_scaling_factors(m.fs.splitter)
    calculate_scaling_factors(m.fs.mixer)
    set_scaling_factor(m.fs.pump.control_volume.work, 1e-3)
    calculate_scaling_factors(m.fs.pump)

    # initialize
    m.fs.splitter.initialize(optarg={'nlp_scaling_method': 'user-scaling'})
    propagate_state(m.fs.s01)
    propagate_state(m.fs.s02)
    # m.fs.NF.initialize(optarg={'nlp_scaling_method': 'user-scaling'})  # IDAES error TODO: discuss with Andrew
    propagate_state(m.fs.s03)
    m.fs.mixer.initialize(optarg={'nlp_scaling_method': 'user-scaling'})
    propagate_state(m.fs.s04)
    m.fs.tb_salt_to_TDS.properties_in[0].mass_frac_phase_comp  # touch a property to create constraint in stateblock
    m.fs.tb_salt_to_TDS.properties_out[0].mass_frac_phase_comp  # touch a property to create constraint in stateblock
    m.fs.tb_salt_to_TDS.initialize(optarg={'nlp_scaling_method': 'user-scaling'})
    propagate_state(m.fs.s05)
    m.fs.pump.initialize(optarg={'nlp_scaling_method': 'user-scaling'})
    propagate_state(m.fs.s06)
    m.fs.RO.initialize(optarg={'nlp_scaling_method': 'user-scaling'})

    return m


def run_flowsheet_simple_example():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    build_flowsheet_simple_example(m)
    property_models.specify_feed(m.fs.splitter.mixed_state[0], base='salt')

    TransformationFactory("network.expand_arcs").apply_to(m)
    check_dof(m)
    solve_with_user_scaling(m)

    m.fs.splitter.report()
    m.fs.NF.inlet.display()
    m.fs.NF.retentate.display()
    m.fs.NF.permeate.display()
    m.fs.mixer.report()
    m.fs.RO.inlet.display()
    m.fs.RO.retentate.display()
    m.fs.RO.permeate.display()

    return m


if __name__ == "__main__":
    run_flowsheet_simple_example()

