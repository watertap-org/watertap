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

"""Simple flowsheet with zero order separators for NF and RO examples"""

from pyomo.environ import ConcreteModel, Constraint, ConstraintList, TransformationFactory
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models.translator import Translator
from proteuslib.flowsheets.full_treatment_train.example_models import unit_separator
from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof
from idaes.core.util.scaling import calculate_scaling_factors


def build_flowsheet_NF_salt_basis_example(m):
    """
    Builds a flowsheet with NF and RO connected without bypass. Both the NF and RO are
    modeled as separators with specified split fractions. The flowsheet includes a translator
    block to convert between the NF's multi-salt property package and the RO's seawater package.
    """
    # build flowsheet
    unit_separator.build_NF_salt_separator_example(m)

    unit_separator.build_RO_separator_example(m)
    m.fs.tb_NF_to_RO = Translator(
        default={"inlet_property_package": m.fs.NF_properties,
                 "outlet_property_package": m.fs.RO_properties})

    # add constraints for translator block (tb)
    m.fs.tb_NF_to_RO.eq_H2O_balance = Constraint(
        expr=m.fs.tb_NF_to_RO.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']
        == m.fs.tb_NF_to_RO.outlet.flow_mass_phase_comp[0, 'Liq', 'H2O'])
    m.fs.tb_NF_to_RO.eq_TDS_balance = Constraint(
        expr=sum(m.fs.tb_NF_to_RO.inlet.flow_mass_phase_comp[0, 'Liq', j] for j in ['NaCl', 'CaSO4', 'MgSO4', 'MgCl2'])
        == m.fs.tb_NF_to_RO.outlet.flow_mass_phase_comp[0, 'Liq', 'TDS'])
    m.fs.tb_NF_to_RO.eq_equal_temperature = Constraint(
        expr=m.fs.tb_NF_to_RO.inlet.temperature[0]
        == m.fs.tb_NF_to_RO.outlet.temperature[0])
    m.fs.tb_NF_to_RO.eq_equal_pressure = Constraint(
        expr=m.fs.tb_NF_to_RO.inlet.pressure[0]
        == m.fs.tb_NF_to_RO.outlet.pressure[0])

    # connect models
    m.fs.S1 = Arc(source=m.fs.NF.permeate, destination=m.fs.tb_NF_to_RO.inlet)
    m.fs.S2 = Arc(source=m.fs.tb_NF_to_RO.outlet, destination=m.fs.RO.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling (unit models are already scaled)
    calculate_scaling_factors(m.fs.tb_NF_to_RO)

    # initialize (default values are close enough)

    # unfix RO inlet
    m.fs.RO.inlet.flow_mass_phase_comp.unfix()
    m.fs.RO.inlet.temperature.unfix()
    m.fs.RO.inlet.pressure.unfix()
    check_dof(m)

    return m

def run_flowsheet_NF_salt_basis_example():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    build_flowsheet_NF_salt_basis_example(m)
    solve_with_user_scaling(m)

    m.fs.NF.inlet.display()
    m.fs.NF.retentate.display()
    m.fs.NF.permeate.display()
    m.fs.RO.inlet.display()
    m.fs.RO.retentate.display()
    m.fs.RO.permeate.display()

    return m


def build_flowsheet_NF_ion_basis_example(m):
    """
    Builds a flowsheet with NF and RO connected without bypass. Both the NF and RO are
    modeled as separators with specified split fractions. The flowsheet includes a translator
    block to convert between the NF's multi-ion property package and the RO's seawater package.
    Note: Cl- is used to balance electro-neutrality.
    """
    # build flowsheet
    unit_separator.build_NF_ion_separator_example(m)

    unit_separator.build_RO_separator_example(m)
    m.fs.tb_NF_to_RO = Translator(
        default={"inlet_property_package": m.fs.NF_properties,
                 "outlet_property_package": m.fs.RO_properties})

    # add constraints for translator block (tb)
    m.fs.tb_NF_to_RO.eq_H2O_balance = Constraint(
        expr=m.fs.tb_NF_to_RO.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']
        == m.fs.tb_NF_to_RO.outlet.flow_mass_phase_comp[0, 'Liq', 'H2O'])
    m.fs.tb_NF_to_RO.eq_TDS_balance = Constraint(
        expr=sum(m.fs.tb_NF_to_RO.inlet.flow_mass_phase_comp[0, 'Liq', j] for j in ['Na', 'Ca', 'Mg', 'SO4', 'Cl'])
        == m.fs.tb_NF_to_RO.outlet.flow_mass_phase_comp[0, 'Liq', 'TDS'])
    m.fs.tb_NF_to_RO.eq_equal_temperature = Constraint(
        expr=m.fs.tb_NF_to_RO.inlet.temperature[0]
        == m.fs.tb_NF_to_RO.outlet.temperature[0])
    m.fs.tb_NF_to_RO.eq_equal_pressure = Constraint(
        expr=m.fs.tb_NF_to_RO.inlet.pressure[0]
        == m.fs.tb_NF_to_RO.outlet.pressure[0])

    # connect models
    m.fs.S1 = Arc(source=m.fs.NF.permeate, destination=m.fs.tb_NF_to_RO.inlet)
    m.fs.S2 = Arc(source=m.fs.tb_NF_to_RO.outlet, destination=m.fs.RO.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling (unit models are already scaled)
    calculate_scaling_factors(m.fs.tb_NF_to_RO)

    # initialize (default values are close enough)

    # unfix RO inlet
    m.fs.RO.inlet.flow_mass_phase_comp.unfix()
    m.fs.RO.inlet.temperature.unfix()
    m.fs.RO.inlet.pressure.unfix()
    check_dof(m)

    return m


def run_flowsheet_NF_ion_basis_example():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    build_flowsheet_NF_ion_basis_example(m)
    solve_with_user_scaling(m)

    m.fs.NF.inlet.display()
    m.fs.NF.retentate.display()
    m.fs.NF.permeate.display()
    m.fs.RO.inlet.display()
    m.fs.RO.retentate.display()
    m.fs.RO.permeate.display()

    return m


if __name__ == "__main__":
    # run_flowsheet_NF_salt_basis_example()
    run_flowsheet_NF_ion_basis_example()

