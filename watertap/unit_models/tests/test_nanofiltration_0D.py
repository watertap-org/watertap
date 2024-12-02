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

import pytest

from pyomo.environ import ConcreteModel, TransformationFactory, units as pyunits, value

from idaes.core import FlowsheetBlock
from idaes.core.scaling import set_scaling_factor
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.solvers import get_solver

from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
)
from watertap.unit_models.nanofiltration0D import (
    Nanofiltration0D,
)


@pytest.mark.component
def test_nf0d():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(
        solute_list=["Ca_2+", "SO4_2-", "Mg_2+", "Na_+", "Cl_-"],
        diffusivity_data={
            ("Liq", "Ca_2+"): 9.2e-10,
            ("Liq", "SO4_2-"): 1.06e-09,
            ("Liq", "Mg_2+"): 7.06e-10,
            ("Liq", "Na_+"): 1.33e-09,
            ("Liq", "Cl_-"): 2.03e-09,
        },
        mw_data={
            "H2O": 0.018,
            "Ca_2+": 0.04,
            "Mg_2+": 0.024,
            "SO4_2-": 0.096,
            "Na_+": 0.023,
            "Cl_-": 0.035,
        },
        stokes_radius_data={
            "Ca_2+": 3.09e-10,
            "Mg_2+": 3.47e-10,
            "SO4_2-": 2.3e-10,
            "Cl_-": 1.21e-10,
            "Na_+": 1.84e-10,
        },
        charge={"Ca_2+": 2, "Mg_2+": 2, "SO4_2-": -2, "Na_+": 1, "Cl_-": -1},
    )

    m.fs.unit = Nanofiltration0D(property_package=m.fs.properties)

    # Fix other inlet state variables
    m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(53.6036)
    m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Ca_2+"].fix(0.00955)
    m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Mg_2+"].fix(0.5808)
    m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "SO4_2-"].fix(0.02225)
    m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(1.61977)
    m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(0.48357)
    m.fs.unit.inlet.temperature[0].fix(298.15)
    m.fs.unit.inlet.pressure[0].fix(4e5)

    m.fs.unit.solvent_recovery.fix()
    m.fs.unit.solute_recovery.fix()
    m.fs.unit.permeate.pressure[0].fix(1e5)

    # Scale model
    set_scaling_factor(m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"], 1e-1)
    set_scaling_factor(m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Ca_2+"], 1e3)
    set_scaling_factor(m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Mg_2+"], 10)
    set_scaling_factor(m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "SO4_2-"], 1e2)
    set_scaling_factor(m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Cl_-"], 1)
    set_scaling_factor(m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "Na_+"], 10)

    scaler = m.fs.unit.default_scaler()
    scaler.scale_model(m.fs.unit)

    initializer = m.fs.unit.default_initializer()
    try:
        initializer.initialize(m.fs.unit)
    except:
        pass

    solver = get_solver("ipopt_v2", writer_config={"scale_model": True, "linear_presolve": True})
    solver.solve(m)

    sm = TransformationFactory("core.scale_model").create_using(m, rename=False)

    dt = DiagnosticsToolbox(sm.fs.unit)
    dt.report_structural_issues()
    dt.report_numerical_issues()

    m.fs.unit.inlet.display()
    m.fs.unit.retentate.display()
    m.fs.unit.permeate.display()

    assert False