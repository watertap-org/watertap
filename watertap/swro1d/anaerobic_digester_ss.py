from pyomo.environ import (
    ConcreteModel,
    Var,
    Expression,
    Constraint,
    TerminationCondition,
    SolverFactory,
    SolverStatus,
    value,
    units as pyunits,
)
from pyomo.network import Port
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from watertap.core.solvers import get_solver
import watertap.property_models.seawater_prop_pack as props

from watertap.unit_models.anaerobic_digester import AD
from watertap.property_models.anaerobic_digestion.adm1_properties import (
    ADM1ParameterBlock,
)
from watertap.property_models.anaerobic_digestion.adm1_properties_vapor import (
    ADM1_vaporParameterBlock,
)
from watertap.property_models.anaerobic_digestion.adm1_reactions import (
    ADM1ReactionParameterBlock,
)

from watertap.unit_models.tests.unit_test_harness import UnitTestHarness
import idaes.core.util.scaling as iscale

solver = get_solver()

def build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props = ADM1ParameterBlock()
    m.fs.props_vap = ADM1_vaporParameterBlock()
    m.fs.rxn_props = ADM1ReactionParameterBlock(property_package=m.fs.props)

    m.fs.unit = AD(
        liquid_property_package=m.fs.props,
        vapor_property_package=m.fs.props_vap,
        reaction_package=m.fs.rxn_props,
        has_heat_transfer=True,
        has_pressure_change=False,
    )

    # Set the operating conditions
    m.fs.unit.inlet.flow_vol.fix(170 / 24 / 3600)
    m.fs.unit.inlet.temperature.fix(308.15)
    m.fs.unit.inlet.pressure.fix(101325)

    m.fs.unit.inlet.conc_mass_comp[0, "S_su"].fix(0.01)
    m.fs.unit.inlet.conc_mass_comp[0, "S_aa"].fix(0.001)
    m.fs.unit.inlet.conc_mass_comp[0, "S_fa"].fix(0.001)
    m.fs.unit.inlet.conc_mass_comp[0, "S_va"].fix(0.001)
    m.fs.unit.inlet.conc_mass_comp[0, "S_bu"].fix(0.001)
    m.fs.unit.inlet.conc_mass_comp[0, "S_pro"].fix(0.001)
    m.fs.unit.inlet.conc_mass_comp[0, "S_ac"].fix(0.001)
    m.fs.unit.inlet.conc_mass_comp[0, "S_h2"].fix(1e-8)
    m.fs.unit.inlet.conc_mass_comp[0, "S_ch4"].fix(1e-5)
    m.fs.unit.inlet.conc_mass_comp[0, "S_IC"].fix(0.48)
    m.fs.unit.inlet.conc_mass_comp[0, "S_IN"].fix(0.14)
    m.fs.unit.inlet.conc_mass_comp[0, "S_I"].fix(0.02)

    m.fs.unit.inlet.conc_mass_comp[0, "X_c"].fix(2)
    m.fs.unit.inlet.conc_mass_comp[0, "X_ch"].fix(5)
    m.fs.unit.inlet.conc_mass_comp[0, "X_pr"].fix(20)
    m.fs.unit.inlet.conc_mass_comp[0, "X_li"].fix(5)
    m.fs.unit.inlet.conc_mass_comp[0, "X_su"].fix(0.0)
    m.fs.unit.inlet.conc_mass_comp[0, "X_aa"].fix(0.010)
    m.fs.unit.inlet.conc_mass_comp[0, "X_fa"].fix(0.010)
    m.fs.unit.inlet.conc_mass_comp[0, "X_c4"].fix(0.010)
    m.fs.unit.inlet.conc_mass_comp[0, "X_pro"].fix(0.010)
    m.fs.unit.inlet.conc_mass_comp[0, "X_ac"].fix(0.010)
    m.fs.unit.inlet.conc_mass_comp[0, "X_h2"].fix(0.010)
    m.fs.unit.inlet.conc_mass_comp[0, "X_I"].fix(25)

    m.fs.unit.inlet.cations[0].fix(0.04)
    m.fs.unit.inlet.anions[0].fix(0.02)

    m.fs.unit.volume_liquid.fix(3400)
    m.fs.unit.volume_vapor.fix(300)
    m.fs.unit.liquid_outlet.temperature.fix(308.15)

    # Set scaling factors for badly scaled variables
    iscale.set_scaling_factor(
        m.fs.unit.liquid_phase.mass_transfer_term[0, "Liq", "S_h2"], 1e7
    )

    iscale.calculate_scaling_factors(m.fs.unit)

    m.fs.unit.pprint()

    solver = SolverFactory('ipopt')
    solver.options['nlp_scaling_method'] = 'user-scaling'
    solver.options['OF_ma57_automatic_scaling'] = 'yes'
    solver.options['halt_on_ampl_error'] = 'no'
    results = solver.solve(m, tee=True)

    m.fs.unit.report()

    return m

m = build()