import matplotlib.pyplot as plt
from pyomo.environ import (
    ConcreteModel,
    Var,
    Expression,
    Constraint,
    TerminationCondition,
    SolverFactory,
    SolverStatus,
    TransformationFactory,
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
import idaes.core.solvers.petsc as petsc  # petsc utilities module
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
    m.fs = FlowsheetBlock(
    dynamic=True,
    time_set=[0, 3],
    time_units=pyunits.s
    )
    # m.fs = FlowsheetBlock(dynamic=False)

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

    for t in m.fs.time:
        # Set the operating conditions
        m.fs.unit.inlet.flow_vol[t].fix(170 / 24 / 3600)
        m.fs.unit.inlet.temperature[t].fix(308.15)
        m.fs.unit.inlet.pressure[t].fix(101325)

        m.fs.unit.inlet.conc_mass_comp[t, "S_su"].fix(0.01)
        m.fs.unit.inlet.conc_mass_comp[t, "S_aa"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[t, "S_fa"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[t, "S_va"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[t, "S_bu"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[t, "S_pro"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[t, "S_ac"].fix(0.001)
        m.fs.unit.inlet.conc_mass_comp[t, "S_h2"].fix(1e-8)
        m.fs.unit.inlet.conc_mass_comp[t, "S_ch4"].fix(1e-5)
        m.fs.unit.inlet.conc_mass_comp[t, "S_IC"].fix(0.48)
        m.fs.unit.inlet.conc_mass_comp[t, "S_IN"].fix(0.14)
        m.fs.unit.inlet.conc_mass_comp[t, "S_I"].fix(0.02)

        m.fs.unit.inlet.conc_mass_comp[t, "X_c"].fix(2)
        m.fs.unit.inlet.conc_mass_comp[t, "X_ch"].fix(5)
        m.fs.unit.inlet.conc_mass_comp[t, "X_pr"].fix(20)
        m.fs.unit.inlet.conc_mass_comp[t, "X_li"].fix(5)
        m.fs.unit.inlet.conc_mass_comp[t, "X_su"].fix(0.0)
        m.fs.unit.inlet.conc_mass_comp[t, "X_aa"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[t, "X_fa"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[t, "X_c4"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[t, "X_pro"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[t, "X_ac"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[t, "X_h2"].fix(0.010)
        m.fs.unit.inlet.conc_mass_comp[t, "X_I"].fix(25)

        m.fs.unit.inlet.cations[t].fix(0.04)
        m.fs.unit.inlet.anions[t].fix(0.02)
    m.fs.unit.liquid_phase.fix_initial_conditions()

    m.fs.unit.volume_liquid.fix(3400)
    m.fs.unit.volume_vapor.fix(300)
    m.fs.unit.liquid_outlet.temperature.fix(308.15)
    
    m.fs.unit.report()

    # Set scaling factors for badly scaled variables
    iscale.set_scaling_factor(
        m.fs.unit.liquid_phase.mass_transfer_term[t, "Liq", "S_h2"], 1e7
    )

    iscale.calculate_scaling_factors(m.fs.unit)

    # solver = SolverFactory('ipopt')
    # # solver.options['max_iter'] = 1000
    # solver.options['nlp_scaling_method'] = 'user-scaling'
    # solver.options['OF_ma57_automatic_scaling'] = 'yes'
    # solver.options['halt_on_ampl_error'] = 'yes'
    # results = solver.solve(m, tee=True)

    TransformationFactory("dae.finite_difference").apply_to(
        m.fs, nfe=1, wrt=m.fs.time, scheme="BACKWARD"
    )

    m.fs.fix_initial_conditions()
    m.fs.unit.initialize()
    m.fs.unit.report()
    # m.fs.unfix_initial_conditions()

    return m

m = build()

result = petsc.petsc_dae_by_time_element(
    m,
    time=m.fs.time,
    keepfiles=True,
    symbolic_solver_labels=True,
    ts_options={
        "--ts_type": "beuler",
        "--ts_dt": 0.01,
        "--ts_rtol": 1e-3,
        # "--ts_adapt_clip":"0.001,3600",
        # "--ksp_monitor":"",
        "--ts_adapt_dt_min": 1e-3,
        "--ts_adapt_dt_max": 3600,
        "--snes_type": "newtontr",
        # "--ts_max_reject": 200,
        # "--snes_monitor":"",
        "--ts_monitor": "",
        "--ts_save_trajectory": 1,
        "--ts_trajectory_type": "visualization",
        "--ts_max_snes_failures": 25,
        # "--show_cl":"",
        "-snes_max_it": 50,
        "-snes_rtol": 0,
        "-snes_stol": 0,
        "-snes_atol": 1e-6,
    },
    skip_initial=False,
    initial_solver="ipopt",
    initial_solver_options={
        "constr_viol_tol": 1e-8,
        "nlp_scaling_method": "user-scaling",
        "linear_solver": "mumps",
        # "OF_ma57_automatic_scaling": "yes",
        "max_iter": 300,
        "tol": 1e-8,
        "halt_on_ampl_error": "no",
    },
)
tj = result.trajectory  # trajectroy data
res = result.results  # solver status list
print(dir(m.fs.unit.vapor_phase))
m.fs.unit.vapor_phase.report(3)
a = plt.plot(tj.time, tj.get_vec(m.fs.unit.vapor_phase.values[3]))
a = plt.ylabel("valve 1 fraction open")
a = plt.xlabel("time (s)")
plt.show()