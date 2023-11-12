from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.unit_models.electrodialysis_0D import (
    ElectricalOperationMode,
    Electrodialysis0D,
    PressureDropMethod,
    FrictionFactorMethod,
    HydraulicDiameterMethod,
    LimitingCurrentDensityMethod,
)

from pyomo.environ import (
    ConcreteModel,
    SolverFactory,
    TerminationCondition,
    value,
    Constraint,
    Var,
    Objective,
    Expression,
    assert_optimal_termination,
    log,
)
from idaes.core import (
    FlowsheetBlock,
    EnergyBalanceType,
    MaterialBalanceType,
    MomentumBalanceType,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
import idaes.core.util.scaling as iscale
from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver
import pyomo.environ as pyo

from idaes.core.util.model_diagnostics import DegeneracyHunter


def main():
    solver = get_solver()
    # create model, flowsheet
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    ion_dict = {
        "solute_list": ["Na_+", "Cl_-"],
        "mw_data": {"H2O": 18e-3, "Na_+": 23e-3, "Cl_-": 35.5e-3},
        "elec_mobility_data": {("Liq", "Na_+"): 5.19e-8, ("Liq", "Cl_-"): 7.92e-8},
        "charge": {"Na_+": 1, "Cl_-": -1},
    }
    m.obj = Objective(expr=0)

    # attach prop pack to flowsheet
    m.fs.properties = MCASParameterBlock(**ion_dict)

    m.fs.EDstack = Electrodialysis0D(
        property_package=m.fs.properties,
        operation_mode=ElectricalOperationMode.Constant_Voltage,
        # has_pressure_change=False,
        # limiting_current_density_method=LimitingCurrentDensityMethod.Theoretical,
        has_pressure_change=True,
        has_nonohmic_potential_membrane=False,
        has_Nernst_diffusion_layer=False,
        pressure_drop_method=PressureDropMethod.Darcy_Weisbach,
        # pressure_drop_method=PressureDropMethod.experimental,
        # pressure_drop_method=PressureDropMethod.none,
        # friction_factor_method=FrictionFactorMethod.Gurreri,
        # friction_factor_method=FrictionFactorMethod.Kuroda,
        friction_factor_method=FrictionFactorMethod.fixed,
        hydraulic_diameter_method=HydraulicDiameterMethod.conventional,
        # hydraulic_diameter_method=HydraulicDiameterMethod.spacer_specific_area_known,
        # hydraulic_diameter_method=HydraulicDiameterMethod.fixed,
    )
    assert_units_consistent(m)

    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 0.1, index=("Liq", "H20")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 5e2, index=("Liq", "Na_+")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 5e2, index=("Liq", "Cl_-")
    )

    # iscale.set_scaling_factor(m.fs.EDstack.cell_width, 5)
    # iscale.set_scaling_factor(m.fs.EDstack.cell_length, 5)
    # iscale.set_scaling_factor(m.fs.EDstack.cell_pair_num, 5)
    # iscale.set_scaling_factor(m.fs.EDstack.voltage, 5)

    iscale.calculate_scaling_factors(m.fs.EDstack)

    # specify the model by fixing parameters
    # fix state variables in diluate and concentrate
    # 3 flow_mol_phase_comp for each chamber, 6

    m.fs.EDstack.inlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.4e-1)
    m.fs.EDstack.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-4)
    m.fs.EDstack.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-4)
    m.fs.EDstack.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.4e-1)
    m.fs.EDstack.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-4)
    m.fs.EDstack.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-4)

    # 2 temperature and pressure for each chamber ,4
    m.fs.EDstack.inlet_diluate.pressure[0].fix(201315)
    m.fs.EDstack.inlet_diluate.temperature.fix(298.15)
    m.fs.EDstack.inlet_concentrate.pressure[0].fix(201315)
    m.fs.EDstack.inlet_concentrate.temperature.fix(298.15)

    m.fs.EDstack.voltage.fix(0.5)
    m.fs.EDstack.cell_width.fix(0.1)
    m.fs.EDstack.cell_length.fix(0.79)
    m.fs.EDstack.cell_pair_num.fix(10)

    m.fs.EDstack.current_utilization.fix(1)
    m.fs.EDstack.channel_height.fix(2.7e-4)

    m.fs.EDstack.electrodes_resistance.fix(0)
    m.fs.EDstack.spacer_porosity.fix(0.83)

    m.fs.EDstack.water_trans_number_membrane["cem"].fix(5.8)
    m.fs.EDstack.water_trans_number_membrane["aem"].fix(4.3)
    m.fs.EDstack.water_permeability_membrane["cem"].fix(2.16e-14)
    m.fs.EDstack.water_permeability_membrane["aem"].fix(1.75e-14)
    m.fs.EDstack.membrane_areal_resistance["cem"].fix(1.89e-4)
    m.fs.EDstack.membrane_areal_resistance["aem"].fix(1.77e-4)
    m.fs.EDstack.membrane_thickness["aem"].fix(1.3e-4)
    m.fs.EDstack.membrane_thickness["cem"].fix(1.3e-4)
    m.fs.EDstack.solute_diffusivity_membrane["cem", "Na_+"].fix(1.8e-10)
    m.fs.EDstack.solute_diffusivity_membrane["aem", "Na_+"].fix(1.25e-10)
    m.fs.EDstack.solute_diffusivity_membrane["cem", "Cl_-"].fix(1.8e-10)
    m.fs.EDstack.solute_diffusivity_membrane["aem", "Cl_-"].fix(1.25e-10)
    m.fs.EDstack.ion_trans_number_membrane["cem", "Na_+"].fix(1)
    m.fs.EDstack.ion_trans_number_membrane["aem", "Na_+"].fix(0)
    m.fs.EDstack.ion_trans_number_membrane["cem", "Cl_-"].fix(0)
    m.fs.EDstack.ion_trans_number_membrane["aem", "Cl_-"].fix(1)

    m.fs.EDstack.diffus_mass.fix(1.6e-9) if hasattr(m.fs.EDstack, "diffus_mass") else 0
    #
    # m.fs.EDstack.spacer_specific_area.fix(10700) if hasattr(
    #     m.fs.EDstack, "spacer_specific_area"
    # ) else 0
    # m.display()

    # m.fs.EDstack.hydraulic_diameter.fix(1e-3) if hasattr(
    #     m.fs.EDstack, "hydraulic_diameter"
    # ) else 0
    m.fs.EDstack.friction_factor.fix(20) if hasattr(
        m.fs.EDstack, "friction_factor"
    ) else 0
    # m.fs.EDstack.pressure_drop.fix(1.32527e3)

    m.display()
    m.fs.EDstack.report()

    print("DOF IS", degrees_of_freedom(m))
    solver = get_solver()
    m.fs.EDstack.initialize(optarg=solver.options)
    m.fs.EDstack.report()

    results = solver.solve(m)
    print(results.solver.termination_condition)

    m.fs.EDstack.report()

    m.fs.EDstack.visc_d.pprint()
    m.fs.EDstack.dens_mass.pprint()
    m.fs.EDstack.pressure_drop_total.pprint() if hasattr(
        m.fs.EDstack, "pressure_drop_total"
    ) else print("No it doesn't have pressure_drop total ")

    m.fs.EDstack.pressure_drop.pprint() if hasattr(
        m.fs.EDstack, "pressure_drop"
    ) else print("No it doesn't have pressure_drop ")
    m.fs.EDstack.hydraulic_diameter.pprint() if hasattr(
        m.fs.EDstack, "hydraulic_diameter"
    ) else print("No it doesn't have hydraulic diameter")

    m.fs.EDstack.current_dens_lim_ioa.pprint() if hasattr(
        m.fs.EDstack, "current_dens_lim_ioa"
    ) else print("No it doesn't have current_dens_lim_ioa")
    # m.fs.EDstack.pressure_drop_total.fix(0)

    m.fs.EDstack.diffus_mass.pprint() if hasattr(
        m.fs.EDstack, "diffus_mass"
    ) else print("No it doesn't have diffus_mass")

    m.fs.EDstack.friction_factor.pprint() if hasattr(
        m.fs.EDstack, "friction_factor"
    ) else print("No it doesn't have friction_factor")
    m.fs.EDstack.N_Re.pprint() if hasattr(m.fs.EDstack, "N_Re") else print(
        "No it doesn't have friction_factor"
    )
    m.fs.EDstack.velocity_diluate.pprint()
    m.fs.EDstack.hydraulic_diameter.pprint() if hasattr(
        m.fs.EDstack, "hydraulic_diameter"
    ) else print("No it doesn't have hydraulic_diameter")
    # print("----------------------------------------------------------------------------------------------")
    # print("display after solving")
    # m.display()
    # Specify Ipopt as the solver
    # opt = pyo.SolverFactory('ipopt')
    # # Specifying an iteration limit of 0 allows us to inspect the initial point
    # opt.options['max_iter'] = 0

    # "Solving" the model with an iteration limit of 0 load the initial point and applies
    # any preprocessors (e.g., enforces bounds)
    # opt.solve(m, tee=True)
    # dh = DegeneracyHunter(m, solver=pyo.SolverFactory('cbc'))
    #
    # dh.check_variable_bounds(tol=1.0)
    # dh.check_residuals(tol=0.1)

    # Create Degeneracy Hunter object
    # The Degeneracy Hunter algorithm needs a MILP solver
    # Here we specify CBC, an open source solver

    # opt.options['max_iter'] = 50
    # opt.solve(m, tee=True)
    # dh.check_residuals(tol=1E-14)
    # ds = dh.find_candidate_equations(verbose=True, tee=True)
    # dh.check_variable_bounds(tol=1E-5)
    # dh.check_rank_equality_constraints()

    return m


if __name__ == "__main__":
    m = main()
