from pyomo.environ import (
    ConcreteModel,
    Objective,
    assert_optimal_termination,
    value,
    units as pyunits,
    check_optimal_termination,
)
from idaes.core.util.exceptions import InitializationError

from idaes.core import FlowsheetBlock, MomentumBalanceType
import idaes.core.util.scaling as iscale
from idaes.core.util.model_statistics import (
    fixed_variables_set,
    unfixed_variables_set,
    unused_variables_set,
)
from pyomo.util.check_units import assert_units_consistent
from watertap.property_models.ion_DSPMDE_prop_pack import (
    DSPMDEParameterBlock,
    ActivityCoefficientModel,
    DensityCalculation,
)
from watertap.unit_models.nanofiltration_DSPMDE_0D import (
    NanofiltrationDSPMDE0D,
    MassTransferCoefficient,
)
from watertap.core.util.initialization import check_dof, check_solve
from watertap.core.util.infeasible import *
from idaes.core.util import get_solver
from watertap.core.util.infeasible import *
from idaes.core.util.model_statistics import *
import idaes.logger as idaeslog
from idaes.core.util.model_diagnostics import DegeneracyHunter
import pyomo.environ as pyo
from pyomo.util.calc_var_value import calculate_variable_from_constraint
import logging
from pyomo.util.infeasible import *
from copy import deepcopy


if __name__ == "__main__":
    solver = get_solver()

    m = ConcreteModel()

    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = DSPMDEParameterBlock(
        default={
            "solute_list": [
                "Ca_2+",
                "SO4_2-",
                "Mg_2+",
                "Na_+",
                "Cl_-",
                # 'X',
                # 'Y',
            ],
            "diffusivity_data": {
                # this value is smaller than in MIT's paper (0.792e-9) according to NanoFiltran (9.20E-10)
                (
                    "Liq",
                    "Ca_2+",
                ): 9.2e-10,
                ("Liq", "SO4_2-"): 1.06e-9,
                ("Liq", "Mg_2+"): 0.707e-9,
                ("Liq", "Na_+"): 1.33e-9,
                ("Liq", "Cl_-"): 2.03e-9,
                # ("Liq", "X"): 2.03e-9,
                # ("Liq", "Y"): 2.03e-9,
            },
            "mw_data": {
                "H2O": 18e-3,
                "Ca_2+": 40e-3,
                "Mg_2+": 24e-3,
                "SO4_2-": 96e-3,
                "Na_+": 23e-3,
                "Cl_-": 35e-3,
                # "X": 20e-3,
                # "Y": 20e-3,
            },
            "stokes_radius_data": {
                "Ca_2+": 0.309e-9,
                "Mg_2+": 0.347e-9,
                "SO4_2-": 0.231e-9,
                "Cl_-": 0.121e-9,
                "Na_+": 0.184e-9,
                # "X": 0.184e-9,
                # "Y": 0.184e-9,
            },
            "charge": {
                "Ca_2+": 2,
                "Mg_2+": 2,
                "SO4_2-": -2,
                "Na_+": 1,
                "Cl_-": -1,
                # "X": -1,
                # "Y": 1
            },
            "activity_coefficient_model": ActivityCoefficientModel.davies,
            "density_calculation": DensityCalculation.constant,
        }
    )

    m.fs.unit = NanofiltrationDSPMDE0D(
        default={
            "property_package": m.fs.properties,
            "mass_transfer_coefficient": MassTransferCoefficient.spiral_wound,
        }
    )
    b = m.fs.unit
    mass_flow_in = 1 * pyunits.kg / pyunits.s
    feed_mass_frac = {
        "Ca_2+": 382e-6,
        "Mg_2+": 1394e-6,
        "SO4_2-": 2136e-6,
        "Cl_-": 20101.6e-6,
        "Na_+": 11122e-6,
        # 'Cl_-': 0.016924782608695656,
        # 'X': 11122e-6,
        # 'Y': 11122e-6,
    }
    print("Check after build; DOF=", degrees_of_freedom(m))

    # Fix mole flow rates of each ion and water
    for ion, x in feed_mass_frac.items():
        mol_comp_flow = (
            x
            * pyunits.kg
            / pyunits.kg
            * mass_flow_in
            / m.fs.unit.feed_side.properties_in[0].mw_comp[ion]
        )
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", ion].fix(mol_comp_flow)
    H2O_mass_frac = 1 - sum(x for x in feed_mass_frac.values())
    H2O_mol_comp_flow = (
        H2O_mass_frac
        * pyunits.kg
        / pyunits.kg
        * mass_flow_in
        / m.fs.unit.feed_side.properties_in[0].mw_comp["H2O"]
    )
    m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(H2O_mol_comp_flow)
    print("Check after fixing inlet conc; DOF=", degrees_of_freedom(m))
    # m.fs.unit.inlet.flow_mol_phase_comp.display()

    # Use assert electroneutrality method from property model to ensure the ion concentrations provided
    # obey electroneutrality condition
    m.fs.unit.feed_side.properties_in[0].assert_electroneutrality(
        defined_state=True, adjust_by_ion="Cl_-", get_property="mass_frac_phase_comp"
    )
    # m.fs.unit.inlet.flow_mol_phase_comp[0, 'Liq', 'Cl_-'].fix(0.5743318840579711)
    # Double-check that no additional variable is fixed by assert_electroneutrality
    print("Check after electroneutrality assert; DOF=", degrees_of_freedom(m))

    # Fix other inlet state variables
    m.fs.unit.inlet.temperature[0].fix(298.15)
    m.fs.unit.inlet.pressure[0].fix(4e5)

    # Fix the membrane variables that are usually fixed for the DSPM-DE model
    m.fs.unit.radius_pore.fix(5e-10)
    m.fs.unit.membrane_thickness_effective.fix(1.33e-6)
    m.fs.unit.membrane_charge_density.fix(-27)
    m.fs.unit.dielectric_constant_pore.fix(41.3)

    # Fix final permeate pressure to be ~atmospheric
    m.fs.unit.mixed_permeate[0].pressure.fix(101325)

    # Fix additional variables for calculating mass transfer coefficient with spiral wound correlation
    m.fs.unit.spacer_porosity.fix(0.85)
    m.fs.unit.spacer_mixing_efficiency.fix()
    m.fs.unit.spacer_mixing_length.fix()
    m.fs.unit.channel_height.fix(1e-3)

    m.fs.unit.velocity[0, 0].fix(0.25)

    # m.fs.unit.recovery_vol_phase[0, "Liq"].fix(0.50)
    # m.fs.unit.width.fix(1.249)
    m.fs.unit.area.fix(50)

    # m.fs.unit.flux_vol_water_avg.fix(1.67e-06)
    # m.fs.unit.rejection_intrinsic_phase_comp[0, 'Liq', 'Ca_2+'].setlb(.2)

    assert_units_consistent(m)

    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e3, index=("Liq", "Ca_2+")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e3, index=("Liq", "SO4_2-")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e2, index=("Liq", "Mg_2+")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e2, index=("Liq", "Cl_-")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e2, index=("Liq", "Na_+")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e0, index=("Liq", "H2O")
    )
    iscale.calculate_scaling_factors(m.fs.unit)

    print(
        "---------------- After calculate_scaling_factors---------------------------------------"
    )
    [print(i[0], i[1]) for i in iscale.badly_scaled_var_generator(m)]

    # deactivate concentration polarization and interface electroneutrality
    # m.fs.unit.eq_electroneutrality_interface.deactivate()
    # m.fs.unit.eq_solute_flux_concentration_polarization.deactivate()

    # deactivate feed electroneutrality constraint
    m.fs.unit.feed_side.eq_electroneutrality_feed.deactivate()
    # deactivate NO concentration polarization
    m.fs.unit.eq_no_concentration_polarization.deactivate()
    try:
        m.fs.unit.initialize(
            automate_rescale=True,
            optarg={
                "halt_on_ampl_error": "yes",
            },
            outlvl=idaeslog.DEBUG,
        )
    except InitializationError:
        pass

    # Use of Degeneracy Hunter for troubleshooting model.
    # m.fs.dummy_objective = Objective(expr=0)
    # solver.options["max_iter"] = 0
    # solver.solve(m, tee=True)
    # dh = DegeneracyHunter(m, solver=pyo.SolverFactory("cbc"))
    # dh.check_residuals(tol=0.1)
    # assert False
    # ##############################################################################################################
    # dh.check_variable_bounds(tol=1e-3)

    #
    # # assert False
    #
    # # for con in m.fs.unit.component_data_objects(Constraint, descend_into=False):
    # #     con.deactivate()
    # # print("DOF with all cons deactivated:", degrees_of_freedom(m.fs.unit))
    # # for con in m.fs.unit.component_data_objects(Constraint, descend_into=False):
    # #     con.activate()

    # solver.options['max_iter'] = 0

    solver.options["max_iter"] = 2000
    results = solver.solve(m, tee=True)
    if check_optimal_termination(results):
        b.report()
        print("SUCCESSFUL FIRST SOLVE!!!!!!!!!!!")
        # Check that electroneutrality is satisfied for feed outlet and mixed permeate- constraints that
        # are deactivated because they lead to failed solve
        # b.feed_side.properties_out[0].assert_electroneutrality(defined_state=False, tee=True, solve=False)
        # b.mixed_permeate[0].assert_electroneutrality(defined_state=False, tee=True, solve=False)
    else:
        print("FIRST SOLVE FAILED..............")
        # print("Badly scaled vars after FIRST failed solve:")
        # [print(i[0], i[1]) for i in iscale.badly_scaled_var_generator(m)]
        # print(
        #     f"Number of badly scaled vars = {len(list(iscale.badly_scaled_var_generator(m)))}"
        # )
        # m.fs.unit._automate_rescale_variables(rescale_factor=1)
        #
        # results = solver.solve(m, tee=True)
        # if check_optimal_termination(results):
        #     b.report()
        #     print("SUCCESSFUL SECOND SOLVE!!!!!!!!!!!")
        #     # Check that electroneutrality is satisfied for feed outlet and mixed permeate- constraints that
        #     # are deactivated because they lead to failed solve
        #     b.feed_side.properties_out[0].assert_electroneutrality(
        #         defined_state=False, tee=True, solve=False
        #     )
        #     b.mixed_permeate[0].assert_electroneutrality(
        #         defined_state=False, tee=True, solve=False
        #     )
        # else:
        #     print("Badly scaled vars after SECOND failed solve:")
        #     [print(i[0], i[1]) for i in iscale.badly_scaled_var_generator(m)]
        #     print(
        #         f"Number of badly scaled vars = {len(list(iscale.badly_scaled_var_generator(m)))}"
        #     )
        #     print("SECOND SOLVE FAILED...............................")
    #
    # More Degeneracy Hunter troubleshooting
    # dh.check_residuals(tol=1e-10)
    # dh.check_variable_bounds(tol=1e-2, relative=True)
    # try:
    #     dh.check_rank_equality_constraints()
    # except:
    #     ds2 = dh.find_candidate_equations(verbose=True, tee=True)
    #     ids = dh.find_irreducible_degenerate_sets(verbose=True, tee=False)

    b.rejection_intrinsic_phase_comp.display()
    b.rejection_observed_phase_comp.display()
    print(f"Degrees of freedom ={degrees_of_freedom(m)} ")
    print_infeasible_constraints(m)
