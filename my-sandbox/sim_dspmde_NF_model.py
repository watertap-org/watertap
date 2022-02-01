from pyomo.environ import ConcreteModel, assert_optimal_termination, value, units as pyunits, check_optimal_termination
from idaes.core import FlowsheetBlock, MomentumBalanceType
import idaes.core.util.scaling as iscale
from idaes.core.util.model_statistics import fixed_variables_set, unfixed_variables_set, unused_variables_set
from pyomo.util.check_units import assert_units_consistent
from watertap.property_models.ion_DSPMDE_prop_pack import DSPMDEParameterBlock, ActivityCoefficientModel, DensityCalculation
from watertap.unit_models.nanofiltration_DSPMDE_0D import NanofiltrationDSPMDE0D
from watertap.core.util.initialization import check_dof, check_solve
from pyomo.util.infeasible import *
from idaes.core.util import get_solver
from watertap.core.util.infeasible import *
from idaes.core.util.model_statistics import *
import idaes.logger as idaeslog

import logging
from pyomo.util.infeasible import *


if __name__ == '__main__':
    solver = get_solver()

    m = ConcreteModel()



    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = DSPMDEParameterBlock(default={
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
            ("Liq", "Ca_2+"): 0.792e-9,
            ("Liq", "SO4_2-"): 1.06e-9,
            ("Liq", "Mg_2+"): 0.706e-9,
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
            "SO4_2-": 0.230e-9,
            "Cl_-": 0.121e-9,
            "Na_+": 0.184e-9,
            # "X": 0.184e-9,
            # "Y": 0.184e-9,
        },
        "density_data": {
            "H2O": 1000,
            "Ca_2+": 1550,
            "Mg_2+": 1738,
            "SO4_2-": 2553,
            "Na_+": 968,
            "Cl_-": 3214,
            # "X": 1000,
            # "Y": 1000

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
        "activity_coefficient_model": ActivityCoefficientModel.ideal,
        "density_calculation": DensityCalculation.seawater
    })

    m.fs.unit = NanofiltrationDSPMDE0D(default={"property_package": m.fs.properties})
    mass_flow_in = 1 * pyunits.kg / pyunits.s
    feed_mass_frac = {
                      'Ca_2+': 382e-6,
                      'Mg_2+': 1394e-6,
                      'SO4_2-': 2136e-6,
                      'Cl_-': 20101.6e-6,
                      'Na_+': 11122e-6,
                      # 'Cl_-': 0.016924782608695656,
                      # 'X': 11122e-6,
                      # 'Y': 11122e-6,

                      }
    check_dof(m, fail_flag=False)

    # Fix mole flow rates of each ion and water
    for ion, x in feed_mass_frac.items():
        mol_comp_flow = x * pyunits.kg / pyunits.kg * mass_flow_in / m.fs.unit.feed_side.properties_in[0].mw_comp[ion]
        m.fs.unit.inlet.flow_mol_phase_comp[0, 'Liq', ion].fix(mol_comp_flow)
    H2O_mass_frac = 1 - sum(x for x in feed_mass_frac.values())
    H2O_mol_comp_flow = H2O_mass_frac * pyunits.kg / pyunits.kg * mass_flow_in / \
                        m.fs.unit.feed_side.properties_in[0].mw_comp['H2O']
    m.fs.unit.inlet.flow_mol_phase_comp[0, 'Liq', 'H2O'].fix(H2O_mol_comp_flow)

    # Use assert electroneutrality method from property model to ensure the ion concentrations provided
    # obey electroneutrality condition
    m.fs.unit.feed_side.properties_in[0].assert_electroneutrality(tol=1e-6)

    # Fix other inlet state variables
    m.fs.unit.inlet.temperature[0].fix(298.15)
    m.fs.unit.inlet.pressure[0].fix(4e5)
    # check_dof(m, fail_flag=False)

    # Fix the membrane variables that are usually fixed for the DSPM-DE model
    m.fs.unit.radius_pore.fix(0.5e-9)
    m.fs.unit.membrane_thickness_effective.fix()
    m.fs.unit.membrane_charge_density.fix(-27)
    m.fs.unit.dielectric_constant_pore.fix(41.3)

    # Fix final permeate pressure to be ~atmospheric
    m.fs.unit.mixed_permeate[0].pressure.fix(101325)
    # m.fs.unit.permeate_side[0,:].pressure.fix(101325)
    # check_dof(m, fail_flag=False)

    # Membrane area and recovery haven't typically been included in the literature for DSPM-DE,
    # but since we include them in our model to simulate/optimize at process level, I chose to fix those here
    m.fs.unit.area.fix(40)
    # m.fs.unit.recovery_vol_phase[0, 'Liq'].fix(0.8)

    #Other things I tried fixing when experimenting with DOF reduction
    # m.fs.unit.flux_mol_phase_comp[0, 0, 'Liq', 'H2O'].fix(5e-5)
    # m.fs.unit.flux_vol_water[0, :].set_value(1e-3 *pyunits.m/pyunits.s)
    # m.fs.unit.flux_vol_water_avg[0].set_value(1e-3 *pyunits.m/pyunits.s)

    # fix the mass transfer coefficient, which is just a var without the associated constraint for now.
    # Planning to add in the associated constraint(s) for this once the model solves at this level
    for i in range(2):
        for j in m.fs.properties.solute_set:
            m.fs.unit.Kf_comp[0, i, j].fix()

    # check_dof(m, fail_flag=True)

    # Set electroneutrality tolerance to 0 (value used in equality constraints for electroneutrality in unit model)
    m.fs.unit.tol_electroneutrality = 0

    print ('----------------BEFORE SCALING---------------------------------------')
    # for ion, i in m.fs.unit.feed_side.properties_in[0].radius_stokes_comp.items():
    #     print(f"{ion}stokes radius:", value(i))
    #     print(value(b.feed_side.properties_in[0].radius_stokes_comp[ion]))
    #     print(f"{ion}lambda:", value(i / m.fs.unit.radius_pore))
    #     print(value(b.lambda_comp[0, ion]))
    #     print(value(b.radius_pore))

    # for c in m.fs.unit.component_data_objects(Constraint, active=True):
    #     print(c)

    # b.eq_solute_flux_concentration_polarization.deactivate()
    # b.Kf_comp.deactivate()
    # assert False
    # solute_flux_concentration_polarization_eq

    # m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1e1, index=('Liq', 'H2O'))
    #
    # for ion in m.fs.properties.solute_set:
    #     if j != 'Na_+' and j != 'Cl_-':
    #         m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1e4, index=('Liq', j))
    #         iscale.set_scaling_factor(m.fs.unit.feed_side.properties_in[0].conc_mol_phase_comp['Liq', j], 1e5)
    #         iscale.set_scaling_factor(m.fs.unit.feed_side.properties_in[0].conc_mol_phase_comp['Liq', j], 1e5)
    #
    #     else:
    #         m.fs.properties.set_default_scaling('flow_mol_phase_comp', 1e2, index=('Liq', j))
    # iscale.set_scaling_factor(m.fs.unit.feed_side.properties_in[0].conc_mol_phase_comp['Liq', 'Ca_2+'], 1e9)
    # iscale.set_scaling_factor(m.fs.unit.feed_side.properties_out[0].conc_mol_phase_comp['Liq', 'Ca_2+'], 1e9)
    iscale.calculate_scaling_factors(m.fs.unit)
    print('---------------- After scaling---------------------------------------')
    [print(i[0], i[1]) for i in iscale.badly_scaled_var_generator(m)]
    m.fs.unit._automate_rescale_variables()
    # iscale.set_scaling_factor(m.fs.unit.feed_side.properties_in[0].conc_mol_phase_comp['Liq', 'Ca_2+'], 1e3)
    # iscale.set_scaling_factor(m.fs.unit.feed_side.properties_in[0].conc_mol_phase_comp['Liq', 'Ca_2+'], 1e3)
    # iscale.set_scaling_factor(m.fs.unit.feed_side.properties_in[0].mole_frac_phase_comp['Liq', 'Ca_2+'], 1e3)
    # iscale.set_scaling_factor(m.fs.unit.feed_side.properties_in[0].mole_frac_phase_comp['Liq', 'SO4_2-'], 1e3)

    iscale.calculate_scaling_factors(m.fs.unit)
    print('---------------- After scaling again---------------------------------------')
    print('List of badly scaled vars:\n')
    [print(i[0], i[1]) for i in iscale.badly_scaled_var_generator(m)]
    print('End of list of badly scaled vars\n')


    # checking state block
    assert_units_consistent(m)

    # check dof = 0
    check_dof(m, fail_flag=False)

    # iscale.calculate_scaling_factors(m.fs.unit)

    # m.fs.unit.eq_electroneutrality_feed.deactivate()
    # m.fs.unit.eq_equal_flow_vol_pore_permeate.deactivate()
    # m.fs.unit.eq_equal_flow_vol_permeate.deactivate()
    # initialize

    m.fs.unit.initialize()

    print ('---------------- AFTER INITIALIZATION---------------------------------------')
    [print(i[0],i[1]) for i in iscale.badly_scaled_var_generator(m)]
    print('\nNUMBER OF badly scaled variables:', len(list(iscale.badly_scaled_var_generator(m))))
    #
    m.fs.unit._automate_rescale_variables(rescale_factor=1)

    # print('---------------- AFTER AUTOMATE RESCALE---------------------------------------')
    # [print(i[0], i[1]) for i in iscale.badly_scaled_var_generator(m)]
    # print('\nNUMBER OF badly scaled variables:', len(list(iscale.badly_scaled_var_generator(m))))

    # b.display()

    # for ion, i in m.fs.unit.feed_side.properties_in[0].radius_stokes_comp.items():
    #     print(f"{ion}stokes radius:", value(i))
    #     print(value(b.feed_side.properties_in[0].radius_stokes_comp[ion]))
    #     print(f"{ion}lambda:", value(i / m.fs.unit.radius_pore))
    #     print(value(b.lambda_comp[0, ion]))
    #     print(value(b.radius_pore))

    # b.display()
    #
    # # # check solve
    # solver.options['max_iter'] = 108
    # solver.options['constr_viol_tol'] = 1e-3
    # log_infeasible_bounds(m.fs.unit)

    # m.fs.unit.area.unfix()
    m.fs.unit.recovery_vol_phase[0, 'Liq'].fix(0.5)
    for con in m.fs.unit.component_data_objects(Constraint, descend_into=False):
        con.deactivate()

    m.fs.unit.eq_water_flux.activate()
    m.fs.unit.eq_solute_solvent_flux.activate()
    m.fs.unit.eq_solute_flux_concentration_polarization.activate()
    m.fs.unit.eq_permeate_isothermal.activate()
    m.fs.unit.eq_permeate_isothermal_mixed.activate()
    m.fs.unit.eq_pressure_permeate_io.activate()
    m.fs.unit.eq_mass_transfer_feed.activate() # when area was unfixed, solve fails with this constraint active
    m.fs.unit.eq_permeate_production.activate()



    m.fs.unit.eq_solute_flux_pore_domain.activate()
    m.fs.unit.eq_recovery_mol_phase_comp.activate() # model solves and gets ~67%!
    m.fs.unit.eq_rejection_phase_comp.activate()
    m.fs.unit.eq_pore_isothermal.activate()
    m.fs.unit.feed_side.eq_feed_interface_isothermal.activate() #- extends number of iterations
    m.fs.unit.feed_side.eq_feed_isothermal.activate() # - extends number of iterations

    m.fs.unit.eq_recovery_vol_phase.activate() # model solves in 27 iterations (instead of 25), but doesn't match recovery_mol_phase_comp
    m.fs.unit.eq_equal_flow_vol_pore_permeate.activate() # brought in after eq_permeate_isothermal - flow rate values were 1e5 but solve worked
    m.fs.unit.eq_equal_flow_vol_permeate.activate() #- also brought flowrate to 1e5




    # Constraints that could mess up the solve



    # results = solver.solve(m, tee=True)
    # if check_optimal_termination(results):
    #     print('SUCCESS!!!!!!!!!!!')


    # Activate later



    # for con in m.fs.unit.component_data_objects(Constraint):
    #     con.activate()
    b = m.fs.unit

    results = solver.solve(m, tee=True)
    if check_optimal_termination(results):
        b.report()
        print('SUCCESS WITH FIRST SOLVE!!!!!!!!!!!')

    results = solver.solve(m, tee=True)
    if check_optimal_termination(results):
        b.report()
        print('SUCCESS WITH SECOND SOLVE!!!!!!!!!!!')
    # Constraints activated that messed up the solve
    m.fs.unit.eq_interfacial_partitioning_feed.activate()
    m.fs.unit.eq_interfacial_partitioning_permeate.activate()
    m.fs.unit.eq_electroneutrality_mixed_permeate.activate() #- extends number of iterations
    m.fs.unit.eq_electroneutrality_interface.activate()
    m.fs.unit.eq_electroneutrality_pore.activate()
    m.fs.unit.eq_electroneutrality_permeate.activate()
    m.fs.unit.eq_electroneutrality_feed.activate() #- extends number of iterations

    # m.fs.unit.eq_recovery_vol_phase.activate() # model solves in 27 iterations (instead of 25), but doesn't match recovery_mol_phase_comp
    # m.fs.unit.eq_equal_flow_vol_pore_permeate.activate() # brought in after eq_permeate_isothermal - flow rate values were 1e5 but solve worked
    # m.fs.unit.eq_equal_flow_vol_permeate.activate() #- also brought flowrate to 1e5
    # results = solver.solve(m, tee=True)
    # if check_optimal_termination(results):
    #     print('SUCCESS!!!!!!!!!!!')

    logging.basicConfig(filename='infeasible6.log', level=logging.INFO)
    log_infeasible_constraints(m, log_expression=True, log_variables=True)
    log_infeasible_bounds(m)
    # else:
    #     for var, sv in iscale.badly_scaled_var_generator(m):
    #         print(var, sv)
    #         m.fs.unit._automate_rescale_variables(rescale_factor=1)
    #         results = solver.solve(m, tee=True)

    #     b.display()
    #     m.fs.unit.report()
        # m.fs.unit.eq_interfacial_partitioning_feed.display()

    # print('---------------------------------')
    #
    # for ion, i in m.fs.unit.feed_side.properties_in[0].radius_stokes_comp.items():
    #     print(f"Back of envelope {ion}lambda:",
    #           value(b.feed_side.properties_in[0].radius_stokes_comp[ion] / m.fs.unit.radius_pore))
    #     print(f"model {ion} lambda result=", value(b.lambda_comp[0, ion]))
    #     print(f"diffusive hindrance factor of {ion} = {value(b.hindrance_factor_diffusive_comp[0, ion])}")
    #     print(f"Pore diffusion coefficient = {value(b.diffus_pore_comp[0, ion])}")
    #     print(f"Bulk diffusion coefficient = {value(b.feed_side.properties_in[0].diffus_phase_comp['Liq', ion])}")
    #     print(f"convective hindrance factor of {ion} = {value(b.hindrance_factor_convective_comp[0, ion])}")
    #     print(f"steric partitioning factor of {ion} = {value(b.partition_factor_steric_comp[0, ion])}")
    #     print(f"Born partitioning factor of {ion} = {value(b.partition_factor_born_solvation_comp[0, ion])}")
    #     print(f"Gibbs energy of solvation {ion} = {value(b.gibbs_solvation_comp[0, ion])}")

    #
    # # # check values

