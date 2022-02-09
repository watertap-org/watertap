from pyomo.environ import ConcreteModel, assert_optimal_termination, value, units as pyunits, check_optimal_termination
from idaes.core import FlowsheetBlock, MomentumBalanceType
import idaes.core.util.scaling as iscale
from idaes.core.util.model_statistics import fixed_variables_set, unfixed_variables_set, unused_variables_set
from pyomo.util.check_units import assert_units_consistent
from watertap.property_models.ion_DSPMDE_prop_pack import DSPMDEParameterBlock, ActivityCoefficientModel, DensityCalculation
from watertap.unit_models.nanofiltration_DSPMDE_0D import NanofiltrationDSPMDE0D, MassTransferCoefficient
from watertap.core.util.initialization import check_dof, check_solve
from watertap.core.util.infeasible import *
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
        "density_calculation": DensityCalculation.constant
    })

    m.fs.stream = m.fs.properties.build_state_block([0], default={'defined_state': True})

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

    # Fix mole flow rates of each ion and water
    for ion, x in feed_mass_frac.items():
        mol_comp_flow = x * pyunits.kg / pyunits.kg * mass_flow_in / m.fs.stream[0].mw_comp[ion]
        m.fs.stream[0].flow_mol_phase_comp['Liq', ion].fix(mol_comp_flow)
    H2O_mass_frac = 1 - sum(x for x in feed_mass_frac.values())
    H2O_mol_comp_flow = H2O_mass_frac * pyunits.kg / pyunits.kg * mass_flow_in / \
                        m.fs.stream[0].mw_comp['H2O']
    m.fs.stream[0].flow_mol_phase_comp['Liq', 'H2O'].fix(H2O_mol_comp_flow)
    m.fs.stream[0].temperature.fix(298.15)
    m.fs.stream[0].pressure.fix(4e5)

    m.fs.stream[0].conc_mol_phase_comp

    iscale.calculate_scaling_factors(m)

    m.fs.stream.initialize()

    res=solver.solve(m,tee=True)
    assert_optimal_termination(res)