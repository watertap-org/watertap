from pyomo.environ import Param,TransformationFactory,exp
from pyomo.network import Arc
from pyomo.environ import units as pyunits
from idaes.core.util.scaling import constraint_scaling_transform
from idaes.core.util.initialization import propagate_state
from idaes.models.unit_models import Product, Feed
from idaes.core import UnitModelCostingBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from watertap.unit_models.pressure_changer import Pump, EnergyRecoveryDevice
from watertap.costing import WaterTAPCosting
from watertap.unit_models.reverse_osmosis_0D import (ReverseOsmosis0D, ConcentrationPolarizationType, MassTransferCoefficient, PressureChangeType)

# Import concrete model from Pyomo
from pyomo.environ import ConcreteModel, Var, Reals, Objective, Constraint, value, units, NonNegativeReals
# Import flowsheet block from IDAES core
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.scaling import calculate_scaling_factors, set_scaling_factor
# Import NaCl property model
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
from util import *

import numpy as np
import matplotlib.pyplot as plt

def build(membrane_type='RO'):
    # flowsheet set up
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock() # or NaClParameterBlock

    # Control volume flow blocks
    m.fs.feed = m.fs.properties.build_state_block([0])
    m.fs.perm = m.fs.properties.build_state_block([0])

    m.fs.feed_mass_flow = Var(
            initialize=0,
            domain=NonNegativeReals,
            units=pyunits.kg * pyunits.s ** -1,
            doc="Mass flow rate of feed stream",
            )
    
    m.fs.perm_mass_flow = Var(
            initialize=0,
            domain=NonNegativeReals,
            units=pyunits.kg * pyunits.s ** -1,
            doc="Mass flow rate of permeate stream",
            )
    
    m.fs.feed_vol_flow = Var(
            initialize=0,
            domain=NonNegativeReals,
            units=pyunits.m ** 3 * pyunits.s ** -1,
            doc="Mass flow rate of feed stream",
            )
    
    m.fs.perm_vol_flow = Var(
            initialize=0,
            domain=NonNegativeReals,
            units=pyunits.m ** 3 * pyunits.s ** -1,
            doc="Mass flow rate of permeate stream",
            )

    m.fs.flux_mass_phase_comp = Var(
                initialize= 0,
                units=pyunits.kg
                * pyunits.m ** -2
                * pyunits.second ** -1,
                doc="Mass flux across membrane at inlet and outlet",
            )
    
    m.fs.salt_flux_mass_phase_comp = Var(
                initialize= 0,
                units=pyunits.kg
                * pyunits.m ** -2
                * pyunits.second ** -1,
                doc="Mass flux across membrane at inlet and outlet",
            )
    
    m.fs.leaky_flux_mass_phase_comp = Var(
                initialize= 0,
                units=pyunits.kg
                * pyunits.m ** -2
                * pyunits.second ** -1,
                doc="Mass flux across membrane at inlet and outlet",
            )
    
    m.fs.leaky_salt_flux_mass_phase_comp = Var(
                initialize= 0,
                units=pyunits.kg
                * pyunits.m ** -2
                * pyunits.second ** -1,
                doc="Mass flux across membrane at inlet and outlet",
            )

    m.fs.A = Var(
            initialize=1e-8,
            domain=NonNegativeReals,
            units=pyunits.m
                    * pyunits.second ** -1
                    * pyunits.kPa ** -1,
            doc="Water permeability coefficient of the membrane",
        )
    
    m.fs.B = Var(
            initialize=1e-6,
            domain=NonNegativeReals,
            units=pyunits.m
                    * pyunits.second ** -1,
            doc="Salt permeability coefficient of the membrane",
        )
    
    m.fs.A_prime = Var(
            initialize=1e-8,
            domain=NonNegativeReals,
            units=pyunits.m
                    * pyunits.second ** -1
                    * pyunits.kPa ** -1,
            doc="Water permeability coefficient of the membrane",
        )
    
    m.fs.B_prime = Var(
            initialize=1e-6,
            domain=NonNegativeReals,
            units=pyunits.m
                    * pyunits.second ** -1,
            doc="Salt permeability coefficient of the membrane",
        )
    
    m.fs.reflect_coeff = Var(
            initialize=1,
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Reflection coefficient of the membrane",
    )
    
    m.fs.alpha = Var(
            initialize=1,
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Alpha coefficient of the membrane",
    )
    
    m.fs.partition_coeff = Var(
            initialize=0.1,
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc='Partition coefficient of the membrane',
    )
    
    m.fs.rejection = Var(
            initialize=0.9,
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc='Solute rejection of the membrane',
    )
    
    m.fs.recovery_vol_phase = Var(
            initialize=0.5,
            domain=NonNegativeReals,
            bounds=(1e-2, 1 - 1e-6),
            units=pyunits.dimensionless,
            doc="Volumetric recovery rate",
    )
    
    m.fs.obs_salt_perm = Var(
            initialize=0.1,
            domain=NonNegativeReals,
            units=pyunits.kg
                    * pyunits.second ** -1 * pyunits.kPa ** -1,
            doc="Reflection coefficient of the membrane",
    )

    if membrane_type == 'Leaky':

        m.fs.alpha_constraint = Constraint(expr=m.fs.alpha == (1 - m.fs.reflect_coeff) / m.fs.obs_salt_perm)

        m.fs.leaky_flux_constraint = Constraint(expr=m.fs.leaky_flux_mass_phase_comp == m.fs.A * (
            (m.fs.feed[0].pressure - m.fs.perm[0].pressure) - m.fs.reflect_coeff*(
            m.fs.feed[0].pressure_osm_phase['Liq'] - m.fs.perm[0].pressure_osm_phase['Liq']))
            )
    
        # m.fs.leaky_salt_flux_constraint = Constraint(expr=m.fs.leaky_salt_flux_mass_phase_comp == m.fs.obs_salt_perm * 
        #         (m.fs.feed[0].flow_mass_phase_comp['Liq', 'TDS'] - m.fs.perm[0].flow_mass_phase_comp['Liq', 'TDS']) + 
        #         (1-m.fs.reflect_coeff)* (m.fs.leaky_flux_mass_phase_comp*m.fs.feed[0].flow_mass_phase_comp['Liq', 'TDS'])
        #         )
        
        m.fs.leaky_salt_flux_constraint = Constraint(expr=m.fs.leaky_salt_flux_mass_phase_comp == m.fs.obs_salt_perm * 
                (m.fs.feed[0].pressure_osm_phase['Liq'] - m.fs.perm[0].pressure_osm_phase['Liq']) + 
                (1-m.fs.reflect_coeff)* (m.fs.leaky_flux_mass_phase_comp*m.fs.feed[0].flow_mass_phase_comp['Liq', 'TDS'])
                )
        
        m.fs.perm_salt_flow_mass_phase_comp_constraint = Constraint(expr=m.fs.perm[0].flow_mass_phase_comp['Liq', 'TDS'] == m.fs.leaky_salt_flux_mass_phase_comp)
        m.fs.perm_flow_mass_phase_comp_constraint = Constraint(expr=m.fs.perm[0].flow_mass_phase_comp['Liq', 'H2O'] == m.fs.leaky_flux_mass_phase_comp)
    else:
        m.fs.flux_constraint = Constraint(expr=m.fs.flux_mass_phase_comp == m.fs.A * (
                (m.fs.feed[0].pressure - m.fs.perm[0].pressure) - (
                m.fs.feed[0].pressure_osm_phase['Liq'] - m.fs.perm[0].pressure_osm_phase['Liq']))
                )
        
        m.fs.salt_flux_constraint = Constraint(expr=m.fs.salt_flux_mass_phase_comp == m.fs.B * (
                (m.fs.feed[0].flow_mass_phase_comp['Liq', 'TDS'] - m.fs.perm[0].flow_mass_phase_comp['Liq', 'TDS']))
                )
    
        m.fs.perm_salt_flow_mass_phase_comp_constraint = Constraint(expr=m.fs.perm[0].flow_mass_phase_comp['Liq', 'TDS'] == m.fs.salt_flux_mass_phase_comp)
        m.fs.perm_flow_mass_phase_comp_constraint = Constraint(expr=m.fs.perm[0].flow_mass_phase_comp['Liq', 'H2O'] == m.fs.flux_mass_phase_comp)

    m.fs.feed_mass_flow_constraint = Constraint(expr=m.fs.feed_mass_flow == (m.fs.feed[0].flow_mass_phase_comp['Liq', 'H2O'] + m.fs.feed[0].flow_mass_phase_comp['Liq', 'TDS']))
    m.fs.perm_mass_flow_constraint = Constraint(expr=m.fs.perm_mass_flow == (m.fs.perm[0].flow_mass_phase_comp['Liq', 'H2O'] + m.fs.perm[0].flow_mass_phase_comp['Liq', 'TDS']))
    m.fs.feed_vol_flow_constraint = Constraint(expr=m.fs.feed_vol_flow == (m.fs.feed_mass_flow) / m.fs.feed[0].dens_mass_solvent)
    m.fs.perm_vol_flow_constraint = Constraint(expr=m.fs.perm_vol_flow == (m.fs.perm_mass_flow) / m.fs.perm[0].dens_mass_solvent)

    m.fs.rejection_constraint = Constraint(expr= m.fs.rejection == 1 - (m.fs.perm[0].flow_mass_phase_comp['Liq', 'TDS'] / m.fs.feed[0].flow_mass_phase_comp['Liq', 'TDS']))
    
    return m

def set_operating_conditions(m, tds=0.035):
    # fix state variables
    m.fs.feed[0].temperature.fix(273 + 25)                      # temperature (K)
    m.fs.feed[0].pressure.fix(75e5)                           # pressure (Pa)
    m.fs.feed[0].flow_mass_phase_comp['Liq', 'H2O'].fix(1-tds)  # mass flowrate of H2O (kg/s)
    m.fs.feed[0].flow_mass_phase_comp['Liq', 'TDS'].fix(tds)  # mass flowrate of TDS (kg/s)
    m.fs.feed[0].pressure_osm_phase

    m.fs.perm[0].temperature.fix(273 + 25)                      # temperature (K)
    m.fs.perm[0].pressure.fix(101325)                           # pressure (Pa)
    m.fs.perm[0].pressure_osm_phase

    m.fs.A.fix(1e-9)
    m.fs.A_prime.fix(1e-9)
    m.fs.B.fix(1e-2)
    # m.fs.obs_salt_perm.fix(1e-3)
    m.fs.obs_salt_perm.fix(5e-11)
#     m.fs.reflect_coeff.fix(0.1)
    m.fs.alpha.fix(2e10)

    return m

def set_scaling(m):
    # set default property values
    m.fs.properties.set_default_scaling("flow_mass_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "TDS")
    )

def param_sweep():
    solver = get_solver()
    m = build()
    fig, ax = plt.subplots(figsize=(7, 6))
    tds_range = np.linspace(0.01, 0.05, 2)
    molar_range = [0.1, 0.25, 0.5]
    flux_results = []
    leaky_flux_results = []
    for idx, mol in enumerate(molar_range):
        set_operating_conditions(m, mol*54.88/1000)
        auto_scale(m)
        results = solver.solve(m)
        print_results(m)
        flux_results.append(value(units.convert(m.fs.flux_mass_phase_comp / m.fs.perm[0].dens_mass_solvent, to_units=pyunits.liter * pyunits.m ** -2 * pyunits.hour ** -1)))
        leaky_flux_results.append(value(units.convert(m.fs.leaky_flux_mass_phase_comp / m.fs.perm[0].dens_mass_solvent, to_units=pyunits.liter * pyunits.m ** -2 * pyunits.hour ** -1)))

    plot_data(molar_range, flux_results, flux_results, ax, label=f'{mol*54.88*100:.0f} g/L')
    plt.show()

def reflection_sensitivity():
    solver = get_solver()
    m = build(membrane_type='Leaky')
    set_operating_conditions(m, 0.035)
    auto_scale(m)
    fig, ax = plt.subplots(figsize=(8, 6))
    ref_coeff_range = np.linspace(0.5, 1, 10)
    alpha_range = np.linspace(1E9, 1E10, 10)
    flux_results = []
    salt_flux_results = []
    rejection_results = []
    reports = []
    for idx, alpha_coeff in enumerate(alpha_range):
        m.fs.alpha.fix(alpha_coeff)
        results = solver.solve(m)
        print_results(m)
        flux_results.append(value(units.convert(m.fs.leaky_flux_mass_phase_comp / m.fs.perm[0].dens_mass_solvent, to_units=pyunits.liter * pyunits.m ** -2 * pyunits.hour ** -1)))
        salt_flux_results.append(value(units.convert(m.fs.leaky_salt_flux_mass_phase_comp, to_units=pyunits.kg * pyunits.m ** -2 * pyunits.hour ** -1)))
        rejection_results.append(value(m.fs.rejection*100))
        reports.append(create_report(m))

    report = pd.concat(reports)
    report.to_csv('/Users/zbinger/leaky_RO_alpha.csv')
    print(report)
    plot_data(alpha_range, flux_results, salt_flux_results, rejection_results, ax, xlabel=r'Alpha $(\frac{1-\sigma}{\omega^\prime})$', ylabel=r'Water Flux $(\frac{L}{m^{2}  hr})$', ylabel2=r'Salt Flux $(\frac{kg}{m^{2}  hr})$')
    fig.tight_layout()
    fig.savefig('/Users/zbinger/leaky_RO_alpha.png', dpi=900)
    plt.show()

def main():
    solver = get_solver()
    m = build(membrane_type='Leaky')
    set_operating_conditions(m, 0.035)
    auto_scale(m)
    results = solver.solve(m)
    print_results(m)
    print(value(pyunits.convert(m.fs.perm[0].flow_mass_phase_comp['Liq', 'TDS'], to_units=pyunits.g / pyunits.s)))
    print(value(pyunits.convert(m.fs.feed[0].flow_mass_phase_comp['Liq', 'TDS'], to_units=pyunits.g / pyunits.s)))
    print(value(m.fs.rejection))
    df = create_report(m)
    print(df)

if __name__ == "__main__":
#     main()
    # param_sweep()
    reflection_sensitivity()