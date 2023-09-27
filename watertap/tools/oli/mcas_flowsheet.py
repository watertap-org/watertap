###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

__author__ = "Paul Vecchiarelli, Adam Atia"

from pyomo.environ import ConcreteModel, assert_optimal_termination
from idaes.core import FlowsheetBlock
from idaes.core.util.scaling import calculate_scaling_factors
from idaes.core.solvers import get_solver
import watertap.property_models.multicomp_aq_sol_prop_pack as props

from watertap.tools.oli.watertap_to_oli_helper_functions import watertap_to_oli

solver = get_solver()

def mcas_flowsheet(source_water: dict):
    """
    Create a simple MCAS Flowsheet from a feedwater dictionary.
    
    :param source_water: Dictionary containing temperature, pressure, and solutes:concentrations.
    """
    
    for key in ["components", "temperature", "pressure"]:
        if key not in source_water.keys():
            raise RuntimeError(f" Key '{key}' missing from source_water. Initialization aborted.")
            
    components = {watertap_to_oli(comp): source_water["components"][comp] for comp in source_water["components"]}    
    
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    mcas_input = create_mcas_input(components)
    m.fs.properties = props.MCASParameterBlock(**mcas_input)
    stream = m.fs.stream = m.fs.properties.build_state_block([0])
    stream[0].temperature.fix(source_water["temperature"])
    stream[0].pressure.fix(source_water["pressure"])
    
    # TODO: equations to control volume and density in calculate_mass_frac_phase_comp
    mass_fracs = calculate_mass_frac_phase_comp(components)
    mass_fracs.update({("flow_vol_phase", "Liq"): 1e-3})
    stream.calculate_state(var_args=mass_fracs, hold_state=True)
    
    stream[0].conc_mass_phase_comp

    # scaling
    setup_scaling(m, stream[0].mass_frac_phase_comp)

    stream.initialize()
    solver = get_solver()
    result = solver.solve(m)
    assert_optimal_termination(result)
            
    return m

def create_mcas_input(components):
    """
    Builds MCAS property package inputs automatically if WaterTAP (or known OLI) names are passed.
    
    :param components: dictionary containing solute concentrations in mg/L.
    """
    
    mcas_input = {"solute_list": [comp.oli_name for comp in components],
                  "mw_data": {comp.oli_name: comp.molar_mass * 1e-3 for comp in components},
                  "charge": {comp.oli_name: comp.charge for comp in components if comp.charge != 0},
                  }
    
    return mcas_input

def calculate_mass_frac_phase_comp(components:dict, volume:float=1, density:float=1e6):
    """
    Calculates mass fraction of a solution from mol/L.
    
    :param components: dictionary containing solute concentrations in mg/L.
    :param volume: volume of solution (in L)
    :param density: density of solution (in mg/L - 1e6 for water)
    """
    
    total_mass = sum(components.values()) + (volume*density)
    
    mass_frac = lambda conc: conc / total_mass
    mass_frac_index = lambda comp: ("mass_frac_phase_comp", ("Liq", comp.oli_name))
    
    return {mass_frac_index(comp): mass_frac(conc) for comp, conc in components.items()}
    
def setup_scaling(m, phase_comp_blk):
    """
    Initializes scaling factors for flow_mol_phase_comp parameter block.
    
    :param m: ConcreteModel object.
    :param phase_comp_blk: object containing all solutes in state block.
    """
    
    comps = set(m.fs.properties.solvent_set | m.fs.properties.solute_set)
    for comp in comps:            
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1, index=("Liq", comp)
        )
        
    calculate_scaling_factors(m)
    return