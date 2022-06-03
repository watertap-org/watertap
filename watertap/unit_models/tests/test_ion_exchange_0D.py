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

import pytest
from pyomo.environ import (
    ConcreteModel,
    Constraint,
    TerminationCondition,
    SolverStatus,
    value,
    Var,
    units as pyunits,
)
from pyomo.network import Port
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from watertap.property_models.IX_prop_pack import IXParameterBlock
from watertap.unit_models.ion_exchange_0D import IonExchange0D

from idaes.core.util import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    unscaled_constraints_generator,
    badly_scaled_var_generator,
)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


#### NOT AN ACTUAL TEST FILE ####
#### JUST MY WORKFLOW ####

m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
ix_in = {
    "solute_list": ["Na_+"],
    "diffusivity_data": {("Liq", "Na_+"): 1.33e-9},
    "mw_data": {"H2O": 18e-3, "Na_+": 23e-3},
    "charge": {"Na_+": 1},
}
m.fs.properties = IXParameterBlock(default=ix_in)
ix = m.fs.unit = IonExchange0D(default={"property_package": m.fs.properties})
m.fs.unit.inlet.pressure[0].fix(101325)
m.fs.unit.inlet.temperature[0].fix(298.15)

mass_flow_in = 10 * pyunits.kg / pyunits.s
feed_mass_frac = {
    "Na_+": 1e-4,
}
for ion, x in feed_mass_frac.items():
    mol_comp_flow = (
        x
        * pyunits.kg
        / pyunits.kg
        * mass_flow_in
        / m.fs.unit.config.property_package.mw_comp[ion]
    )
    m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", ion].fix(mol_comp_flow)

H2O_mass_frac = 1 - sum(x for x in feed_mass_frac.values())
H2O_mol_comp_flow = (
    H2O_mass_frac
    * pyunits.kg
    / pyunits.kg
    * mass_flow_in
    / m.fs.unit.config.property_package.mw_comp["H2O"]
)

ix.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(H2O_mol_comp_flow)
ix.sphericity.fix()
ix.resin_diam.fix()
ix.K_eq["Na_+"].fix(2.5)
ix.resin_max_capacity.fix(5)
ix.bed_porosity.fix()
ix.resin_bulk_dens.fix()
ix.dimensionless_time.fix(1)
ix.lh.fix(0)
ix.sfr.fix(5)
ix.Re.fix(5.3)

calculate_scaling_factors(m)
m.fs.unit.initialize()
print(f"\nDOF = {degrees_of_freedom(ix)}\n")

solver = get_solver()
results = solver.solve(m)

print(f"\nMODEL SOLVE = {results.solver.termination_condition}")
