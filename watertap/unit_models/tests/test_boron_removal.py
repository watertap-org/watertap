#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
Tests for anaerobic digestor example.

Verified against results from:

Rosen, C. and Jeppsson, U., 2006.
Aspects on ADM1 Implementation within the BSM2 Framework.
Department of Industrial Electrical Engineering and Automation, Lund University, Lund, Sweden, pp.1-35.

"""

import pytest
from pyomo.environ import (
    ConcreteModel,
)

from idaes.core import (
    FlowsheetBlock,
)

from idaes.core.solvers import get_solver

from watertap.unit_models.anaerobic_digestor import AD
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

from idaes.core import UnitModelCostingBlock
from watertap.costing import WaterTAPCosting

# Imports from idaes core
from idaes.core import AqueousPhase
from idaes.core.base.components import Solvent, Solute, Cation, Anion
from idaes.core.base.phases import PhaseType as PT

# Imports from idaes generic models
from idaes.models.properties.modular_properties.pure.ConstantProperties import Constant
from idaes.models.properties.modular_properties.state_definitions import FpcTP
from idaes.models.properties.modular_properties.eos.ideal import Ideal

# Importing the enum for concentration unit basis used in the 'get_concentration_term' function
from idaes.models.properties.modular_properties.base.generic_reaction import (
    ConcentrationForm,
)

# Import the object/function for heat of reaction
from idaes.models.properties.modular_properties.reactions.dh_rxn import constant_dh_rxn

# Import safe log power law equation
from idaes.models.properties.modular_properties.reactions.equilibrium_forms import (
    log_power_law_equil,
)

# Import k-value functions
from idaes.models.properties.modular_properties.reactions.equilibrium_constant import (
    van_t_hoff,
)

# Import the idaes objects for Generic Properties and Reactions
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)

from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock

from watertap.unit_models.boron_removal import BoronRemoval
from pyomo.environ import (
    ConcreteModel,
    assert_optimal_termination,
    value,
    Param,
    Var,
    units as pyunits,
    Constraint,
)
from idaes.core import (
    FlowsheetBlock,
)
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
import idaes.core.util.scaling as iscale
from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver
import re

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
class TestBoronRemoval_IonPropPack_Min(UnitTestHarness):
    @pytest.mark.unit
    def configure(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        ion_dict = {
            "solute_list": ["B[OH]3", "B[OH]4_-", "Na_+"],
            "mw_data": {
                "H2O": 18e-3,
                "B[OH]3": 61.83e-3,
                "B[OH]4_-": 78.83e-3,
                "Na_+": 23e-3,
            },
            "charge": {
                "B[OH]4_-": -1,
                "Na_+": 1,
                "B[OH]3": 0,
            },
        }

        m.fs.properties = MCASParameterBlock(**ion_dict)

        map = {
            "boron_name": "B[OH]3",  # [is required]
            "borate_name": "B[OH]4_-",  # [is required]
            "caustic_additive": {
                "cation_name": "Na_+",  # [is optional]
                "mw_additive": (40, pyunits.g / pyunits.mol),  # [is required]
                "moles_cation_per_additive": 1,  # [is required]
            },
        }
        m.fs.unit = BoronRemoval(
            property_package=m.fs.properties, chemical_mapping_data=map
        )

        # Set the operating conditions
        m.fs.unit.inlet.pressure.fix(101325)
        m.fs.unit.inlet.temperature.fix(298.15)

        state = {
            "H2O": 100,
            "H_+": 1e-7,
            "OH_-": 1e-7,
            "B[OH]3": 2e-4,
            "B[OH]4_-": 1e-6,
            "Na_+": 1e-4,
            "HCO3_-": 1e-4,
        }

        for j in state:
            idx = (0, "Liq", j)
            if idx in m.fs.unit.inlet.flow_mol_phase_comp:
                m.fs.unit.inlet.flow_mol_phase_comp[idx].fix(state[j])
        m.fs.unit.caustic_dose_rate.fix(1.8e-5)
        m.fs.unit.reactor_volume.fix(1)

        # Set scaling factors for badly scaled variables
        for j in state:
            idx = (0, "Liq", j)
            if idx in m.fs.unit.inlet.flow_mol_phase_comp:
                m.fs.properties.set_default_scaling(
                    "flow_mol_phase_comp", 1 / state[j], index=("Liq", j)
                )

        iscale.calculate_scaling_factors(m.fs)

        self.unit_model_block = m.fs.unit

        self.unit_solutions[
            m.fs.unit.outlet.flow_mol_phase_comp[0, "Liq", m.fs.unit.boron_name_id]
        ] = 1.5062434e-5
        self.unit_solutions[
            m.fs.unit.outlet.flow_mol_phase_comp[0, "Liq", m.fs.unit.borate_name_id]
        ] = 1.859391e-4
        self.unit_solutions[m.fs.unit.caustic_dose_rate[0]] = 1.8e-5
        self.unit_solutions[m.fs.unit.outlet_pH()] = 7
