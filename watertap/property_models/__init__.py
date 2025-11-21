#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
# IDAES imports
from idaes.core import MaterialBalanceType, EnergyBalanceType, MaterialFlowBasis

from .multicomp_aq_sol_prop_pack import *
from .NaCl_prop_pack import NaClParameterBlock, NaClParameterData
from .NaCl_T_dep_prop_pack import NaClTDepParameterBlock, NaClTDepParameterData
from .seawater_prop_pack import SeawaterParameterBlock, SeawaterParameterData
from .water_prop_pack import WaterParameterBlock, WaterParameterData
from .zero_order_properties import ZOParameterBlock, ZOParameterData, ZOStateBlock

# Unit Specific
from .unit_specific.coagulation_prop_pack import (
    CoagulationParameterBlock,
    CoagulationParameterData,
    CoagulationStateBlock,
)
from .unit_specific.cryst_prop_pack import (
    CrystallizerParameterBlock,
    CrystallizerParameterData,
    HeatOfCrystallizationModel,
)
from .unit_specific.NDMA_prop_pack import NDMAParameterBlock, NDMAParameterData

# Activated Sludge
from .unit_specific.activated_sludge.asm1_properties import (
    ASM1ParameterBlock,
    ASM1StateBlock,
    ASM1PropertiesScaler,
    ASM1ParameterData,
)
from .unit_specific.activated_sludge.asm1_reactions import (
    ASM1ReactionParameterBlock,
    ASM1ReactionScaler,
    ASM1ReactionBlock,
)
from .unit_specific.activated_sludge.asm2d_properties import (
    ASM2dParameterBlock,
    ASM2dStateBlock,
)
from .unit_specific.activated_sludge.asm2d_reactions import (
    ASM2dReactionParameterBlock,
    ASM2dReactionBlock,
)
from .unit_specific.activated_sludge.asm3_properties import (
    ASM3ParameterBlock,
    ASM3PropertiesScaler,
    ASM3StateBlock,
)
from .unit_specific.activated_sludge.asm3_reactions import (
    ASM3ReactionParameterBlock,
    ASM3ReactionScaler,
    ASM3ReactionBlock,
)
from .unit_specific.activated_sludge.modified_asm2d_properties import (
    ModifiedASM2dParameterBlock,
    ModifiedASM2dPropertiesScaler,
    ModifiedASM2dStateBlock,
)
from .unit_specific.activated_sludge.modified_asm2d_reactions import (
    ModifiedASM2dReactionParameterBlock,
    ModifiedASM2dReactionScaler,
    ModifiedASM2dReactionBlock,
)

# Anaerobic Digestion
from .unit_specific.anaerobic_digestion.adm1_properties_vapor import (
    ADM1_vaporParameterBlock,
    ADM1VaporPropertiesScaler,
    ADM1_vaporStateBlock,
)
from .unit_specific.anaerobic_digestion.adm1_properties import (
    ADM1ParameterBlock,
    ADM1PropertiesScaler,
    ADM1StateBlock,
)
from .unit_specific.anaerobic_digestion.adm1_reactions import (
    ADM1ReactionParameterBlock,
    ADM1ReactionScaler,
    ADM1ReactionBlock,
)
from .unit_specific.anaerobic_digestion.modified_adm1_properties import (
    ModifiedADM1ParameterBlock,
    ModifiedADM1PropertiesScaler,
    ModifiedADM1StateBlock,
)
from .unit_specific.anaerobic_digestion.modified_adm1_reactions import (
    ModifiedADM1ReactionParameterBlock,
    ModifiedADM1ReactionScaler,
    ModifiedADM1ReactionBlock,
)
