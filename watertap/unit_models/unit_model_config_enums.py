from enum import Enum, auto
from idaes.core.util.misc import StrEnum


##############################
# CSTR Injection
##############################
class ElectricityConsumption(Enum):
    """
    none: no electricity consumption
    fixed: assume electricity intensity
    calculated: calculate based on aeration energy equation from BSM2 documentation
    """

    none = auto()
    fixed = auto()
    calculated = auto()

##############################
# Dewatering
##############################
class ActivatedSludgeModelType(Enum):
    """
    ASM1: ASM1 model
    ASM2D: ASM2D model
    modified_ASM2D: modified ASM2D model for ADM1 compatibility
    """

    ASM1 = auto()
    ASM2D = auto()
    modified_ASM2D = auto()

##############################
# Electrodialysis 0D
# Electrodialysis 1D
# Electrodialysis Bipolar 1D
##############################
class LimitingCurrentDensityMethod(Enum):
    InitialValue = 0
    Empirical = 1
    Theoretical = 2


class ElectricalOperationMode(Enum):
    Constant_Current = 0
    Constant_Voltage = 1


class PressureDropMethod(Enum):
    none = 0
    experimental = 1
    Darcy_Weisbach = 2


class FrictionFactorMethod(Enum):
    fixed = 0
    Gurreri = 1
    Kuroda = 2


class HydraulicDiameterMethod(Enum):
    fixed = 0
    spacer_specific_area_known = 1
    conventional = 2

##############################
# Electrodialysis Bipolar 1D
##############################

class LimitingCurrentDensitybpmMethod(Enum):
    InitialValue = 0
    Empirical = 1

##############################
# GAC
##############################

class FilmTransferCoefficientType(Enum):
    # liquid phase film transfer coefficient is a user specified value
    fixed = auto()
    # calculate liquid phase film transfer coefficient from the Gnielinski correlation
    calculated = auto()


class SurfaceDiffusionCoefficientType(Enum):
    fixed = auto()  # surface diffusion coefficient is a user specified value
    calculated = auto()  # calculate surface diffusion coefficient

##############################
# Ion Exchange
##############################

class IonExchangeType(StrEnum):
    anion = "anion"
    cation = "cation"


class RegenerantChem(StrEnum):
    HCl = "HCl"
    NaOH = "NaOH"
    H2SO4 = "H2SO4"
    NaCl = "NaCl"
    MeOH = "MeOH"
    single_use = "single_use"


class IsothermType(StrEnum):
    langmuir = "langmuir"
    freundlich = "freundlich"

##############################
# Nanofiltration DSPMDE 0D
##############################

class MassTransferCoefficient(Enum):
    none = auto()
    # mass transfer coefficient is a user specified value
    fixed = auto()
    # mass transfer coefficient is calculated using spiral wound correlation
    spiral_wound = auto()


class ConcentrationPolarizationType(Enum):
    none = auto()
    calculated = auto()


class ConcentrationPolarizationType(Enum):
    """
    none: no concentration polarization
    fixed: concentration polarization modulus is a user specified value
    calculated: calculate concentration polarization (concentration at membrane interface)
    """

    none = auto()
    fixed = auto()
    calculated = auto()

##############################
# OARO 0D
# OARO 1D
# RO 0D
# RO 1D
##############################


class MassTransferCoefficient(Enum):
    """
    none: mass transfer coefficient not utilized for concentration polarization effect
    fixed: mass transfer coefficient is a user specified value
    calculated: mass transfer coefficient is calculated
    """

    none = auto()
    fixed = auto()
    calculated = auto()
    # TODO: add option for users to define their own relationship?


class TransportModel(Enum):
    """
    SD: Solvent and solute mass flux is calculated using the Solution-Diffusion model
    SKK: Solvent and solute mass flux is calculated using the Spiegler-Kedem-Katchalsky model
    """

    SD = auto()
    SKK = auto()


class ModuleType(Enum):
    """
    flat_sheet: flat-sheet module configuration
    spiral_wound: spiral-wound module configuration
    """

    flat_sheet = auto()
    spiral_wound = auto()


class PressureChangeType(Enum):
    """
    fixed_per_stage: pressure drop across membrane channel is a user-specified value
    fixed_per_unit_length: pressure drop per unit length across membrane channel is a user-specified value
    calculated: pressure drop across membrane channel is calculated
    """

    fixed_per_stage = auto()
    fixed_per_unit_length = auto()
    calculated = auto()


class FrictionFactor(Enum):
    """
    default_by_module_type: Will revert FrictionFactor to either the Darcy's friction factor correlation
    by Guillen & Hoek (flat-plate) or by Schock & Miquel (spiral-wound)
    """

    default_by_module_type = auto()

##############################
# Pump
##############################

class VariableEfficiency(Enum):
    none = auto()  # default is constant efficiency
    flow = auto()  # flow-only correlation
    flow_head = auto()  # flow and head correlation

##############################
# Pressure Exchanger
##############################

class PressureExchangeType(Enum):
    efficiency = auto()
    high_pressure_difference = auto()

##############################
# Steam Heater 0D
##############################

class Mode(Enum):
    HEATER = auto()
    CONDENSER = auto()