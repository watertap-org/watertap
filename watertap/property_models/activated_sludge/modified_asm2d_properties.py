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
Thermophysical property package to be used in conjunction with modified ASM2d reactions.

Reference:
X. Flores-Alsina, K. Solon, C.K. Mbamba, S. Tait, K.V. Gernaey, U. Jeppsson, D.J. Batstone,
Modelling phosphorus (P), sulfur (S) and iron (Fe) interactions fordynamic simulations of anaerobic digestion processes,
Water Research. 95 (2016) 370-382. https://www.sciencedirect.com/science/article/pii/S0043135416301397

"""

# Import Pyomo libraries
import pyomo.environ as pyo

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    MaterialFlowBasis,
    PhysicalParameterBlock,
    StateBlockData,
    StateBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    LiquidPhase,
    Solute,
    Solvent,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import fix_state_vars, revert_state_vars
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale

# Some more information about this module
__author__ = "Marcus Holly, Adam Atia, Xinhong Liu"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("ModifiedASM2dParameterBlock")
class ModifiedASM2dParameterData(PhysicalParameterBlock):
    """
    Property Parameter Block Class
    """

    def build(self):
        """
        Callable method for Block construction.
        """
        super().build()

        self._state_block_class = ModifiedASM2dStateBlock

        # Add Phase objects
        self.Liq = LiquidPhase()

        # Add Component objects
        self.H2O = Solvent()

        # Soluble species
        self.S_A = Solute(
            doc="Fermentation products, considered to be acetate. [kg COD/m^3]"
        )
        self.S_F = Solute(
            doc="Fermentable, readily bio-degradable organic substrates. [kg COD/m^3]"
        )
        self.S_I = Solute(doc="Inert soluble organic material. [kg COD/m^3]")
        self.S_N2 = Solute(
            doc="Dinitrogen, N2. SN2 is assumed to be the only nitrogenous product of denitrification. [kg N/m^3]"
        )
        self.S_NH4 = Solute(doc="Ammonium plus ammonia nitrogen. [kg N/m^3]")
        self.S_NO3 = Solute(
            doc="Nitrate plus nitrite nitrogen (N03' + N02' -N). SN03 is assumed to include nitrate as well as nitrite nitrogen. [kg N/m^3]"
        )
        self.S_O2 = Solute(doc="Dissolved oxygen. [kg O2/m^3]")
        self.S_PO4 = Solute(
            doc="Inorganic soluble phosphorus, primarily ortho-phosphates. [kg P/m^3]"
        )
        self.S_K = Solute(doc="Potassium, [kg K/m^3]")
        self.S_Mg = Solute(doc="Magnesium, [kg Mg/m^3]")
        self.S_IC = Solute(doc="Inorganic carbon, [kg C/m^3]")

        # Particulate species
        self.X_AUT = Solute(doc="Autotrophic nitrifying organisms. [kg COD/m^3]")
        self.X_H = Solute(doc="Heterotrophic organisms. [kg COD/m^3]")
        self.X_I = Solute(doc="Inert particulate organic material. [kg COD/m^3]")
        self.X_PAO = Solute(doc="Phosphate-accumulating organisms. [kg COD/m^3]")
        self.X_PHA = Solute(
            doc="A cell internal storage product of phosphorus-accumulating organisms, primarily comprising poly-hydroxy-alkanoates (PHA). [kg COD/m^3]"
        )
        self.X_PP = Solute(doc="Poly-phosphate. [kg P/m^3]")
        self.X_S = Solute(doc="Slowly biodegradable substrates. [kg COD/m^3]")

        # Create sets for use across ASM models and associated unit models
        self.non_particulate_component_set = pyo.Set(
            initialize=[
                "S_A",
                "S_F",
                "S_I",
                "S_N2",
                "S_NH4",
                "S_NO3",
                "S_O2",
                "S_PO4",
                "S_K",
                "S_Mg",
                "S_IC",
                "H2O",
            ]
        )
        self.particulate_component_set = pyo.Set(
            initialize=[
                "X_AUT",
                "X_H",
                "X_I",
                "X_PAO",
                "X_PHA",
                "X_PP",
                "X_S",
            ]
        )
        self.tss_component_set = pyo.Set(
            initialize=[
                "X_AUT",
                "X_H",
                "X_I",
                "X_PAO",
                "X_PHA",
                "X_PP",
                "X_S",
            ]
        )

        # Heat capacity of water
        self.cp_mass = pyo.Param(
            mutable=False,
            initialize=4182,
            doc="Specific heat capacity of water",
            units=pyo.units.J / pyo.units.kg / pyo.units.K,
        )
        # Density of water
        self.dens_mass = pyo.Param(
            mutable=False,
            initialize=997,
            doc="Density of water",
            units=pyo.units.kg / pyo.units.m**3,
        )

        # Thermodynamic reference state
        self.pressure_ref = pyo.Param(
            within=pyo.PositiveReals,
            mutable=True,
            default=101325.0,
            doc="Reference pressure",
            units=pyo.units.Pa,
        )
        self.temperature_ref = pyo.Param(
            within=pyo.PositiveReals,
            mutable=True,
            default=298.15,
            doc="Reference temperature",
            units=pyo.units.K,
        )

        # COD to VSS coefficients
        self.CODtoVSS_XI = pyo.Var(
            initialize=1.5686,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="mass COD per mass VSS of XI",
        )
        self.CODtoVSS_XS = pyo.Var(
            initialize=1.5686,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="mass COD per mass VSS of XS",
        )
        self.CODtoVSS_XBM = pyo.Var(
            initialize=1.3072,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="mass COD per mass VSS of biomass",
        )
        self.CODtoVSS_XPHA = pyo.Var(
            initialize=1.9608,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="mass COD per mass VSS of XPHA",
        )
        # Inorganic solids parameters
        self.ISS_P = pyo.Var(
            initialize=3.23,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="mass ISS per mass P",
        )
        self.f_ISS_BM = pyo.Var(
            initialize=0.15,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="ISS fractional content of biomass",
        )

        # Fix Vars that are treated as Params
        for v in self.component_objects(pyo.Var):
            v.fix()

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "flow_vol": {"method": None},
                "pressure": {"method": None},
                "temperature": {"method": None},
                "conc_mass_comp": {"method": None},
            }
        )
        obj.define_custom_properties(
            {
                "VSS": {"method": "_VSS"},
                "ISS": {"method": "_ISS"},
                "TSS": {"method": "_TSS"},
            }
        )
        obj.add_default_units(
            {
                "time": pyo.units.s,
                "length": pyo.units.m,
                "mass": pyo.units.kg,
                "amount": pyo.units.kmol,
                "temperature": pyo.units.K,
            }
        )


class _ModifiedASM2dStateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """

    def initialize(
        self,
        state_args=None,
        state_vars_fixed=False,
        hold_state=False,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        Initialization routine for property package.

        Keyword Arguments:
            state_args : Dictionary with initial guesses for the state vars
                         chosen. Note that if this method is triggered
                         through the control volume, and if initial guesses
                         were not provided at the unit model level, the
                         control volume passes the inlet values as initial
                         guess.The keys for the state_args dictionary are:
            flow_vol : value at which to initialize total volumetric flow (default=None)
            conc_mass_comp : value at which to initialize component concentrations (default=None)
            pressure : value at which to initialize pressure (default=None)
            temperature : value at which to initialize temperature (default=None)
            outlvl : sets output level of initialization routine
            state_vars_fixed: Flag to denote if state vars have already been fixed.
                              True - states have already been fixed and
                              initialization does not need to worry
                              about fixing and unfixing variables.
                              False - states have not been fixed. The state
                              block will deal with fixing/unfixing.
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)
            hold_state : flag indicating whether the initialization routine
                         should unfix any state variables fixed during
                         initialization (default=False).
                         True - states variables are not unfixed, and
                         a dict of returned containing flags for
                         which states were fixed during initialization.
                         False - state variables are unfixed after
                         initialization by calling the release_state method.

        Returns:
            If hold_states is True, returns a dict containing flags for
            which states were fixed during initialization.
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")

        if state_vars_fixed is False:
            # Fix state variables if not already fixed
            flags = fix_state_vars(self, state_args)

        else:
            # Check when the state vars are fixed already result in dof 0
            for k in self.keys():
                if degrees_of_freedom(self[k]) != 0:
                    raise Exception(
                        "State vars fixed but degrees of freedom "
                        "for state block is not zero during "
                        "initialization."
                    )

        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                self.release_state(flags)

        init_log.info("Initialization Complete.")

    def release_state(self, flags, outlvl=idaeslog.NOTSET):
        """
        Method to relase state variables fixed during initialization.

        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state=True.
            outlvl : sets output level of of logging
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")

        if flags is None:
            return
        # Unfix state variables
        revert_state_vars(self, flags)
        init_log.info("State Released.")


@declare_process_block_class(
    "ModifiedASM2dStateBlock", block_class=_ModifiedASM2dStateBlock
)
class ModifiedASM2dStateBlockData(StateBlockData):
    """
    StateBlock for calculating thermophysical proeprties associated with the ASM2d
    reaction system.
    """

    def build(self):
        """
        Callable method for Block construction
        """
        super().build()

        # Create state variables
        self.flow_vol = pyo.Var(
            initialize=1.0,
            domain=pyo.NonNegativeReals,
            doc="Total volumetric flowrate",
            units=pyo.units.m**3 / pyo.units.s,
        )
        self.pressure = pyo.Var(
            domain=pyo.NonNegativeReals,
            initialize=101325.0,
            bounds=(1e3, 1e6),
            doc="Pressure",
            units=pyo.units.Pa,
        )
        self.temperature = pyo.Var(
            domain=pyo.NonNegativeReals,
            initialize=298.15,
            bounds=(273.15, 323.15),
            doc="Temperature",
            units=pyo.units.K,
        )
        self.conc_mass_comp = pyo.Var(
            self.params.solute_set,
            domain=pyo.NonNegativeReals,
            initialize=0.1,
            doc="Component mass concentrations",
            units=pyo.units.kg / pyo.units.m**3,
        )

    # On-demand properties
    def _VSS(self):
        self.VSS = pyo.Var(
            initialize=1,
            domain=pyo.NonNegativeReals,
            doc="Volatile suspended solids",
            units=pyo.units.kg / pyo.units.m**3,
        )

        # TODO: X_SRB not included yet in biomass term summation
        def rule_VSS(b):
            return (
                b.VSS
                == b.conc_mass_comp["X_I"] / b.params.CODtoVSS_XI
                + b.conc_mass_comp["X_S"] / b.params.CODtoVSS_XS
                + (
                    b.conc_mass_comp["X_H"]
                    + b.conc_mass_comp["X_PAO"]
                    + b.conc_mass_comp["X_AUT"]
                )
                / b.params.CODtoVSS_XBM
                + b.conc_mass_comp["X_PHA"] / b.params.CODtoVSS_XPHA
            )

        self.eq_VSS = pyo.Constraint(rule=rule_VSS)

    def _ISS(self):
        self.ISS = pyo.Var(
            initialize=1,
            domain=pyo.NonNegativeReals,
            doc="Inorganic suspended solids",
            units=pyo.units.kg / pyo.units.m**3,
        )

        # TODO: Several HFO and other terms omitted since not included yet.
        def rule_ISS(b):
            return (
                b.ISS
                == b.params.f_ISS_BM
                * (
                    b.conc_mass_comp["X_H"]
                    + b.conc_mass_comp["X_PAO"]
                    + b.conc_mass_comp["X_AUT"]
                )
                / b.params.CODtoVSS_XBM
                + b.params.ISS_P * b.conc_mass_comp["X_PP"]
            )

        self.eq_ISS = pyo.Constraint(rule=rule_ISS)

    def _TSS(self):
        self.TSS = pyo.Var(
            initialize=1,
            domain=pyo.NonNegativeReals,
            doc="Total suspended solids",
            units=pyo.units.kg / pyo.units.m**3,
        )

        def rule_TSS(b):
            return b.TSS == b.VSS + b.ISS

        self.eq_TSS = pyo.Constraint(rule=rule_TSS)

    def get_material_flow_terms(self, p, j):
        if j == "H2O":
            return self.flow_vol * self.params.dens_mass
        else:
            return self.flow_vol * self.conc_mass_comp[j]

    def get_enthalpy_flow_terms(self, p):
        return (
            self.flow_vol
            * self.params.dens_mass
            * self.params.cp_mass
            * (self.temperature - self.params.temperature_ref)
        )

    def get_material_density_terms(self, p, j):
        if j == "H2O":
            return self.params.dens_mass
        else:
            return self.conc_mass_comp[j]

    def get_energy_density_terms(self, p):
        return (
            self.params.dens_mass
            * self.params.cp_mass
            * (self.temperature - self.params.temperature_ref)
        )

    def default_material_balance_type(self):
        return MaterialBalanceType.componentPhase

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def define_state_vars(self):
        return {
            "flow_vol": self.flow_vol,
            "conc_mass_comp": self.conc_mass_comp,
            "temperature": self.temperature,
            "pressure": self.pressure,
        }

    def define_display_vars(self):
        return {
            "Volumetric Flowrate": self.flow_vol,
            "Mass Concentration": self.conc_mass_comp,
            "Temperature": self.temperature,
            "Pressure": self.pressure,
        }

    def get_material_flow_basis(self):
        return MaterialFlowBasis.mass

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # TODO: revisit scaling of these new on-demand props
        if self.is_property_constructed("VSS"):
            if iscale.get_scaling_factor(self.VSS) is None:
                iscale.set_scaling_factor(self.VSS, 1)

        if self.is_property_constructed("ISS"):
            if iscale.get_scaling_factor(self.ISS) is None:
                iscale.set_scaling_factor(self.ISS, 1)

        if self.is_property_constructed("TSS"):
            if iscale.get_scaling_factor(self.TSS) is None:
                iscale.set_scaling_factor(self.TSS, 1)
