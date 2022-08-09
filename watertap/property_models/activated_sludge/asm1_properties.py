#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Thermophysical property package to be used in conjunction with ASM1 reactions.
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
    Component,
    Solute,
    Solvent,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import fix_state_vars, revert_state_vars
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale

# Some more information about this module
__author__ = "Andrew Lee"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("ASM1ParameterBlock")
class ASM1ParameterData(PhysicalParameterBlock):
    """
    Property Parameter Block Class
    """

    def build(self):
        """
        Callable method for Block construction.
        """
        super().build()

        self._state_block_class = ASM1StateBlock

        # Add Phase objects
        self.Liq = LiquidPhase()

        # Add Component objects
        self.H2O = Solvent()

        self.S_I = Solute(doc="Soluble inert organic matter, S_I")
        self.S_S = Solute(doc="Readily biodegradable substrate, S_S")
        self.X_I = Solute(doc="Particulate inert organic matter, X_I")
        self.X_S = Solute(doc="Slowly biodegradable substrate, X_S")
        self.X_BH = Solute(doc="Active heterotrophic biomass, X_B,H")
        self.X_BA = Solute(doc="Active autotrophic biomass, X_B,A")
        self.X_P = Solute(doc="Particulate products arising from biomass decay, X_P")
        self.S_O = Solute(doc="Oxygen, S_O")
        self.S_NO = Solute(doc="Nitrate and nitrite nitrogen, S_NO")
        self.S_NH = Solute(doc="NH4+ + NH3 nitrogen, S_NH")
        self.S_ND = Solute(doc="Soluble biodegradable organic nitrogen, S_ND")
        self.X_ND = Solute(doc="Particulate biodegradable organic nitrogen, X_ND")

        self.S_ALK = Component(doc="Alkalinity, S_ALK")

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

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "flow_vol": {"method": None},
                "pressure": {"method": None},
                "temperature": {"method": None},
                "conc_mass_comp": {"method": None},
                "alkalinity": {"method": None},
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


class _ASM1StateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """

    def initialize(
        blk,
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
                     were not provied at the unit model level, the
                     control volume passes the inlet values as initial
                     guess.The keys for the state_args dictionary are:

                     flow_mol_comp : value at which to initialize component
                                     flows (default=None)
                     pressure : value at which to initialize pressure
                                (default=None)
                     temperature : value at which to initialize temperature
                                  (default=None)
            outlvl : sets output level of initialization routine
            state_vars_fixed: Flag to denote if state vars have already been
                              fixed.
                              - True - states have already been fixed and
                                       initialization does not need to worry
                                       about fixing and unfixing variables.
                             - False - states have not been fixed. The state
                                       block will deal with fixing/unfixing.
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)
            hold_state : flag indicating whether the initialization routine
                         should unfix any state variables fixed during
                         initialization (default=False).
                         - True - states varaibles are not unfixed, and
                                 a dict of returned containing flags for
                                 which states were fixed during
                                 initialization.
                        - False - state variables are unfixed after
                                 initialization by calling the
                                 relase_state method

        Returns:
            If hold_states is True, returns a dict containing flags for
            which states were fixed during initialization.
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="properties")

        if state_vars_fixed is False:
            # Fix state variables if not already fixed
            flags = fix_state_vars(blk, state_args)

        else:
            # Check when the state vars are fixed already result in dof 0
            for k in blk.keys():
                if degrees_of_freedom(blk[k]) != 0:
                    raise Exception(
                        "State vars fixed but degrees of freedom "
                        "for state block is not zero during "
                        "initialization."
                    )

        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                blk.release_state(flags)

        init_log.info("Initialization Complete.")

    def release_state(blk, flags, outlvl=idaeslog.NOTSET):
        """
        Method to relase state variables fixed during initialization.

        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state=True.
            outlvl : sets output level of of logging
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="properties")

        if flags is None:
            return
        # Unfix state variables
        revert_state_vars(blk, flags)
        init_log.info("State Released.")


@declare_process_block_class("ASM1StateBlock", block_class=_ASM1StateBlock)
class ASM1StateBlockData(StateBlockData):
    """
    StateBlock for calculating thermophysical proeprties associated with the ASM1
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
            doc="Total volumentric flowrate",
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
            bounds=(293.15, 323.15),
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
        self.alkalinity = pyo.Var(
            domain=pyo.NonNegativeReals,
            initialize=1,
            doc="Alkalinity in molar concentration",
            units=pyo.units.kmol / pyo.units.m**3,
        )

        # Material and energy flow and density expressions
        def material_flow_expression(self, j):
            if j == "H2O":
                return self.flow_vol * self.params.dens_mass
            elif j == "S_ALK":
                # Convert moles of alkalinity to mass of C assuming all is HCO3-
                return (
                    self.flow_vol
                    * self.alkalinity
                    * (12 * pyo.units.kg / pyo.units.kmol)
                )
            else:
                return self.flow_vol * self.conc_mass_comp[j]

        self.material_flow_expression = pyo.Expression(
            self.component_list,
            rule=material_flow_expression,
            doc="Material flow terms",
        )

        def enthalpy_flow_expression(self):
            return (
                self.flow_vol
                * self.params.dens_mass
                * self.params.cp_mass
                * (self.temperature - self.params.temperature_ref)
            )

        self.enthalpy_flow_expression = pyo.Expression(
            rule=enthalpy_flow_expression, doc="Enthalpy flow term"
        )

        def material_density_expression(self, j):
            if j == "H2O":
                return self.params.dens_mass
            elif j == "S_ALK":
                # Convert moles of alkalinity to mass of C assuming all is HCO3-
                return self.alkalinity * (12 * pyo.units.kg / pyo.units.kmol)
            else:
                return self.conc_mass_comp[j]

        self.material_density_expression = pyo.Expression(
            self.component_list,
            rule=material_density_expression,
            doc="Material density terms",
        )

        def energy_density_expression(self):
            return (
                self.params.dens_mass
                * self.params.cp_mass
                * (self.temperature - self.params.temperature_ref)
            )

        self.energy_density_expression = pyo.Expression(
            rule=energy_density_expression, doc="Energy density term"
        )

    def get_material_flow_terms(self, p, j):
        return self.material_flow_expression[j]

    def get_enthalpy_flow_terms(self, p):
        return self.enthalpy_flow_expression

    def get_material_density_terms(self, p, j):
        return self.material_density_expression[j]

    def get_energy_density_terms(self, p):
        return self.energy_density_expression

    def default_material_balance_type(self):
        return MaterialBalanceType.componentPhase

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def define_state_vars(b):
        return {
            "flow_vol": b.flow_vol,
            "alkalinity": b.alkalinity,
            "conc_mass_comp": b.conc_mass_comp,
            "temperature": b.temperature,
            "pressure": b.pressure,
        }

    def define_display_vars(b):
        return {
            "Volumetric Flowrate": b.flow_vol,
            "Molar Alkalinity": b.alkalinity,
            "Mass Concentration": b.conc_mass_comp,
            "Temperature": b.temperature,
            "Pressure": b.pressure,
        }

    def get_material_flow_basis(b):
        return MaterialFlowBasis.mass

    def calculate_scaling_factors(self):
        # Get default scale factors and do calculations from base classes
        super().calculate_scaling_factors()

        # No constraints in this model as yet, just need to set scaling factors
        # for expressions
        sf_F = iscale.get_scaling_factor(self.flow_vol, default=1e2, warning=True)
        sf_T = iscale.get_scaling_factor(self.temperature, default=1e-2, warning=True)

        # Mass flow and density terms
        for j in self.component_list:
            if j == "H2O":
                sf_C = pyo.value(1 / self.params.dens_mass)
            elif j == "S_ALK":
                sf_C = 1e-1 * iscale.get_scaling_factor(
                    self.alkalinity, default=1, warning=True
                )
            else:
                sf_C = iscale.get_scaling_factor(
                    self.conc_mass_comp[j], default=1e2, warning=True
                )

            iscale.set_scaling_factor(self.material_flow_expression[j], sf_F * sf_C)
            iscale.set_scaling_factor(self.material_density_expression[j], sf_C)

        # Enthalpy and energy terms
        sf_rho_cp = pyo.value(1 / (self.params.dens_mass * self.params.cp_mass))
        iscale.set_scaling_factor(
            self.enthalpy_flow_expression, sf_F * sf_rho_cp * sf_T
        )
        iscale.set_scaling_factor(self.energy_density_expression, sf_rho_cp * sf_T)
