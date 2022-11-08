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
Thermophysical property package to be used in conjunction with ADM1 reactions.
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

# Some more information about this module
__author__ = "Adam Atia"
# Using Andrew Lee's formulation of ASM1 as a template

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("ADM1ParameterBlock")
class ADM1ParameterData(PhysicalParameterBlock):
    """
    Property Parameter Block Class
    """

    def build(self):
        """
        Callable method for Block construction.
        """
        super().build()

        self._state_block_class = ADM1StateBlock

        # Add Phase objects
        self.Liq = LiquidPhase()

        # Add Component objects
        self.H2O = Solvent()

        # All soluble components on kg COD/m^3 basis
        self.S_su = Solute(doc="Monosaccharides")
        self.S_aa = Solute(doc="Amino acids")
        self.S_fa = Solute(doc="Long chain fatty acids")
        self.S_va = Solute(doc="Total valerate")
        self.S_bu = Solute(doc="Total butyrate")
        self.S_pro = Solute(doc="Total propionate")
        self.S_ac = Solute(doc="Total acetate")
        self.S_h2 = Solute(doc="Hydrogen gas")
        self.S_ch4 = Solute(doc="Methane gas")
        self.S_IC = Solute(doc="Inorganic carbon")
        self.S_IN = Solute(doc="Inorganic nitrogen")
        self.S_I = Solute(doc="Soluble inerts")

        self.X_c = Solute(doc="Composites")
        self.X_ch = Solute(doc="Carbohydrates")
        self.X_pr = Solute(doc="Proteins")
        self.X_li = Solute(doc="Lipids")
        self.X_su = Solute(doc="Sugar degraders")
        self.X_aa = Solute(doc="Amino acid degraders")
        self.X_fa = Solute(doc="Long chain fatty acid (LCFA) degraders")
        self.X_c4 = Solute(doc="Valerate and butyrate degraders")
        self.X_pro = Solute(doc="Propionate degraders")
        self.X_ac = Solute(doc="Acetate degraders")
        self.X_h2 = Solute(doc="Hydrogen degraders")
        self.X_I = Solute(doc="Particulate inerts")

        # TODO: Additional components not in Table but referred to in text
        self.S_cat = Component(doc="Total cation equivalents concentration")
        self.S_an = Component(doc="Total anion equivalents concentration")

        # TODO: checkpoint

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
                "anions": {"method": None},
                "cations": {"method": None},
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


class _ADM1StateBlock(StateBlock):
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
                     were not provided at the unit model level, the
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
                         - True - states variables are not unfixed, and
                                 a dict of returned containing flags for
                                 which states were fixed during
                                 initialization.
                        - False - state variables are unfixed after
                                 initialization by calling the
                                 release_state method

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
        Method to release state variables fixed during initialization.

        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state=True.
            outlvl : sets output level of logging
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="properties")

        if flags is None:
            return
        # Unfix state variables
        revert_state_vars(blk, flags)
        init_log.info("State Released.")


@declare_process_block_class("ADM1StateBlock", block_class=_ADM1StateBlock)
class ADM1StateBlockData(StateBlockData):
    """
    StateBlock for calculating thermophysical properties associated with the ADM1
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
            bounds=(298.15, 323.15),
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
        self.anions = pyo.Var(
            domain=pyo.NonNegativeReals,
            initialize=0.1,
            doc="Ions in molar concentration",
            units=pyo.units.kmol / pyo.units.m**3,
        )
        self.cations = pyo.Var(
            domain=pyo.NonNegativeReals,
            initialize=0.1,
            doc="Ions in molar concentration",
            units=pyo.units.kmol / pyo.units.m**3,
        )

    def get_material_flow_terms(b, p, j):
        if j == "H2O":
            return b.flow_vol * b.params.dens_mass
        elif j == "S_cat":
            # Convert moles of anions to mass assuming all is CL-
            return b.flow_vol * b.anions * (35 * pyo.units.kg / pyo.units.kmol)
        elif j == "S_an":
            # Convert moles of cations to mass assuming all is Na+
            return b.flow_vol * b.cations * (23 * pyo.units.kg / pyo.units.kmol)
        else:
            return b.flow_vol * b.conc_mass_comp[j]

    def get_enthalpy_flow_terms(b, p):
        return (
            b.flow_vol
            * b.params.dens_mass
            * b.params.cp_mass
            * (b.temperature - b.params.temperature_ref)
        )

    def get_material_density_terms(b, p, j):
        if j == "H2O":
            return b.params.dens_mass
        else:
            return b.conc_mass_comp[j]

    def get_energy_density_terms(b, p):
        return (
            b.params.dens_mass
            * b.params.cp_mass
            * (b.temperature - b.params.temperature_ref)
        )

    def default_material_balance_type(self):
        return MaterialBalanceType.componentPhase

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def define_state_vars(b):
        return {
            "flow_vol": b.flow_vol,
            "anions": b.anions,
            "cations": b.cations,
            "conc_mass_comp": b.conc_mass_comp,
            "temperature": b.temperature,
            "pressure": b.pressure,
        }

    def define_display_vars(b):
        return {
            "Volumetric Flowrate": b.flow_vol,
            "Molar anions": b.anions,
            "Molar cations": b.cations,
            "Mass Concentration": b.conc_mass_comp,
            "Temperature": b.temperature,
            "Pressure": b.pressure,
        }

    def get_material_flow_basis(b):
        return MaterialFlowBasis.mass
