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
__author__ = "Andrew Lee, Adam Atia, Xinhong Liu"

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

        # Create sets for use across ASM models and associated unit models (e.g., thickener, dewaterer)
        self.non_particulate_component_set = pyo.Set(
            initialize=[
                "S_I",
                "S_S",
                "S_O",
                "S_NO",
                "S_NH",
                "S_ND",
                "H2O",
                "S_ALK",
            ]
        )
        self.particulate_component_set = pyo.Set(
            initialize=["X_I", "X_S", "X_P", "X_BH", "X_BA", "X_ND"]
        )
        self.tss_component_set = pyo.Set(
            initialize=["X_I", "X_S", "X_P", "X_BH", "X_BA"]
        )

        # Heat capacity of water
        self.cp_mass = pyo.Param(
            initialize=4182,
            doc="Specific heat capacity of water",
            units=pyo.units.J / pyo.units.kg / pyo.units.K,
        )
        # Density of water
        self.dens_mass = pyo.Param(
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
        self.f_p = pyo.Var(
            initialize=0.08,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Fraction of biomass yielding particulate products, f_p",
        )
        self.i_xb = pyo.Var(
            initialize=0.08,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Mass fraction of N per COD in biomass, i_xb",
        )
        self.i_xp = pyo.Var(
            initialize=0.06,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Mass fraction of N per COD in particulates, i_xp",
        )
        self.COD_to_SS = pyo.Var(
            initialize=0.75,
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Conversion factor applied for TSS calculation",
        )
        self.BOD5_factor = pyo.Var(
            ["raw", "effluent"],
            initialize={"raw": 0.65, "effluent": 0.25},
            units=pyo.units.dimensionless,
            domain=pyo.PositiveReals,
            doc="Conversion factor for BOD5",
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
                "alkalinity": {"method": None},
                "TSS": {"method": "_TSS"},
                "BOD5": {"method": "_BOD5"},
                "TKN": {"method": "_TKN"},
                "Total_N": {"method": "_Total_N"},
                "COD": {"method": "_COD"},
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
            alkalinity: value of alkalinity expressed as molar concentration
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

        def _TSS(self):
            tss = (
                self.conc_mass_comp["X_S"]
                + self.conc_mass_comp["X_I"]
                + self.conc_mass_comp["X_BH"]
                + self.conc_mass_comp["X_BA"]
                + self.conc_mass_comp["X_P"]
            )
            return self.params.COD_to_SS * tss

        self.TSS = pyo.Expression(
            rule=_TSS,
            doc="Total suspended solids (TSS)",
        )

        def _BOD5(self, i):
            bod5 = (
                self.conc_mass_comp["X_S"]
                + self.conc_mass_comp["X_S"]
                + (1 - self.params.f_p)
                * (self.conc_mass_comp["X_BH"] + self.conc_mass_comp["X_BA"])
            )
            # TODO: 0.25 should be a parameter instead as it changes by influent/effluent
            return self.params.BOD5_factor[i] * bod5

        self.BOD5 = pyo.Expression(
            ["raw", "effluent"],
            rule=_BOD5,
            doc="Five-day Biological Oxygen Demand (BOD5)",
        )

        def _COD(self):
            cod = (
                self.conc_mass_comp["S_S"]
                + self.conc_mass_comp["S_I"]
                + self.conc_mass_comp["X_S"]
                + self.conc_mass_comp["X_S"]
                + self.conc_mass_comp["X_I"]
                + self.conc_mass_comp["X_BH"]
                + self.conc_mass_comp["X_BA"]
                + self.conc_mass_comp["X_P"]
            )
            return cod

        self.COD = pyo.Expression(
            rule=_COD,
            doc="Chemical Oxygen Demand",
        )

        def _TKN(self):
            tkn = (
                self.conc_mass_comp["S_NH"]
                + self.conc_mass_comp["S_ND"]
                + self.conc_mass_comp["X_ND"]
                + self.params.i_xb
                * (self.conc_mass_comp["X_BH"] + self.conc_mass_comp["X_BA"])
                + self.params.i_xp
                * (self.conc_mass_comp["X_P"] + self.conc_mass_comp["X_I"])
            )
            return tkn

        self.TKN = pyo.Expression(
            rule=_TKN,
            doc="Total Kjeldahl Nitrogen",
        )

        def _Total_N(self):
            totaln = self.TKN + self.conc_mass_comp["S_NO"]
            return totaln

        self.Total_N = pyo.Expression(
            rule=_Total_N,
            doc="Total Nitrogen",
        )

        iscale.set_scaling_factor(self.flow_vol, 1e1)
        iscale.set_scaling_factor(self.temperature, 1e-1)
        iscale.set_scaling_factor(self.pressure, 1e-6)
        iscale.set_scaling_factor(self.conc_mass_comp, 1e1)

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

    def define_state_vars(self):
        return {
            "flow_vol": self.flow_vol,
            "alkalinity": self.alkalinity,
            "conc_mass_comp": self.conc_mass_comp,
            "temperature": self.temperature,
            "pressure": self.pressure,
        }

    def define_display_vars(self):
        return {
            "Volumetric Flowrate": self.flow_vol,
            "Molar Alkalinity": self.alkalinity,
            "Mass Concentration": self.conc_mass_comp,
            "Temperature": self.temperature,
            "Pressure": self.pressure,
        }

    def get_material_flow_basis(self):
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
