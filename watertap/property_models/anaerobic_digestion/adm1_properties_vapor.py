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
    VaporPhase,
    Component,
    Solute,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.constants import Constants
from idaes.core.util.initialization import fix_state_vars, revert_state_vars
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale

# Some more information about this module
__author__ = "Alejandro Garciadiego, Xinhong Liu"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("ADM1_vaporParameterBlock")
class ADM1_vaporParameterData(PhysicalParameterBlock):
    """
    Property Parameter Block Class
    """

    def build(self):
        """
        Callable method for Block construction.
        """
        super().build()

        self._state_block_class = ADM1_vaporStateBlock

        # Add Phase objects
        self.Vap = VaporPhase()

        # Add Component objects
        self.H2O = Component()

        # All soluble components on kg COD/m^3 basis
        self.S_h2 = Solute(doc="Hydrogen gas")
        self.S_ch4 = Solute(doc="Methane gas")
        self.S_co2 = Solute(doc="Carbon dioxide")

        # Heat capacity of water
        self.cp_mass = pyo.Param(
            mutable=False,
            initialize=1.996,
            doc="Specific heat capacity of water",
            units=pyo.units.J / pyo.units.kg / pyo.units.K,
        )
        # Density of water
        self.dens_mass = pyo.Param(
            mutable=False,
            # initialize=0.927613356,
            initialize=0.01,
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
                "pressure_sat": {"method": None},
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


class _ADM1_vaporStateBlock(StateBlock):
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
            flow_mol_comp : value at which to initialize component flows (default=None)
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
                         initialization by calling the
                         release_state method

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
        Method to release state variables fixed during initialization.

        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state=True.
            outlvl : sets output level of logging
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")

        if flags is None:
            return
        # Unfix state variables
        revert_state_vars(self, flags)
        init_log.info("State Released.")


@declare_process_block_class("ADM1_vaporStateBlock", block_class=_ADM1_vaporStateBlock)
class ADM1_vaporStateBlockData(StateBlockData):
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
            initialize=1,
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

        Comp_dict = {"S_ch4": 1.6256, "S_co2": 0.01415 * 12, "S_h2": 1e-5}

        self.conc_mass_comp = pyo.Var(
            self.params.solute_set,
            domain=pyo.NonNegativeReals,
            initialize=Comp_dict,
            doc="Component mass concentrations",
            units=pyo.units.kg / pyo.units.m**3,
        )

        init = {"S_ch4": 65077, "S_co2": 36255, "S_h2": 1.639, "H2O": 5570}

        self.pressure_sat = pyo.Var(
            self.params.component_list,
            domain=pyo.NonNegativeReals,
            initialize=init,
            doc="Component pressure",
            units=pyo.units.Pa,
        )

        def pressure_sat_rule(b, j):
            if j == "S_h2":
                return b.pressure_sat[j] == pyo.units.convert(
                    b.conc_mass_comp[j]
                    * (1000 * pyo.units.g / pyo.units.kg)
                    * Constants.gas_constant
                    * b.temperature
                    / (16 * pyo.units.g / pyo.units.mole),
                    to_units=pyo.units.Pa,
                )
            elif j == "S_ch4":
                return b.pressure_sat[j] == pyo.units.convert(
                    b.conc_mass_comp[j]
                    * (1000 * pyo.units.g / pyo.units.kg)
                    * Constants.gas_constant
                    * b.temperature
                    / (64 * pyo.units.g / pyo.units.mole),
                    to_units=pyo.units.Pa,
                )
            elif j == "H2O":
                return pyo.log(b.pressure_sat[j] / pyo.units.Pa) == (
                    pyo.log(0.0313)
                    + 5290
                    * pyo.units.K
                    * ((1 / b.params.temperature_ref) - (1 / b.temperature))
                    + pyo.log(1e5)
                )
            elif j == "S_co2":
                return b.pressure_sat[j] == pyo.units.convert(
                    b.conc_mass_comp[j]
                    * (1000 * pyo.units.g / pyo.units.kg)
                    * Constants.gas_constant
                    * b.temperature
                    / (12 * pyo.units.g / pyo.units.mole),
                    to_units=pyo.units.Pa,
                )
            else:
                raise Exception("Vapor component not implemented.")

        self._pressure_sat = pyo.Constraint(
            self.params.component_list,
            rule=pressure_sat_rule,
            doc="Saturation pressure for components",
        )

        def material_flow_expression(self, j):
            if j == "H2O":
                return self.flow_vol * self.params.dens_mass
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

        iscale.set_scaling_factor(self.flow_vol, 1e5)
        iscale.set_scaling_factor(self.temperature, 1e-1)
        iscale.set_scaling_factor(self.pressure, 1e-3)
        iscale.set_scaling_factor(self.conc_mass_comp, 1e2)
        iscale.set_scaling_factor(self.conc_mass_comp["S_h2"], 1e3)
        iscale.set_scaling_factor(self.pressure_sat, 1e-3)
        iscale.set_scaling_factor(self.pressure_sat["S_h2"], 1e-2)

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

        for t, v in self._pressure_sat.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    self.pressure_sat,
                    default=1,
                    warning=True,
                ),
            )
