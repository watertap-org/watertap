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
Thermophysical property package to be used in conjunction with modified ASM2d reactions.

Reference:
[1] X. Flores-Alsina, K. Solon, C.K. Mbamba, S. Tait, K.V. Gernaey, U. Jeppsson, D.J. Batstone,
Modelling phosphorus (P), sulfur (S) and iron (Fe) interactions for dynamic simulations of anaerobic digestion processes,
Water Research. 95 (2016) 370-382. https://www.sciencedirect.com/science/article/pii/S0043135416301397

[2] K. Solon, X. Flores-Alsina, C. Kazadi Mbamba, D. Ikumi, E.I.P. Volcke, C. Vaneeckhaute, G. Ekama,
P.A. Vanrolleghem, D.J. Batstone, K.V. Gernaey, U. Jeppsson, Plant-wide modelling of phosphorus transformations in
wastewater treatment systems: Impacts of control and operational strategies, Water Research. 113 (2017) 97-110
https://www.sciencedirect.com/science/article/pii/S0043135417300829

"""

# Import Pyomo libraries
import pyomo.environ as pyo

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    MaterialFlowBasis,
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
from idaes.core.scaling import CustomScalerBase
from idaes.core.base.property_base import PhysicalParameterBlock

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

        # Effluent Quality Index (EQI) parameters [2]
        self.i_NSF = pyo.Var(
            initialize=0.03352,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="N content of fermentable substrate, S_F, [kg N/kg COD]",
        )
        self.i_NSI = pyo.Var(
            initialize=0.06003,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="N content of inert soluble COD S_I, [kg N/kg COD]",
        )
        self.i_NXI = pyo.Var(
            initialize=0.06003,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="N content of inert particulate COD X_I, [kg N/kg COD]",
        )
        self.i_NXS = pyo.Var(
            initialize=0.03352,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="N content of slowly biodegradable substrate X_S, [kg N/kg COD]",
        )
        self.i_NBM = pyo.Var(
            initialize=0.08615,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="N content of biomass, X_H, X_PAO, X_AUT, [kg N/kg COD]",
        )
        self.f_SI = pyo.Var(
            initialize=0.00,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Production of S_I in hydrolysis, [kg COD/kg COD]",
        )
        self.f_XIH = pyo.Var(
            initialize=0.1,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Fraction of inert COD generated in lysis of X_H, [kg COD/kg COD]",
        )
        self.f_XIP = pyo.Var(
            initialize=0.1,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Fraction of inert COD generated in lysis of X_PAO and X_PHA, [kg COD/kg COD]",
        )
        self.f_XIA = pyo.Var(
            initialize=0.1,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="Fraction of inert COD generated in lysis of X_AUT, [kg COD/kg COD]",
        )
        self.i_PSF = pyo.Var(
            initialize=0.00559,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="P content of fermentable substrate, S_F, [kg P/kg COD]",
        )
        self.i_PSI = pyo.Var(
            initialize=0.00649,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="P content of inert soluble COD S_I, [kg P/kg COD]",
        )
        self.i_PXI = pyo.Var(
            initialize=0.00649,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="P content of inert particulate COD X_I, [kg P/kg COD]",
        )
        self.i_PXS = pyo.Var(
            initialize=0.00559,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="P content of slowly biodegradable substrate X_S, [kg P/kg COD]",
        )
        self.i_PBM = pyo.Var(
            initialize=0.02154,
            units=pyo.units.dimensionless,
            domain=pyo.NonNegativeReals,
            doc="P content of biomass, X_H, X_PAO, X_AUT, [kg P/kg COD]",
        )
        self.BOD5_factor = pyo.Param(
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
                "VSS": {"method": "_VSS"},
                "ISS": {"method": "_ISS"},
                "TSS": {"method": "_TSS"},
                "COD": {"method": "_COD"},
                "TKN": {"method": "_TKN"},
                "SNOX": {"method": "_SNOX"},
                "BOD5": {"method": "_BOD5"},
                "SP_organic": {"method": "_SP_organic"},
                "SP_inorganic": {"method": "_SP_inorganic"},
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


class ModifiedASM2dPropertiesScaler(CustomScalerBase):
    """
    Scaler for the Activated Sludge Model No.2d property package.
    Flow and temperature are scaled by the default value (if no user input provided), and
    pressure is scaled assuming an order of magnitude of 1e5 Pa.
    """

    UNIT_SCALING_FACTORS = {
        # "QuantityName: (reference units, scaling factor)
        "Pressure": (pyo.units.Pa, 1e-5),
    }

    DEFAULT_SCALING_FACTORS = {
        "flow_vol": 1e1,
        "temperature": 1e-2,
    }

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        self.scale_variable_by_default(model.temperature, overwrite=overwrite)
        self.scale_variable_by_default(model.flow_vol, overwrite=overwrite)
        self.scale_variable_by_units(model.pressure, overwrite=overwrite)

    # There are currently no constraints in this model
    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        pass


class _ModifiedASM2dStateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """

    default_scaler = ModifiedASM2dPropertiesScaler

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

        # TODO: X_SRB not included yet in biomass term summation
        def _VSS(self):
            vss = (
                self.conc_mass_comp["X_I"] / self.params.CODtoVSS_XI
                + self.conc_mass_comp["X_S"] / self.params.CODtoVSS_XS
                + (
                    self.conc_mass_comp["X_H"]
                    + self.conc_mass_comp["X_PAO"]
                    + self.conc_mass_comp["X_AUT"]
                )
                / self.params.CODtoVSS_XBM
                + self.conc_mass_comp["X_PHA"] / self.params.CODtoVSS_XPHA
            )
            return vss

        self.VSS = pyo.Expression(rule=_VSS, doc="Volatile suspended solids")

        def _ISS(self):
            iss = (
                self.params.f_ISS_BM
                * (
                    self.conc_mass_comp["X_H"]
                    + self.conc_mass_comp["X_PAO"]
                    + self.conc_mass_comp["X_AUT"]
                )
                / self.params.CODtoVSS_XBM
                + self.params.ISS_P * self.conc_mass_comp["X_PP"]
            )
            return iss

        self.ISS = pyo.Expression(rule=_ISS, doc="Inorganic suspended solids")

        def _TSS(self):
            tss = self.VSS + self.ISS
            return tss

        self.TSS = pyo.Expression(rule=_TSS, doc="Total suspended solids")

        def _COD(self):
            cod = (
                self.conc_mass_comp["S_F"]
                + self.conc_mass_comp["S_A"]
                + self.conc_mass_comp["S_I"]
                + self.conc_mass_comp["X_I"]
                + self.conc_mass_comp["X_S"]
                + self.conc_mass_comp["X_H"]
                + self.conc_mass_comp["X_PAO"]
                + self.conc_mass_comp["X_PHA"]
                + self.conc_mass_comp["X_AUT"]
            )
            return cod

        self.COD = pyo.Expression(rule=_COD, doc="Chemical oxygen demand")

        def _TKN(self):
            tkn = (
                self.conc_mass_comp["S_NH4"]
                + self.params.i_NSF * self.conc_mass_comp["S_F"]
                + self.params.i_NSI * self.conc_mass_comp["S_I"]
                + self.params.i_NXI * self.conc_mass_comp["X_I"]
                + self.params.i_NXS * self.conc_mass_comp["X_S"]
                + self.params.i_NBM
                * (
                    self.conc_mass_comp["X_H"]
                    + self.conc_mass_comp["X_PAO"]
                    + self.conc_mass_comp["X_AUT"]
                )
            )
            return tkn

        self.TKN = pyo.Expression(rule=_TKN, doc="Kjeldahl nitrogen")

        def _SNOX(self):
            snox = self.conc_mass_comp["S_NO3"]

            return snox

        self.SNOX = pyo.Expression(rule=_SNOX, doc="Nitrogen oxide")

        def _BOD5(self, i):
            bod5 = (
                self.conc_mass_comp["S_F"]
                + self.conc_mass_comp["S_A"]
                + (1 - self.params.f_SI) * self.conc_mass_comp["X_S"]
                + (1 - self.params.f_XIH) * self.conc_mass_comp["X_H"]
                + (1 - self.params.f_XIP)
                * (self.conc_mass_comp["X_PAO"] + self.conc_mass_comp["X_PHA"])
                + (1 - self.params.f_XIA) * self.conc_mass_comp["X_AUT"]
            )

            return self.params.BOD5_factor[i] * bod5

        self.BOD5 = pyo.Expression(
            ["raw", "effluent"], rule=_BOD5, doc="Five-day biological oxygen demand"
        )

        def _SP_organic(self):
            sp_organic = (
                self.conc_mass_comp["X_PP"]
                + self.params.i_PSF * self.conc_mass_comp["S_F"]
                + self.params.i_PSI * self.conc_mass_comp["S_I"]
                + self.params.i_PXI * self.conc_mass_comp["X_I"]
                + self.params.i_PXS * self.conc_mass_comp["X_S"]
                + self.params.i_PBM
                * (
                    self.conc_mass_comp["X_H"]
                    + self.conc_mass_comp["X_PAO"]
                    + self.conc_mass_comp["X_AUT"]
                )
            )
            return sp_organic

        self.SP_organic = pyo.Expression(rule=_SP_organic, doc="Organic phosphorus")

        def _SP_inorganic(self):
            sp_inorganic = self.conc_mass_comp["S_PO4"]
            return sp_inorganic

        self.SP_inorganic = pyo.Expression(
            rule=_SP_inorganic, doc="Inorganic phosphorus"
        )

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
