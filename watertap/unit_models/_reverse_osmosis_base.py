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

from copy import deepcopy
from enum import Enum, auto
from pyomo.environ import (
    Block,
    exp,
    NonNegativeReals,
    Param,
    Suffix,
    Var,
    units as pyunits,
)
from pyomo.common.collections import ComponentSet
from pyomo.common.config import ConfigBlock, ConfigValue, In
from idaes.core import (
    UnitModelBlockData,
    useDefault,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from idaes.core.util import scaling as iscale
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.tables import create_stream_table_dataframe
import idaes.logger as idaeslog

from watertap.core.membrane_channel_base import validate_membrane_config_args

_log = idaeslog.getLogger(__name__)


class ConcentrationPolarizationType(Enum):
    none = auto()  # simplified assumption: no concentration polarization
    fixed = (
        auto()
    )  # simplified assumption: concentration polarization modulus is a user specified value
    calculated = (
        auto()
    )  # calculate concentration polarization (concentration at membrane interface)


class MassTransferCoefficient(Enum):
    none = (
        auto()
    )  # mass transfer coefficient not utilized for concentration polarization effect
    fixed = auto()  # mass transfer coefficient is a user specified value
    calculated = auto()  # mass transfer coefficient is calculated
    # TODO: add option for users to define their own relationship?


class PressureChangeType(Enum):
    fixed_per_stage = (
        auto()
    )  # pressure drop across membrane channel is a user-specified value
    fixed_per_unit_length = (
        auto()
    )  # pressure drop per unit length across membrane channel is a user-specified value
    calculated = auto()  # pressure drop across membrane channel is calculated


class _ReverseOsmosisBaseData(UnitModelBlockData):
    """
    Reverse Osmosis base class
    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not.
    **default** = False. Membrane units do not yet support dynamic
    behavior.""",
        ),
    )

    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False. Membrane units do not have defined volume, thus
    this must be False.""",
        ),
    )

    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
    **default** - useDefault.
    **Valid values:** {
    **useDefault** - use default package from parent model or flowsheet,
    **PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )

    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
    and used when constructing these.
    **default** - None.
    **Valid values:** {
    see property package for documentation.}""",
        ),
    )

    CONFIG.declare(
        "material_balance_type",
        ConfigValue(
            default=MaterialBalanceType.useDefault,
            domain=In(MaterialBalanceType),
            description="Material balance construction flag",
            doc="""Indicates what type of mass balance should be constructed,
    **default** - MaterialBalanceType.useDefault.
    **Valid values:** {
    **MaterialBalanceType.useDefault - refer to property package for default
    balance type
    **MaterialBalanceType.none** - exclude material balances,
    **MaterialBalanceType.componentPhase** - use phase component balances,
    **MaterialBalanceType.componentTotal** - use total component balances,
    **MaterialBalanceType.elementTotal** - use total element balances,
    **MaterialBalanceType.total** - use total material balance.}""",
        ),
    )

    CONFIG.declare(
        "energy_balance_type",
        ConfigValue(
            default=EnergyBalanceType.useDefault,
            domain=In(EnergyBalanceType),
            description="Energy balance construction flag",
            doc="""Indicates what type of energy balance should be constructed.
    **default** - EnergyBalanceType.useDefault.
    **Valid values:** {
    **EnergyBalanceType.useDefault - refer to property package for default
    balance type
    **EnergyBalanceType.none** - exclude energy balances,
    **EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
    **EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
    **EnergyBalanceType.energyTotal** - single energy balance for material,
    **EnergyBalanceType.energyPhase** - energy balances for each phase.}""",
        ),
    )

    CONFIG.declare(
        "momentum_balance_type",
        ConfigValue(
            default=MomentumBalanceType.pressureTotal,
            domain=In(MomentumBalanceType),
            description="Momentum balance construction flag",
            doc="""Indicates what type of momentum balance should be constructed,
    **default** - MomentumBalanceType.pressureTotal.
    **Valid values:** {
    **MomentumBalanceType.none** - exclude momentum balances,
    **MomentumBalanceType.pressureTotal** - single pressure balance for material,
    **MomentumBalanceType.pressurePhase** - pressure balances for each phase,
    **MomentumBalanceType.momentumTotal** - single momentum balance for material,
    **MomentumBalanceType.momentumPhase** - momentum balances for each phase.}""",
        ),
    )

    CONFIG.declare(
        "concentration_polarization_type",
        ConfigValue(
            default=ConcentrationPolarizationType.calculated,
            domain=In(ConcentrationPolarizationType),
            description="External concentration polarization effect in RO",
            doc="""
            Options to account for concentration polarization.

            **default** - ``ConcentrationPolarizationType.calculated``

        .. csv-table::
            :header: "Configuration Options", "Description"

            "``ConcentrationPolarizationType.none``", "Simplifying assumption to ignore concentration polarization"
            "``ConcentrationPolarizationType.fixed``", "Specify an estimated value for the concentration polarization modulus"
            "``ConcentrationPolarizationType.calculated``", "Allow model to perform calculation of membrane-interface concentration"
        """,
        ),
    )

    CONFIG.declare(
        "mass_transfer_coefficient",
        ConfigValue(
            default=MassTransferCoefficient.calculated,
            domain=In(MassTransferCoefficient),
            description="Mass transfer coefficient in RO feed channel",
            doc="""
            Options to account for mass transfer coefficient.

            **default** - ``MassTransferCoefficient.calculated``

        .. csv-table::
            :header: "Configuration Options", "Description"

            "``MassTransferCoefficient.none``", "Mass transfer coefficient not used in calculations"
            "``MassTransferCoefficient.fixed``", "Specify an estimated value for the mass transfer coefficient in the feed channel"
            "``MassTransferCoefficient.calculated``", "Allow model to perform calculation of mass transfer coefficient"
        """,
        ),
    )

    CONFIG.declare(
        "has_pressure_change",
        ConfigValue(
            default=False,
            domain=In([True, False]),
            description="Pressure change term construction flag",
            doc="""Indicates whether terms for pressure change should be
    constructed,
    **default** - False.
    **Valid values:** {
    **True** - include pressure change terms,
    **False** - exclude pressure change terms.}""",
        ),
    )

    CONFIG.declare(
        "pressure_change_type",
        ConfigValue(
            default=PressureChangeType.fixed_per_stage,
            domain=In(PressureChangeType),
            description="Pressure change term construction flag",
            doc="""
        Indicates what type of pressure change calculation will be made. To use any of the
        ``pressure_change_type`` options to account for pressure drop, the configuration keyword
        ``has_pressure_change`` must also be set to ``True``. Also, if a value is specified for pressure
        change, it should be negative to represent pressure drop.

        **default** - ``PressureChangeType.fixed_per_stage``


    .. csv-table::
        :header: "Configuration Options", "Description"

        "``PressureChangeType.fixed_per_stage``", "Specify an estimated value for pressure drop across the membrane feed channel"
        "``PressureChangeType.fixed_per_unit_length``", "Specify an estimated value for pressure drop per unit length across the membrane feed channel"
        "``PressureChangeType.calculated``", "Allow model to perform calculation of pressure drop across the membrane feed channel"

""",
        ),
    )

    CONFIG.declare(
        "has_full_reporting",
        ConfigValue(
            default=False,
            domain=In([True, False]),
            description="Level of reporting results",
            doc="""Level of reporting results.
            **default** - False.
            **Valid values:** {
            **False** - include minimal reporting of results,
            **True** - report additional properties of interest that aren't constructed by
            the unit model by default. Also, report averaged expression values""",
        ),
    )

    def build(self):
        """
        Common variables and constraints for an RO unit model
        """
        super().build()

        if len(self.config.property_package.solvent_set) > 1:
            raise ConfigurationError(
                "Membrane models only support one solvent component,"
                "the provided property package has specified {} solvent components".format(
                    len(self.config.property_package.solvent_set)
                )
            )

        if len(self.config.property_package.phase_list) > 1 or "Liq" not in [
            p for p in self.config.property_package.phase_list
        ]:
            raise ConfigurationError(
                "Membrane models only support one liquid phase ['Liq'],"
                "the property package has specified the following phases {}".format(
                    [p for p in self.config.property_package.phase_list]
                )
            )

        validate_membrane_config_args(self)

        self._add_membrane_channel()

        self.membrane_channel.add_geometry()

        self.membrane_channel.add_state_blocks(has_phase_equilibrium=False)

        self.membrane_channel.add_material_balances(
            balance_type=self.config.material_balance_type, has_mass_transfer=True
        )

        self.membrane_channel.add_energy_balances(
            balance_type=self.config.energy_balance_type, has_enthalpy_transfer=True
        )

        self.membrane_channel.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            pressure_change_type=self.config.pressure_change_type,
            has_pressure_change=self.config.has_pressure_change,
        )

        self.membrane_channel.add_mass_transfer()

        self.membrane_channel.add_isothermal_conditions()

        self.membrane_channel.add_volumetric_flowrate_balance()

        self.membrane_channel.add_flux_balance()

        self.membrane_channel.add_concentration_polarization(
            concentration_polarization_type=self.config.concentration_polarization_type,
            mass_transfer_coefficient=self.config.mass_transfer_coefficient,
        )

        self.membrane_channel.add_recovery_vol_phase()

        self.membrane_channel.add_expressions()

        self.membrane_channel.apply_transformation()

        # Add Ports
        self.add_inlet_port(name="inlet", block=self.membrane_channel)
        self.add_outlet_port(name="retentate", block=self.membrane_channel)
        self.add_port(name="permeate", block=self.membrane_channel.mixed_permeate)

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Feed Inlet": self.inlet,
                "Feed Outlet": self.retentate,
                "Permeate Outlet": self.permeate,
            },
            time_point=time_point,
        )

    def _get_state_args(
        self, source, mixed_permeate_properties, initialize_guess, state_args
    ):
        """
        Arguments:
            source : property model containing inlet feed
            mixed_permeate_properties : mixed permeate property block
            initialize_guess : a dict of guesses for solvent_recovery, solute_recovery,
                               and cp_modulus. These guesses offset the initial values
                               for the retentate, permeate, and membrane interface
                               state blocks from the inlet feed
                               (default =
                               {'deltaP': -1e4,
                               'solvent_recovery': 0.5,
                               'solute_recovery': 0.01,
                               'cp_modulus': 1.1})
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for the inlet
                         feed side state block (see documentation of the specific
                         property package).
        """

        # assumptions
        if initialize_guess is None:
            initialize_guess = {}
        if "deltaP" not in initialize_guess:
            initialize_guess["deltaP"] = -1e4
        if "solvent_recovery" not in initialize_guess:
            initialize_guess["solvent_recovery"] = 0.5
        if "solute_recovery" not in initialize_guess:
            initialize_guess["solute_recovery"] = 0.01
        if "cp_modulus" not in initialize_guess:
            if (
                self.config.concentration_polarization_type
                != ConcentrationPolarizationType.none
            ):
                initialize_guess["cp_modulus"] = 1.1
            else:
                initialize_guess["cp_modulus"] = 1

        if state_args is None:
            state_args = {}
            state_dict = source.define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        if "flow_mass_phase_comp" not in state_args.keys():
            raise ConfigurationError(
                f"{self.__class__.__name__} initialization routine expects "
                "flow_mass_phase_comp as a state variable. Check "
                "that the property package supports this state "
                "variable or that the state_args provided to the "
                "initialize call includes this state variable"
            )

        # slightly modify initial values for other state blocks
        state_args_retentate = deepcopy(state_args)
        state_args_permeate = deepcopy(state_args)

        state_args_retentate["pressure"] += initialize_guess["deltaP"]
        state_args_permeate["pressure"] = mixed_permeate_properties.pressure.value
        for j in self.config.property_package.solvent_set:
            state_args_retentate["flow_mass_phase_comp"][("Liq", j)] *= (
                1 - initialize_guess["solvent_recovery"]
            )
            state_args_permeate["flow_mass_phase_comp"][("Liq", j)] *= initialize_guess[
                "solvent_recovery"
            ]
        for j in self.config.property_package.solute_set:
            state_args_retentate["flow_mass_phase_comp"][("Liq", j)] *= (
                1 - initialize_guess["solute_recovery"]
            )
            state_args_permeate["flow_mass_phase_comp"][("Liq", j)] *= initialize_guess[
                "solute_recovery"
            ]

        state_args_interface_in = deepcopy(state_args)
        state_args_interface_out = deepcopy(state_args_retentate)

        for j in self.config.property_package.solute_set:
            state_args_interface_in["flow_mass_phase_comp"][
                ("Liq", j)
            ] *= initialize_guess["cp_modulus"]
            state_args_interface_out["flow_mass_phase_comp"][
                ("Liq", j)
            ] *= initialize_guess["cp_modulus"]

        # TODO: I think this is what we'd like to do, but IDAES needs to be changed
        state_args_interface = {}
        for t in self.flowsheet().config.time:
            for x in self.length_domain:
                assert 0.0 <= x <= 1.0
                state_args_tx = {}
                for k in state_args_interface_in:
                    if isinstance(state_args_interface_in[k], dict):
                        if k not in state_args_tx:
                            state_args_tx[k] = {}
                        for index in state_args_interface_in[k]:
                            state_args_tx[k][index] = (
                                1.0 - x
                            ) * state_args_interface_in[k][
                                index
                            ] + x * state_args_interface_out[
                                k
                            ][
                                index
                            ]
                    else:
                        state_args_tx[k] = (1.0 - x) * state_args_interface_in[
                            k
                        ] + x * state_args_interface_out[k]
                state_args_interface[t, x] = state_args_tx

        x = 0.5
        state_args_tx = {}
        for k in state_args_interface_in:
            if isinstance(state_args_interface_in[k], dict):
                if k not in state_args_tx:
                    state_args_tx[k] = {}
                for index in state_args_interface_in[k]:
                    state_args_tx[k][index] = (1.0 - x) * state_args_interface_in[k][
                        index
                    ] + x * state_args_interface_out[k][index]
            else:
                state_args_tx[k] = (1.0 - x) * state_args_interface_in[
                    k
                ] + x * state_args_interface_out[k]
        state_args_interface = state_args_tx

        return {
            "feed_side": state_args,
            "retentate": state_args_retentate,
            "permeate": state_args_permeate,
            "interface": state_args_interface,
        }


    # permeate properties need to rescale solute values by 100
    def _rescale_permeate_variable(self, var, factor=100):
        if var not in self._permeate_scaled_properties:
            sf = iscale.get_scaling_factor(var)
            iscale.set_scaling_factor(var, sf * factor)
            self._permeate_scaled_properties.add(var)

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # these variables should have user input, if not there will be a warning
        if iscale.get_scaling_factor(self.area) is None:
            sf = iscale.get_scaling_factor(self.area, default=10, warning=True)
            iscale.set_scaling_factor(self.area, sf)

        if iscale.get_scaling_factor(self.A_comp) is None:
            iscale.set_scaling_factor(self.A_comp, 1e12)

        if iscale.get_scaling_factor(self.B_comp) is None:
            iscale.set_scaling_factor(self.B_comp, 1e8)

        if iscale.get_scaling_factor(self.recovery_vol_phase) is None:
            iscale.set_scaling_factor(self.recovery_vol_phase, 1)

        for (t, p, j), v in self.recovery_mass_phase_comp.items():
            if j in self.config.property_package.solvent_set:
                sf = 1
            elif j in self.config.property_package.solute_set:
                sf = 100
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)

        for v in self.rejection_phase_comp.values():
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, 1)

        if hasattr(self, "channel_height"):
            if iscale.get_scaling_factor(self.channel_height) is None:
                iscale.set_scaling_factor(self.channel_height, 1e3)

        if hasattr(self, "spacer_porosity"):
            if iscale.get_scaling_factor(self.spacer_porosity) is None:
                iscale.set_scaling_factor(self.spacer_porosity, 1)

        if hasattr(self, "dh"):
            if iscale.get_scaling_factor(self.dh) is None:
                iscale.set_scaling_factor(self.dh, 1e3)

        if hasattr(self, "Kf"):
            for v in self.Kf.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e4)

        if hasattr(self, "N_Re"):
            for t, x in self.N_Re.keys():
                if iscale.get_scaling_factor(self.N_Re[t, x]) is None:
                    iscale.set_scaling_factor(self.N_Re[t, x], 1e-2)

        if hasattr(self, "N_Sc"):
            for t, x in self.N_Sc.keys():
                if iscale.get_scaling_factor(self.N_Sc[t, x]) is None:
                    iscale.set_scaling_factor(self.N_Sc[t, x], 1e-2)

        if hasattr(self, "N_Sh"):
            for t, x in self.N_Sh.keys():
                if iscale.get_scaling_factor(self.N_Sh[t, x]) is None:
                    iscale.set_scaling_factor(self.N_Sh[t, x], 1e-2)

        if hasattr(self, "velocity"):
            for v in self.velocity.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1)

        if hasattr(self, "friction_factor_darcy"):
            for v in self.friction_factor_darcy.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1)

        if hasattr(self, "cp_modulus"):
            for v in self.cp_modulus.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1)

        for sb in (self.permeate_side, self.mixed_permeate):
            for blk in sb.values():
                for j in self.config.property_package.solute_set:
                    self._rescale_permeate_variable(blk.flow_mass_phase_comp["Liq", j])
                    if blk.is_property_constructed("mass_frac_phase_comp"):
                        self._rescale_permeate_variable(
                            blk.mass_frac_phase_comp["Liq", j]
                        )
                    if blk.is_property_constructed("conc_mass_phase_comp"):
                        self._rescale_permeate_variable(
                            blk.conc_mass_phase_comp["Liq", j]
                        )
                    if blk.is_property_constructed("mole_frac_phase_comp"):
                        self._rescale_permeate_variable(blk.mole_frac_phase_comp[j])
                    if blk.is_property_constructed("molality_phase_comp"):
                        self._rescale_permeate_variable(
                            blk.molality_phase_comp["Liq", j]
                        )
                if blk.is_property_constructed("pressure_osm_phase"):
                    self._rescale_permeate_variable(blk.pressure_osm_phase["Liq"])

        for (t, x, p, j), v in self.flux_mass_phase_comp.items():
            if iscale.get_scaling_factor(v) is None:
                comp = self.config.property_package.get_component(j)
                if comp.is_solvent():  # scaling based on solvent flux equation
                    sf = (
                        iscale.get_scaling_factor(self.A_comp[t, j])
                        * iscale.get_scaling_factor(self.dens_solvent)
                        * iscale.get_scaling_factor(
                            self.feed_side.properties[t, x].pressure
                        )
                    )
                    iscale.set_scaling_factor(v, sf)
                elif comp.is_solute():  # scaling based on solute flux equation
                    sf = iscale.get_scaling_factor(
                        self.B_comp[t, j]
                    ) * iscale.get_scaling_factor(
                        self.feed_side.properties[t, x].conc_mass_phase_comp[p, j]
                    )
                    iscale.set_scaling_factor(v, sf)
