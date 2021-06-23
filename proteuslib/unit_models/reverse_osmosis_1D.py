###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################


from enum import Enum
# Import Pyomo libraries
from pyomo.environ import (Var,
                           Set,
                           Param,
                           SolverFactory,
                           Suffix,
                           NonNegativeReals,
                           NegativeReals,
                           Reference,
                           Block,
                           units as pyunits,
                           exp,
                           value)
from pyomo.common.config import ConfigBlock, ConfigValue, In
# Import IDAES cores
from idaes.core import (ControlVolume1DBlock,
                        declare_process_block_class,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        UnitModelBlockData,
                        useDefault,
                        FlowDirection)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util import get_solver, scaling as iscale

import idaes.logger as idaeslog


__author__ = "Adam Atia"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("ReverseOsmosis1D")
class ReverseOsmosis1DData(UnitModelBlockData):
    """Standard 1D Reverse Osmosis Unit Model Class."""

    CONFIG = UnitModelBlockData.CONFIG()

    # Template for config arguments for feed and permeate side
    _SideTemplate = ConfigBlock()

    _SideTemplate.declare("dynamic", ConfigValue(
        default=False,
        domain=In([False]),
        description="Dynamic model flag - must be False",
        doc="""Indicates whether this model will be dynamic or not.
    **default** = False. RO units do not yet support dynamic
    behavior."""))

    _SideTemplate.declare("has_holdup", ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False. RO units do not have defined volume, thus
    this must be False."""))

    _SideTemplate.declare("material_balance_type", ConfigValue(
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
    **MaterialBalanceType.total** - use total material balance.}"""))

    _SideTemplate.declare("energy_balance_type", ConfigValue(
            default=EnergyBalanceType.useDefault,
            domain=In(EnergyBalanceType),
            description="Energy balance construction flag",
            doc="""Indicates what type of energy balance should be constructed,
    **default** - EnergyBalanceType.useDefault.
    **Valid values:** {
    **EnergyBalanceType.useDefault - refer to property package for default
    balance type
    **EnergyBalanceType.none** - exclude energy balances,
    **EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
    **EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
    **EnergyBalanceType.energyTotal** - single energy balance for material,
    **EnergyBalanceType.energyPhase** - energy balances for each phase.}"""))

    _SideTemplate.declare("momentum_balance_type", ConfigValue(
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
    **MomentumBalanceType.momentumPhase** - momentum balances for each phase.}"""))

    _SideTemplate.declare("has_pressure_change", ConfigValue(
            default=False,
            domain=In([True, False]),
            description="Pressure change term construction flag",
            doc="""Indicates whether terms for pressure change should be
    constructed,
    **default** - False.
    **Valid values:** {
    **True** - include pressure change terms,
    **False** - exclude pressure change terms.}"""))

    _SideTemplate.declare("property_package", ConfigValue(
            default=None,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations
    **default** - useDefault.
    **Valid values:** {
    **useDefault** - use default package from parent model or flowsheet,
    **PhysicalParameterObject** - a PhysicalParameterBlock object.}"""))

    _SideTemplate.declare("property_package_args", ConfigValue(
            default={},
            description="Arguments for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
    and used when constructing these.
    **default** - None.
    **Valid values:** {
    see property package for documentation.}"""))

    _SideTemplate.declare(
        "transformation_method",
        ConfigValue(
            default=useDefault,
            description="Discretization method to use for DAE transformation",
            doc="""Discretization method to use for DAE transformation. See Pyomo
    documentation for supported transformations."""))

    _SideTemplate.declare("transformation_scheme", ConfigValue(
            default=useDefault,
            description="Discretization scheme to use for DAE transformation",
            doc="""Discretization scheme to use when transforming domain. See
    Pyomo documentation for supported schemes."""))

    # Create individual config blocks for feed and permeate side
    CONFIG.declare("feed_side", _SideTemplate(doc="feed side config arguments"))
    CONFIG.declare("permeate_side", _SideTemplate(doc="permeate side config arguments"))

    # Common config args for both sides
    CONFIG.declare("finite_elements", ConfigValue(
            default=20,
            domain=int,
            description="Number of finite elements in length domain",
            doc="""Number of finite elements to use when discretizing length 
            domain (default=20)"""))

    CONFIG.declare("collocation_points", ConfigValue(
            default=5,
            domain=int,
            description="Number of collocation points per finite element",
            doc="""Number of collocation points to use per finite element when
            discretizing length domain (default=5)"""))

    def _process_config(self):
        pass  #TODO: add config errors here

    def build(self):
        """
        Build 1D RO model (pre-DAE transformation).

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        self._process_config()

        units_meta = self.config.property_package.get_metadata().get_derived_units

        # ==========================================================================
        """ Build 1D Control volume for feed side"""
        self.feed_side = ControlVolume1DBlock(default={
            "dynamic": self.config.feed_side.dynamic,
            "has_holdup": self.config.feed_side.has_holdup,
            "property_package": self.config.feed_side.property_package,
            "property_package_args": self.config.feed_side.property_package_args,
            "transformation_method": self.config.transformation_method,
            "transformation_scheme": self.config.transformation_scheme,
            "finite_elements": self.config.finite_elements,
            "collocation_points" self.config.collocation_points
        })


        # ==========================================================================
        """ Build 1D Control volume for permeate side"""
        self.permeate_side = ControlVolume1DBlock(default={
            "dynamic": self.config.permeate_side.dynamic,
            "has_holdup": self.config.permeate_side.has_holdup,
            "property_package": self.config.permeate_side.property_package,
            "property_package_args": self.config.permeate_side.property_package_args,
            "transformation_method": self.config.transformation_method,
            "transformation_scheme": self.config.transformation_scheme,
            "finite_elements": self.config.finite_elements,
            "collocation_points": self.config.collocation_points
        })

        feed_side = self.feed_side
        permeate_side = self.permeate_side

        # ==========================================================================
        """ Add geometry for each control volume"""
        feed_side.add_geometry()
        permeate_side.add_geometry()

        # ==========================================================================
        """ Add state blocks for each control volume"""
        feed_side.add_state_blocks(has_phase_equilibrium=False)
        permeate_side.add_state_blocks(has_phase_equilibrium=False)

        # ==========================================================================
        """ Populate feed side"""
        feed_side.add_material_balances(balance_type=self.config.material_balance_type,
                                        has_mass_transfer=True)
        feed_side.add_energy_balances(balance_type=self.config.energy_balance_type,
                                           has_enthalpy_transfer=True)
        feed_side.add_momentum_balances(balance_type=self.config.momentum_balance_type,
                                        has_pressure_change=self.config.has_pressure_change)

        # ==========================================================================
        """ Only enable mass transfer for permeate side"""
        permeate_side.add_material_balances(balance_type=self.config.material_balance_type.none,
                                            has_mass_transfer=True)

        # ==========================================================================
        """ Apply transformation to feed and permeate sides"""
        feed_side.apply_transformation()
        permeate_side.apply_transformation()

        # ==========================================================================
        """ Add inlet/outlet ports for feed side and only an outlet port for permeate side"""
        feed_side.add_inlet_port(name="feed_inlet", block=feed_side)
        feed_side.add_outlet_port(name="feed_outlet", block=feed_side)
        permeate_side.add_outlet_port(name="permeate_outlet", block=permeate_side)