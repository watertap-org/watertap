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

# Import Pyomo libraries
from pyomo.environ import (
    Block,
    Set,
    Var,
    Suffix,
    NonNegativeReals,
    Reals,
    Reference,
    units as pyunits,
    log,
    value,
    Expr_if,
    exp,
    check_optimal_termination,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.common.collections import ComponentSet

# Import IDAES cores
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    MaterialBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
    MaterialFlowBasis,
)
from idaes.core.solvers import get_solver
from idaes.core.util.math import smooth_min
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError, InitializationError
import idaes.core.util.scaling as iscale
from idaes.core.util.constants import Constants
from watertap.core.util.initialization import check_solve, check_dof

import idaes.logger as idaeslog


_log = idaeslog.getLogger(__name__)

# TODO:
# - Further refinement of initialization routine to make more robust
# - Further refinement of scaling to make more robust
# - Add more tests to test robustness
# - Add constraints for computing pressure drop in spiral wound membrane


class MassTransferCoefficient(Enum):
    none = auto()
    # mass transfer coefficient is a user specified value
    fixed = auto()
    # mass transfer coefficient is calculated using spiral wound correlation
    spiral_wound = auto()


class ConcentrationPolarizationType(Enum):
    none = auto()
    calculated = auto()


@declare_process_block_class("NanofiltrationDSPMDE0D")
class NanofiltrationData(UnitModelBlockData):
    """
    Nanofiltration model based on Donnan Steric Pore Model with Dielectric Exclusion (DSPM-DE).

    Assumptions
        - Membrane electric potential at membrane interface is taken as reference (i.e., equal to 0)

    References:
        Geraldes and Alves, 2008 (https://doi.org/10.1016/j.memsci.2008.04.054)
        Roy et al., 2015 (http://dx.doi.org/10.1016/j.memsci.2015.06.030)
        Labban et al., 2017 (http://dx.doi.org/10.1016/j.memsci.2016.08.062)
    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False. NF units do not support dynamic
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
    **default** - False. NF units do not have defined volume, thus
    this must be False.""",
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
    and used when constructing these,
    **default** - None.
    **Valid values:** {
    see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "mass_transfer_coefficient",
        ConfigValue(
            default=MassTransferCoefficient.spiral_wound,
            domain=In(MassTransferCoefficient),
            description="Mass transfer coefficient in feed channel",
            doc="""
            Options to account for mass transfer coefficient.

            **default** - ``MassTransferCoefficient.fixed``

        .. csv-table::
            :header: "Configuration Options", "Description"

            "``MassTransferCoefficient.none``", "Simplifying assumption to ignore mass transfer coefficient"
            "``MassTransferCoefficient.fixed``", "Specify an estimated value for the mass transfer coefficient in the feed channel"
            "``MassTransferCoefficient.spiral_wound``", "Allow model to perform calculation of mass transfer coefficient based on
            spiral wound module correlation"
        """,
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
            "``ConcentrationPolarizationType.calculated``", "Allow model to perform calculation of membrane-interface concentration"
        """,
        ),
    )

    def _process_config(self):
        if len(self.config.property_package.solvent_set) > 1:
            raise ConfigurationError(
                "NF model only supports one solvent component,"
                "the provided property package has specified {} solvent components".format(
                    len(self.config.property_package.solvent_set)
                )
            )

        if len(self.config.property_package.solvent_set) == 0:
            raise ConfigurationError(
                "The NF model was expecting a solvent and did not receive it."
            )

        if (
            hasattr(self.config.property_package, "ion_set")
            and len(self.config.property_package.ion_set) == 0
        ) and (
            hasattr(self.config.property_package, "solute_set")
            and len(self.config.property_package.solute_set) == 0
        ):
            raise ConfigurationError(
                "This NF model was expecting ions/solutes and did not receive any."
            )
        if (
            self.config.concentration_polarization_type
            != ConcentrationPolarizationType.calculated
            and self.config.mass_transfer_coefficient != MassTransferCoefficient.none
        ) or (
            self.config.concentration_polarization_type
            == ConcentrationPolarizationType.calculated
            and self.config.mass_transfer_coefficient == MassTransferCoefficient.none
        ):
            raise ConfigurationError(
                "\nConflict between configuration options:\n"
                "'mass_transfer_coefficient' cannot be set to {} "
                "while 'concentration_polarization_type' is set to {}.\n\n"
                "'mass_transfer_coefficient' must be set to MassTransferCoefficient.none\nor "
                "'concentration_polarization_type' must be set to ConcentrationPolarizationType.calculated".format(
                    self.config.mass_transfer_coefficient,
                    self.config.concentration_polarization_type,
                )
            )

    def build(self):
        super().build()

        self._process_config()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # For stateblock-specific scaling in calculate_scaling_factors
        self._sb_scaled_properties = ComponentSet()

        units_meta = self.config.property_package.get_metadata().get_derived_units

        self.io_list = io_list = Set(initialize=[0, 1])  # inlet/outlet set

        # These two sets should always be present
        #   and their union is the full set of dissolved species
        if hasattr(self.config.property_package, "ion_set") and hasattr(
            self.config.property_package, "solute_set"
        ):
            self.solute_set = solute_set = (
                self.config.property_package.ion_set
                | self.config.property_package.solute_set
            )
        else:
            raise ConfigurationError(
                "This NF model was expecting an "
                "ion_set and solute_set and did not them."
            )

        solvent_set = self.config.property_package.solvent_set
        solvent_solute_set = solvent_set | solute_set
        phase_list = self.config.property_package.phase_list

        # Build control volume for feed side
        self.feed_side = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        self.feed_side.add_state_blocks(has_phase_equilibrium=False)

        self.feed_side.add_material_balances(
            balance_type=self.config.material_balance_type, has_mass_transfer=True
        )

        self.feed_side.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change,
        )

        # Make indexed stateblock and separate stateblock for permeate-side and permeate outlet, respectively.
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False  # these blocks are not inlets

        # Add permeate block
        self.permeate_side = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            io_list,
            doc="Material properties of permeate along permeate channel",
            **tmp_dict,
        )
        self.mixed_permeate = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of mixed permeate exiting the module",
            **tmp_dict,
        )

        # Add Ports
        self.add_inlet_port(name="inlet", block=self.feed_side)
        self.add_outlet_port(name="retentate", block=self.feed_side)
        self.add_port(name="permeate", block=self.mixed_permeate)

        # Membrane interface: indexed state block
        self.feed_side.properties_interface = (
            self.config.property_package.state_block_class(
                self.flowsheet().config.time,
                io_list,
                doc="Material properties of feed-side membrane interface",
                **tmp_dict,
            )
        )
        # Pore entrance: indexed state block
        self.pore_entrance = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            io_list,
            doc="Fluid properties within the membrane pore entrance",
            **tmp_dict,
        )
        # Pore exit: indexed state block
        self.pore_exit = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            io_list,
            doc="Fluid properties within the membrane pore exit",
            **tmp_dict,
        )

        # References for control volume
        # pressure change
        if (
            self.config.has_pressure_change is True
            and self.config.momentum_balance_type != MomentumBalanceType.none
        ):
            self.deltaP = Reference(self.feed_side.deltaP)

        ###############################################################################################################
        # Variables
        ###############################################################################################################
        # 1. Component mole flux, J, DOF=Nj*2 for inlet/outlet
        self.flux_mol_phase_comp = Var(
            self.flowsheet().config.time,
            io_list,
            phase_list,
            solvent_solute_set,
            initialize=(lambda b, t, x, p, j: 2.78e-2 if j in solvent_set else 4e-6),
            # bounds=lambda b, t, x, p, j: (0, 1.5)
            # if j in solvent_set
            # else (0, 1e-2),  # TODO: keep checking these
            domain=NonNegativeReals,
            units=units_meta("amount")
            * units_meta("length") ** -2
            * units_meta("time") ** -1,
            doc="Component molar flux at inlet and outlet of membrane",
        )

        # 2. Pore radius, rp, DOF = 1
        self.radius_pore = Var(
            initialize=0.5e-9,  # TODO: revisit
            domain=NonNegativeReals,
            units=units_meta("length"),
            doc="Membrane pore radius",
        )

        # 3. Effective thickness of membrane, Ak, DOF=1
        self.membrane_thickness_effective = Var(
            initialize=1.33e-6,  # Value used by Labban et al., 2017
            domain=NonNegativeReals,
            units=units_meta("length"),
            doc="Effective membrane thickness",
        )

        # 4. Effective thickness of membrane, Ak, DOF=1
        self.membrane_charge_density = Var(
            self.flowsheet().config.time,
            initialize=-50,  # near value used in Roy et al.
            domain=Reals,
            units=pyunits.mol * pyunits.m**-3,
            doc="Membrane charge density",
        )
        self.dielectric_constant_pore = Var(
            self.flowsheet().config.time,
            initialize=42,  # near value used in Roy et al.
            bounds=(1, None),
            units=pyunits.dimensionless,  # TODO: revisit bounds/domain
            doc="Pore dielectric constant",
        )
        self.electric_potential = Var(
            self.flowsheet().config.time,
            io_list,
            ["pore_entrance", "pore_exit", "permeate"],
            initialize=-1e-3,  # TODO:revisit
            domain=Reals,
            units=pyunits.V,
            doc="Electric potential of pore entrance/exit, and permeate",
        )
        self.rejection_intrinsic_phase_comp = Var(
            self.flowsheet().config.time,
            phase_list,
            solute_set,
            initialize=0.1,
            bounds=(-1.001, 1.001),
            units=pyunits.dimensionless,
            doc="Intrinsic solute rejection",
        )
        self.area = Var(
            initialize=50,
            bounds=(0, 1e3),
            domain=NonNegativeReals,
            units=units_meta("length") ** 2,
            doc="Membrane area",
        )

        self.recovery_vol_phase = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            initialize=0.5,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Volumetric-based recovery",
        )
        if (
            self.config.concentration_polarization_type
            == ConcentrationPolarizationType.calculated
        ):
            self.Kf_comp = Var(
                self.flowsheet().config.time,
                self.io_list,
                solute_set,
                initialize=5e-5,
                bounds=(0, 1e-3),
                domain=NonNegativeReals,
                units=units_meta("length") * units_meta("time") ** -1,
                doc="Component mass transfer coefficient in feed channel at inlet and outlet",
            )
            self.electric_potential_grad_feed_interface = Var(
                self.flowsheet().config.time,
                io_list,
                initialize=-1e-8,  # TODO: revisit
                domain=Reals,
                # TODO: revisit- Geraldes and Alves give unitless while Roy et al. give V/m
                units=pyunits.V * pyunits.m**-1,
                doc="Electric potential gradient of feed-membrane interface",
            )
        else:
            pass

        if (
            self.config.mass_transfer_coefficient
            == MassTransferCoefficient.spiral_wound
        ):

            self.N_Sc_comp = Var(
                self.flowsheet().config.time,
                self.io_list,
                solute_set,
                initialize=5e2,
                bounds=(1e2, 2e3),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Schmidt number at inlet and outlet",
            )
            self.N_Pe_comp = Var(
                self.flowsheet().config.time,
                self.io_list,
                solute_set,
                initialize=1e5,
                # bounds=(5e3, None),  # TODO:unsure of value ranges at the moment
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Peclet number at inlet and outlet",
            )
            self.spacer_mixing_efficiency = Var(
                initialize=0.5,
                # bounds=(1e2, 2e3),  # TODO:unsure of value ranges at the moment- since this is efficiency, assuming 0 -1
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Mixing efficiency of spacer net",
            )
            self.spacer_mixing_length = Var(
                initialize=0.6,
                # bounds=(1e2, 5e3),  #TODO:unsure of value ranges at the moment
                domain=NonNegativeReals,
                units=units_meta("length"),
                doc="Characteristic length of spacer",
            )
        else:
            pass

        self.length = Var(
            initialize=10,
            bounds=(0, 5e2),
            domain=NonNegativeReals,
            units=units_meta("length"),
            doc="Effective membrane length",
        )
        self.width = Var(
            initialize=5,
            bounds=(0, 5e2),
            domain=NonNegativeReals,
            units=units_meta("length"),
            doc="Effective feed-channel width",
        )
        self.channel_height = Var(
            initialize=1e-3,
            domain=NonNegativeReals,
            bounds=(0, 5e-3),
            units=units_meta("length"),
            doc="Feed channel height",
        )
        self.spacer_porosity = Var(
            initialize=0.95,
            bounds=(0.1, 1.001),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Feed-channel spacer porosity",
        )
        self.velocity = Var(
            self.flowsheet().config.time,
            self.io_list,
            initialize=0.5,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("length") / units_meta("time"),
            doc="Crossflow velocity in feed channel at inlet and outlet",
        )
        ###############################################################################################################
        # Expressions
        ###############################################################################################################
        # Make expressions that don't depend on any variables
        self._make_expressions()

        ################################################################################################################
        # Constraints
        ################################################################################################################

        # 0. Membrane area
        @self.Constraint(doc="Membrane area")
        def eq_area(b):
            return b.area == b.length * b.width

        # 1. Feed-solution/membrane equilibrium, DOF= Nj * 2 for inlet/outlet
        @self.Constraint(
            self.flowsheet().config.time,
            io_list,
            phase_list,
            solute_set,
            doc="Interfacial partitioning at feed side of membrane",
        )
        def eq_interfacial_partitioning_feed(b, t, x, p, j):
            return (
                b.pore_entrance[t, x].act_coeff_phase_comp[p, j]
                * b.pore_entrance[t, x].conc_mol_phase_comp[p, j]
                / (
                    b.feed_side.properties_interface[t, x].act_coeff_phase_comp[p, j]
                    * b.feed_side.properties_interface[t, x].conc_mol_phase_comp[p, j]
                )
                == b.partition_factor_steric_comp[t, j]
                * b.partition_factor_born_solvation_comp[t, j]
                * b.partition_factor_donnan_comp_feed[t, x, j]
            )

        # 2. Permeate solution/membrane equilibrium, DOF= Nj * 2 for inlet/outlet
        @self.Constraint(
            self.flowsheet().config.time,
            io_list,
            phase_list,
            solute_set,
            doc="Interfacial partitioning at permeate side of membrane",
        )
        def eq_interfacial_partitioning_permeate(b, t, x, p, j):
            return (
                b.pore_exit[t, x].act_coeff_phase_comp[p, j]
                * b.pore_exit[t, x].conc_mol_phase_comp[p, j]
                / (
                    b.permeate_side[t, x].act_coeff_phase_comp[p, j]
                    * b.permeate_side[t, x].conc_mol_phase_comp[p, j]
                )
                == b.partition_factor_steric_comp[t, j]
                * b.partition_factor_born_solvation_comp[t, j]
                * b.partition_factor_donnan_comp_permeate[t, x, j]
            )

        # 4. Charge balance inside the membrane, DOF=N nodes across membrane thickness *2 for inlet/outlet: N=2, DOF=4
        @self.Constraint(
            self.flowsheet().config.time,
            io_list,
            ["pore_entrance", "pore_exit"],
            phase_list,
            doc="Electroneutrality within membrane pore",
        )
        def eq_electroneutrality_pore(b, t, x, y, p):
            if y == "pore_entrance":
                pore_loc = b.pore_entrance[t, x]
            elif y == "pore_exit":
                pore_loc = b.pore_exit[t, x]
            return (
                sum(
                    pore_loc.conc_mol_phase_comp[p, j] * pore_loc.charge_comp[j]
                    for j in solute_set
                )
                + b.membrane_charge_density[t]
                == 0
            )

        # 4. Permeate electroneutrality, DOF=1 *2 for inlet/outlet:  DOF=2
        @self.Constraint(
            self.flowsheet().config.time,
            io_list,
            phase_list,
            doc="Electroneutrality in permeate",
        )
        def eq_electroneutrality_permeate(b, t, x, p):
            return (
                sum(
                    b.permeate_side[t, x].conc_mol_phase_comp[p, j]
                    * b.permeate_side[t, x].charge_comp[j]
                    for j in solute_set
                )
                == 0
            )

        # 5. Water flux via Hagen-Poiseuille relationship, DOF= 1 * 2 for inlet/outlet: DOF= 2
        @self.Constraint(
            self.flowsheet().config.time,
            io_list,
            phase_list,
            doc="Hagen-Poiseuille relationship for water flux across membrane",
        )
        def eq_water_flux(b, t, x, p):
            if x == 0:
                prop_feed = b.feed_side.properties_in[t]
            elif x == 1:
                prop_feed = b.feed_side.properties_out[t]
            prop_perm = b.permeate_side[t, x]
            prop_feed_inter = b.feed_side.properties_interface[t, x]
            return b.flux_vol_water[t, x] == (
                prop_feed.pressure
                - prop_perm.pressure
                - (
                    prop_feed_inter.pressure_osm_phase["Liq"]
                    - prop_perm.pressure_osm_phase["Liq"]
                )
            ) * (b.radius_pore**2) / (
                8 * prop_feed.visc_d_phase[p] * b.membrane_thickness_effective
            )

        # 6. Unhindered mass transfer; Js,i=Jw*cp,i; DOF= Nj * 2 for inlet/outlet
        @self.Constraint(
            self.flowsheet().config.time,
            io_list,
            phase_list,
            solute_set,
            doc="Solute flux as function of solvent flux",
        )
        def eq_solute_solvent_flux(b, t, x, p, j):
            return (
                b.flux_mol_phase_comp[t, x, p, j]
                == b.flux_vol_water[t, x]
                * b.permeate_side[t, x].conc_mol_phase_comp[p, j]
            )

        # TESTING PROBLEMATIC CONSTRAINT RESULTING IN ERRONEOUSLY LOW REJECTION:
        @self.Expression(
            self.flowsheet().config.time,
            io_list,
            phase_list,
            solute_set,
            doc="Diffusive transport across membrane pore",
        )
        def diffusive_term(b, t, x, p, j):
            return (
                -b.diffus_pore_comp[t, j]
                * (
                    b.pore_exit[t, x].conc_mol_phase_comp[p, j]
                    - b.pore_entrance[t, x].conc_mol_phase_comp[p, j]
                )
                / b.membrane_thickness_effective
            )

        @self.Expression(
            self.flowsheet().config.time,
            io_list,
            phase_list,
            solute_set,
            doc="Convective transport across membrane pore",
        )
        def convective_term(b, t, x, p, j):
            return (
                b.hindrance_factor_convective_comp[t, j]
                * b.conc_mol_phase_comp_pore_avg[t, x, p, j]
                * b.flux_vol_water[t, x]
            )

        @self.Expression(
            self.flowsheet().config.time,
            io_list,
            phase_list,
            solute_set,
            doc="Electromigrative transport across membrane pore",
        )
        def electromigration_term(b, t, x, p, j):
            return (
                -b.feed_side.properties_in[t].charge_comp[j]
                * b.conc_mol_phase_comp_pore_avg[t, x, p, j]
                * b.diffus_pore_comp[t, j]
                * Constants.faraday_constant
                / (Constants.gas_constant * b.feed_side.properties_in[t].temperature)
                * (
                    b.electric_potential[t, x, "pore_exit"]
                    - b.electric_potential[t, x, "pore_entrance"]
                )
                / b.membrane_thickness_effective
            )

        # 7. Extended Nernst Planck equation, DOF= Nj * 2 for inlet/outlet
        @self.Constraint(
            self.flowsheet().config.time,
            io_list,
            phase_list,
            solute_set,
            doc="Solute flux within pore domain",
        )
        def eq_solute_flux_pore_domain(b, t, x, p, j):
            return (
                b.flux_mol_phase_comp[t, x, p, j]
                == b.diffusive_term[t, x, p, j]
                + b.convective_term[t, x, p, j]
                + b.electromigration_term[t, x, p, j]
            )

        # 8. Feed-solution/membrane mass transfer resistance, DOF= Nj * 2 for inlet/outlet
        @self.Constraint(
            self.flowsheet().config.time,
            io_list,
            phase_list,
            solute_set,
            doc="Feed-interface mass transfer resistance accounting for concentration polarization",
        )
        def eq_solute_flux_concentration_polarization(b, t, x, p, j):
            if x == 0:
                bulk = b.feed_side.properties_in[t]
            elif x:
                bulk = b.feed_side.properties_out[t]
            interface = b.feed_side.properties_interface[t, x]
            if (
                self.config.concentration_polarization_type
                == ConcentrationPolarizationType.calculated
            ):
                return (
                    b.flux_mol_phase_comp[t, x, p, j]
                    == -b.Kf_comp[t, x, j]
                    * (
                        interface.conc_mol_phase_comp[p, j]
                        - bulk.conc_mol_phase_comp[p, j]
                    )
                    + b.flux_vol_water[t, x] * interface.conc_mol_phase_comp[p, j]
                    - interface.charge_comp[j]
                    * interface.conc_mol_phase_comp[p, j]
                    * interface.diffus_phase_comp[p, j]
                    * Constants.faraday_constant
                    / Constants.gas_constant
                    / interface.temperature
                    * b.electric_potential_grad_feed_interface[t, x]
                )
            elif (
                self.config.concentration_polarization_type
                == ConcentrationPolarizationType.none
            ):
                return (
                    interface.conc_mol_phase_comp[p, j]
                    == bulk.conc_mol_phase_comp[p, j]
                )
            else:
                raise ConfigurationError(
                    "Provide a valid value for concentration_polarization_type"
                )

        # 9. Isothermal conditions at permeate inlet/outlet, DOF= 1*2 for inlet/outlet
        @self.Constraint(
            self.flowsheet().config.time,
            io_list,
            doc="Isothermal assumption for permeate",
        )
        def eq_permeate_isothermal(b, t, x):
            return (
                b.feed_side.properties_in[t].temperature
                == b.permeate_side[t, x].temperature
            )

        # 10. Isothermal conditions at feed/membrane interface, DOF= 1*2 for inlet/outlet
        @self.feed_side.Constraint(
            self.flowsheet().config.time,
            io_list,
            doc="Isothermal assumption for feed-membrane interface",
        )
        def eq_feed_interface_isothermal(b, t, x):
            return (
                b.properties_in[t].temperature
                == b.properties_interface[t, x].temperature
            )

        # 11. Isobaric conditions at permeate side, DOF= 1*2 for inlet/outlet
        # TOdo: mixed permeate pressure is currently unused variable, but fixing its value satisfies this constraint
        @self.Constraint(
            self.flowsheet().config.time, io_list, doc="Isobaric permeate pressure"
        )
        def eq_pressure_permeate_io(b, t, x):
            return b.permeate_side[t, x].pressure == b.mixed_permeate[t].pressure

        # 12. Javg * area = -dMf
        @self.Constraint(
            self.flowsheet().config.time,
            phase_list,
            solvent_solute_set,
            doc="Component mass transfer from feed",
        )
        def eq_mass_transfer_feed(b, t, p, j):
            if (
                b.feed_side.properties_in[0].get_material_flow_basis()
                == MaterialFlowBasis.molar
            ):
                return (
                    b.flux_mol_phase_comp_avg[t, p, j] * b.area
                    == -b.feed_side.mass_transfer_term[t, p, j]
                )

        # 13. Mass transfer equal to permeate flow terms; mole_flow,perm final = -dMf = Javg * area
        @self.Constraint(
            self.flowsheet().config.time,
            phase_list,
            solvent_solute_set,
            doc="Permeate production/average mass transfer constraint",
        )
        def eq_permeate_production(b, t, p, j):
            if b.mixed_permeate[0].get_material_flow_basis() == MaterialFlowBasis.molar:
                return (
                    b.mixed_permeate[t].get_material_flow_terms(p, j)
                    == b.flux_mol_phase_comp_avg[t, p, j] * b.area
                )

        # 14. Volumetric recovery rate
        @self.Constraint(self.flowsheet().config.time, phase_list)
        def eq_recovery_vol_phase(b, t, p):
            return (
                b.recovery_vol_phase[t, p]
                * b.feed_side.properties_in[t].flow_vol_phase[p]
                == b.mixed_permeate[t].flow_vol_phase[p]
            )

        # 15. Intrinsic rejection rate
        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            solute_set,
            doc="Intrinsic solute rejection",
        )
        def eq_rejection_intrinsic_phase_comp(b, t, p, j):
            return (
                b.rejection_intrinsic_phase_comp[t, p, j]
                == 1
                - b.mixed_permeate[t].conc_mol_phase_comp[p, j]
                / b.feed_side.properties_interface[t, 0].conc_mol_phase_comp[p, j]
            )

        if (
            self.config.concentration_polarization_type
            == ConcentrationPolarizationType.calculated
        ):
            # 3. Feed-solution/membrane electroneutrality, DOF=1 *2 for inlet/outlet: DOF= 2
            @self.Constraint(
                self.flowsheet().config.time,
                io_list,
                phase_list,
                doc="Electroneutrality at feed-side membrane interface",
            )
            def eq_electroneutrality_interface(b, t, x, p):
                return (
                    sum(
                        b.feed_side.properties_interface[t, x].conc_mol_phase_comp[p, j]
                        * b.feed_side.properties_interface[t, x].charge_comp[j]
                        for j in solute_set
                    )
                    == 0
                )

            if (
                self.config.mass_transfer_coefficient
                == MassTransferCoefficient.spiral_wound
            ):
                # 16. Mass transfer coefficient
                @self.Constraint(
                    self.flowsheet().config.time,
                    io_list,
                    solute_set,
                    doc="Mass transfer coefficient",
                )
                def eq_Kf_comp(b, t, x, j):
                    bulk_diff = b.feed_side.properties_in[t].diffus_phase_comp["Liq", j]
                    return (
                        b.Kf_comp[t, x, j]
                        == 0.753
                        * (
                            b.spacer_mixing_efficiency
                            / (2 - b.spacer_mixing_efficiency)
                        )
                        ** 0.5
                        * (2 * bulk_diff / b.channel_height)
                        * b.N_Sc_comp[t, x, j] ** (-1 / 6)
                        * (
                            b.N_Pe_comp[t, x, j]
                            * b.channel_height
                            / (2 * b.spacer_mixing_length)
                        )
                        ** 0.5
                    )
                    # TODO: NOTE--- error in MIT paper; 1/2 of channel height should be in numerator

                # 17. Schmidt number calculation
                @self.Constraint(
                    self.flowsheet().config.time,
                    io_list,
                    solute_set,
                    doc="Schmidt number equation",
                )
                def eq_N_Sc_comp(b, t, x, j):
                    if not x:
                        prop_io = b.feed_side.properties_in[t]
                    elif x:
                        prop_io = b.feed_side.properties_out[t]
                    return (
                        b.N_Sc_comp[t, x, j]
                        * prop_io.dens_mass_phase["Liq"]
                        * prop_io.diffus_phase_comp["Liq", j]
                        == prop_io.visc_d_phase["Liq"]
                    )

                # 18. Peclet number calculation
                @self.Constraint(
                    self.flowsheet().config.time,
                    io_list,
                    solute_set,
                    doc="Peclet number equation",
                )
                def eq_N_Pe_comp(b, t, x, j):
                    bulk_diff = b.feed_side.properties_in[t].diffus_phase_comp["Liq", j]
                    return (
                        b.N_Pe_comp[t, x, j]
                        == 2 * b.channel_height * b.velocity[t, x] / bulk_diff
                    )

        # 19. Crossflow velocity at inlet and outlet
        @self.Constraint(
            self.flowsheet().config.time,
            self.io_list,
            doc="Crossflow velocity constraint",
        )
        def eq_velocity(b, t, x):
            if not x:
                prop_io = b.feed_side.properties_in[t]
            elif x:
                prop_io = b.feed_side.properties_out[t]
            return b.velocity[t, x] * b.area_cross == prop_io.flow_vol_phase["Liq"]

        # TODO: seems stale since temperature unused at pore entrance/exit- confirm+remove;
        #  1/17/22: after including temp variables for pore in interfacial equilib eqns, this is relevant
        @self.Constraint(
            self.flowsheet().config.time,
            io_list,
            ["pore_entrance", "pore_exit"],
            doc="Isothermal assumption for pore inlet/outlet",
        )
        def eq_pore_isothermal(b, t, x, y):
            if y == "pore_entrance":
                prop = b.pore_entrance[t, x]
            elif y == "pore_exit":
                prop = b.pore_exit[t, x]
            return b.feed_side.properties_in[t].temperature == prop.temperature

        # Experimental Constraint with new density calculation in prop package-- temp equality in permeate
        @self.Constraint(
            self.flowsheet().config.time,
            doc="Isothermal assumption for mixed permeate",
        )
        def eq_permeate_isothermal_mixed(b, t):
            return (
                b.feed_side.properties_in[t].temperature
                == b.mixed_permeate[t].temperature
            )

        # Experimental constraint: noticed feed outlet temp didn't match inlet
        @self.feed_side.Constraint(
            self.flowsheet().config.time, doc="Isothermal assumption for feed-outlet"
        )
        def eq_feed_isothermal(b, t):
            return b.properties_in[t].temperature == b.properties_out[t].temperature

        # Experimental constraint
        @self.feed_side.Constraint(
            self.flowsheet().config.time,
            io_list,
            doc="Equal volumetric flow for bulk feed and interface",
        )
        def eq_equal_flow_vol_feed_interface(b, t, x):
            if not x:
                bulk = b.properties_in[t]
            else:
                bulk = b.properties_out[t]
            return (
                bulk.flow_vol_phase["Liq"]
                == b.properties_interface[t, x].flow_vol_phase["Liq"]
            )

        # Experimental constraint
        @self.Constraint(
            self.flowsheet().config.time,
            io_list,
            doc="Equal volumetric flow for pore entrance/exit at inlet and outlet",
        )
        def eq_equal_flow_vol_pore(b, t, x):
            return (
                b.pore_entrance[t, x].flow_vol_phase["Liq"]
                == b.pore_exit[t, x].flow_vol_phase["Liq"]
            )

        # Experimental constraint
        @self.Constraint(
            self.flowsheet().config.time,
            io_list,
            doc="Volumetric flow at pore exit inlet and outlet equal to mixed permeate",
        )
        def eq_equal_flow_vol_pore_exit_perm(b, t, x):
            return (
                b.pore_exit[t, x].flow_vol_phase["Liq"]
                == b.mixed_permeate[t].flow_vol_phase["Liq"]
            )

        # Experimental constraint
        @self.Constraint(
            self.flowsheet().config.time,
            io_list,
            doc="Volumetric flow at permeate of inlet and outlet equal to mixed permeate",
        )
        def eq_equal_flow_vol_permeate(b, t, x):
            return (
                b.permeate_side[t, x].flow_vol_phase["Liq"]
                == b.mixed_permeate[t].flow_vol_phase["Liq"]
            )

    def _make_expressions(self):
        solute_set = self.solute_set
        io_list = self.io_list

        # Stokes radius to membrane pore radius ratio (for each solute)
        @self.Expression(
            self.flowsheet().config.time,
            solute_set,
            doc="Ratio of stokes radius to membrane pore radius equation",
        )
        def lambda_comp(b, t, j):
            return smooth_min(
                1, b.feed_side.properties_in[t].radius_stokes_comp[j] / b.radius_pore
            )

        @self.Expression(
            self.flowsheet().config.time,
            solute_set,
            doc="Diffusive hindered transport coefficient",
        )
        def hindrance_factor_diffusive_comp(b, t, j):
            eps = 1e-8
            return Expr_if(
                b.lambda_comp[t, j] > 0.95,
                0.984 * ((1 - b.lambda_comp[t, j]) / b.lambda_comp[t, j]) ** (5 / 2),
                (
                    1
                    + 9.0 / 8.0 * b.lambda_comp[t, j] * log(b.lambda_comp[t, j])
                    - 1.56034 * b.lambda_comp[t, j]
                    + 0.528155 * b.lambda_comp[t, j] ** 2
                    + 1.91521 * b.lambda_comp[t, j] ** 3
                    - 2.81903 * b.lambda_comp[t, j] ** 4
                    + 0.270788 * b.lambda_comp[t, j] ** 5
                    + 1.10115 * b.lambda_comp[t, j] ** 6
                    - 0.435933 * b.lambda_comp[t, j] ** 7
                )
                / (1 - b.lambda_comp[t, j] + eps) ** 2,
            )
            # Relationship used by Geraldes & Alves
            # (1 - 2.3 * b.lambda_comp[t, j]
            #     + 1.154 * b.lambda_comp[t, j] ** 2
            #     + 0.224 * b.lambda_comp[t, j] ** 3
            #     )

        @self.Expression(
            self.flowsheet().config.time, solute_set, doc="Pore diffusion coefficient"
        )
        def diffus_pore_comp(b, t, j):
            return (
                b.hindrance_factor_diffusive_comp[t, j]
                * b.feed_side.properties_in[t].diffus_phase_comp["Liq", j]
            )

        @self.Expression(
            self.flowsheet().config.time,
            solute_set,
            doc="Convective hindered transport coefficient",
        )
        def hindrance_factor_convective_comp(b, t, j):
            return (
                1
                + 3.867 * b.lambda_comp[t, j]
                - 1.907 * b.lambda_comp[t, j] ** 2
                - 0.834 * b.lambda_comp[t, j] ** 3
            ) / (1 + 1.867 * b.lambda_comp[t, j] - 0.741 * b.lambda_comp[t, j] ** 2)

        # TODO: some conflict between studies in the literature on whether it should be
        # 1 / b.feed_side.properties_in[t].dielectric_constant - 1 / b.dielectric_constant_pore[t]
        # OR
        # 1 / b.dielectric_constant_pore[t] -1 / b.feed_side.properties_in[t].dielectric_constant
        # Choosing the more commonly used form for now.
        @self.Expression(
            self.flowsheet().config.time, solute_set, doc="Steric partitioning factor"
        )
        def partition_factor_steric_comp(b, t, j):
            return (1 - b.lambda_comp[t, j]) ** 2

        @self.Expression(
            self.flowsheet().config.time,
            solute_set,
            doc="Gibbs free energy of solvation for each ion",
        )
        def gibbs_solvation_comp(b, t, j):
            return (
                b.feed_side.properties_in[t].charge_comp[j] ** 2
                * Constants.elemental_charge**2
                / (
                    8
                    * Constants.pi
                    * Constants.vacuum_electric_permittivity
                    * b.feed_side.properties_in[t].radius_stokes_comp[j]
                )
                * (
                    -1 / b.feed_side.properties_in[t].dielectric_constant
                    + 1 / b.dielectric_constant_pore[t]
                )
            )

        @self.Expression(
            self.flowsheet().config.time,
            solute_set,
            doc="Born solvation contribution to partitioning",
        )
        def partition_factor_born_solvation_comp(b, t, j):
            return exp(
                -b.gibbs_solvation_comp[t, j]
                / (
                    Constants.boltzmann_constant
                    * b.feed_side.properties_in[t].temperature
                )
            )

        @self.Expression(
            self.flowsheet().config.time,
            io_list,
            solute_set,
            doc="Donnan exclusion contribution to partitioning on feed side",
        )
        def partition_factor_donnan_comp_feed(b, t, x, j):
            return exp(
                -b.feed_side.properties_in[t].charge_comp[j]
                * Constants.faraday_constant
                / (Constants.gas_constant * b.pore_entrance[t, x].temperature)
                * b.electric_potential[t, x, "pore_entrance"]
            )

        @self.Expression(
            self.flowsheet().config.time,
            io_list,
            solute_set,
            doc="Donnan exclusion contribution to partitioning on permeate side",
        )
        def partition_factor_donnan_comp_permeate(b, t, x, j):
            return exp(
                -b.feed_side.properties_in[t].charge_comp[j]
                * Constants.faraday_constant
                / (Constants.gas_constant * b.pore_exit[t, x].temperature)
                * (
                    b.electric_potential[t, x, "pore_exit"]
                    - b.electric_potential[t, x, "permeate"]
                )
            )

        # Volumetric Water Flux at inlet and outlet ------------------------------------#
        @self.Expression(
            self.flowsheet().config.time,
            io_list,
            doc="Volumetric water flux at inlet and outlet",
        )
        def flux_vol_water(b, t, x):
            if not x:
                prop = b.feed_side.properties_in[t]
            elif x:
                prop = b.feed_side.properties_out[t]
            return (
                b.flux_mol_phase_comp[t, x, "Liq", "H2O"]
                * prop.mw_comp["H2O"]
                / prop.dens_mass_solvent
            )

        # Average Volumetric Water Flux ------------------------------------#
        @self.Expression(
            self.flowsheet().config.time, doc="Average volumetric water flux"
        )
        def flux_vol_water_avg(b, t):
            return sum(b.flux_vol_water[t, x] for x in io_list) / len(io_list)

        # Average mole flux of each component ------------------------------------#
        @self.Expression(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Average molar component flux",
        )
        def flux_mol_phase_comp_avg(b, t, p, j):
            return sum(b.flux_mol_phase_comp[t, x, p, j] for x in io_list) / len(
                io_list
            )

        # Average concentration inside the membrane------------------------------------#
        @self.Expression(
            self.flowsheet().config.time,
            io_list,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Average molar concentration inside the membrane",
        )
        def conc_mol_phase_comp_pore_avg(b, t, x, p, j):
            return (
                b.pore_entrance[t, x].conc_mol_phase_comp[p, j]
                + b.pore_exit[t, x].conc_mol_phase_comp[p, j]
            ) / len(io_list)

        # OBSERVED rejection of each ion ------------------------------------#
        @self.Expression(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            solute_set,
            doc="Observed solute rejection",
        )
        def rejection_observed_phase_comp(b, t, p, j):
            return (
                1
                - b.mixed_permeate[t].conc_mol_phase_comp[p, j]
                / b.feed_side.properties_in[t].conc_mol_phase_comp[p, j]
            )

        # TODO - no relationship described between mixing length and spacer mixing efficiency with spacer porosity.
        #  Need effective cross-sectional area for velocity at inlet AND outlet. Assuming spacer porosity as an
        #  additional variable in the model that is independent of aforementioned parameters which are used in
        #  the mass transfer coefficient calculation for spiral wound modules. Revisit later.

        # Cross sectional area ------------------------------------#
        @self.Expression(doc="Cross-sectional area")
        def area_cross(b):
            return b.channel_height * b.width * b.spacer_porosity

    def initialize_build(
        self,
        initialize_guess=None,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
        fail_on_warning=False,
        ignore_dof=False,
        automate_rescale=True,
    ):

        """
        General wrapper for pressure changer initialization routines

        Keyword Arguments:
            initialize_guess : a dict of guesses for ....
                   #TODO: revise as appropriate
                   solvent_recovery, solute_recovery,
                   and cp_modulus. These guesses offset the initial values
                   for the retentate, permeate, and membrane interface
                   state blocks from the inlet feed
                   (default =
                   {'deltaP': -1e4,
                   'solvent_recovery': 0.5,
                   'solute_recovery': 0.01,
                   'cp_modulus': 1.1})
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None)
            solver : str indicating which solver to use during
                     initialization (default = None)
            fail_on_warning : boolean argument to fail or only produce  warning upon unsuccessful solve (default=False)
            ignore_dof : boolean argument to ignore when DOF != 0 (default=False)
            automate_rescale: boolean argument to automatically rescale poorly scaled vars
        Returns:
            None
        """

        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")
        # Set solver options
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Extract initial state of inlet feed
        source = self.feed_side.properties_in[self.flowsheet().config.time.first()]
        state_args = self._get_state_args(
            source, self.mixed_permeate[0], initialize_guess, state_args
        )

        # ---------------------------------------------------------------------
        # Initialize holdup block
        flags_feed_side = self.feed_side.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args["feed_side"],
            hold_state=True,
        )
        init_log.info_high("Initialization Step 1 Complete.")
        if not ignore_dof:
            check_dof(self, fail_flag=fail_on_warning, logger=init_log)
        # ---------------------------------------------------------------------
        # Initialize other state blocks based on properties at
        # inlet state block
        self.feed_side.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args["retentate"],
        )
        self.feed_side.properties_interface.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args["interface_in"],
        )
        self.permeate_side.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args["permeate"],
        )
        self.mixed_permeate.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args["permeate"],
        )
        self.pore_entrance.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args["interface_in"],
        )
        self.pore_exit.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args["permeate"],
        )
        init_log.info_high("Initialization Step 2 Complete.")

        # Provide better guesses for unit model variables
        # N_Pe_comp
        if (
            self.config.mass_transfer_coefficient
            == MassTransferCoefficient.spiral_wound
        ):
            for (t, x, j), _ in self.eq_N_Pe_comp.items():
                if not self.N_Pe_comp[t, x, j].is_fixed():
                    self.N_Pe_comp[t, x, j].set_value(
                        value(
                            2
                            * self.channel_height
                            * self.velocity[t, x]
                            / self.config.property_package.diffus_phase_comp["Liq", j]
                        )
                    )
                if not self.Kf_comp[t, x, j].is_fixed():
                    self.Kf_comp[t, x, j].set_value(
                        value(
                            0.753
                            * (
                                self.spacer_mixing_efficiency
                                / (2 - self.spacer_mixing_efficiency)
                            )
                            ** 0.5
                            * (
                                2
                                * self.config.property_package.diffus_phase_comp[
                                    "Liq", j
                                ]
                                / self.channel_height
                            )
                            * self.N_Sc_comp[t, x, j] ** (-1 / 6)
                            * (
                                self.N_Pe_comp[t, x, j]
                                * self.channel_height
                                / (2 * self.spacer_mixing_length)
                            )
                            ** 0.5
                        )
                    )
        for (t, x), _ in self.eq_pressure_permeate_io.items():
            if not self.permeate_side[0, 0].pressure.is_fixed():
                self.permeate_side[0, 0].pressure.set_value(
                    value(self.mixed_permeate[0].pressure)
                )
            if not self.permeate_side[0, 1].pressure.is_fixed():
                self.permeate_side[0, 1].pressure.set_value(
                    value(self.mixed_permeate[0].pressure)
                )
        for j in self.config.property_package.component_list:
            if not self.feed_side.mass_transfer_term[0.0, "Liq", j].is_fixed():
                self.feed_side.mass_transfer_term[0.0, "Liq", j].set_value(
                    value(
                        self.feed_side.properties_out[0].flow_mol_phase_comp["Liq", j]
                        - self.feed_side.properties_in[0].flow_mol_phase_comp["Liq", j]
                    )
                )
            if not self.flux_mol_phase_comp[0, 0, "Liq", j].is_fixed():
                self.flux_mol_phase_comp[0, 0, "Liq", j].set_value(
                    value(-self.feed_side.mass_transfer_term[0.0, "Liq", j] / self.area)
                )
            if not self.flux_mol_phase_comp[0, 1, "Liq", j].is_fixed():
                self.flux_mol_phase_comp[0, 1, "Liq", j].set_value(
                    value(-self.feed_side.mass_transfer_term[0.0, "Liq", j] / self.area)
                )
        # self.report()
        # Double-check for poorly scaled variables after state block initialization
        # and rescale them so that scaled variable values = 1:
        if automate_rescale:
            badly_scaled_vars = list(iscale.badly_scaled_var_generator(self))
            if len(badly_scaled_vars) > 0:
                init_log.warning(
                    f"{len(badly_scaled_vars)} poorly scaled "
                    f"variable(s) will be rescaled so that each scaled variable value = 1"
                )
                [print(i[0], i[1]) for i in badly_scaled_vars]
                self._automate_rescale_variables()
        # ---------------------------------------------------------------------
        # Solve unit attempt 1
        # self.eq_solute_solvent_flux.deactivate()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
            if not check_optimal_termination(res):
                init_log.warning(
                    "Trouble solving NanofiltrationDSPMDE0D unit model with deactivated constraint."
                )
                if automate_rescale:
                    badly_scaled_vars = list(iscale.badly_scaled_var_generator(self))
                    if len(badly_scaled_vars) > 0:
                        init_log.warning(
                            f"{len(badly_scaled_vars)} poorly scaled "
                            f"variable(s) will be rescaled so that each scaled variable value = 1"
                        )
                        [print(i[0], i[1]) for i in badly_scaled_vars]
                        self._automate_rescale_variables()
                # Solve unit attempt 2
                res = opt.solve(self, tee=slc.tee)
                if not check_optimal_termination(res):
                    init_log.warning(
                        "Trouble solving NanofiltrationDSPMDE0D unit model. Trying one more time."
                    )
                    if automate_rescale:
                        badly_scaled_vars = list(
                            iscale.badly_scaled_var_generator(self)
                        )
                        if len(badly_scaled_vars) > 0:
                            init_log.warning(
                                f"{len(badly_scaled_vars)} poorly scaled "
                                f"variable(s) will be rescaled so that each scaled variable value = 1"
                            )
                            [print(i[0], i[1]) for i in badly_scaled_vars]
                            self._automate_rescale_variables()
                    # Solve unit attempt 3
                    res = opt.solve(self, tee=slc.tee)
                    if not check_optimal_termination(res):
                        raise InitializationError(
                            "The property package failed to solve during initialization."
                        )
                    else:
                        init_log.info_high(
                            "Initialization Solve successful on 3rd attempt."
                        )
                else:
                    init_log.info_high(
                        "Initialization Solve successful on 2nd attempt."
                    )
            else:
                init_log.info_high("Initialization Solve successful on 1st attempt.")

        # Release Inlet state
        self.feed_side.release_state(flags_feed_side, outlvl)
        # # Rescale any badly scaled vars
        # if automate_rescale:
        #     badly_scaled_vars = list(iscale.badly_scaled_var_generator(self))
        #     if len(badly_scaled_vars) > 0:
        #         init_log.warn(
        #             f"After solve: {len(badly_scaled_vars)} poorly scaled "
        #             f"variable(s) will be rescaled so that each scaled variable value = 1"
        #         )
        #     self._automate_rescale_variables()
        init_log.info(f"Initialization Complete: {idaeslog.condition(res)}")

    def _get_performance_contents(self, time_point=0):
        # TODO: replace 0 with time_point
        var_dict = {}
        expr_dict = {}
        var_dict["Volumetric Recovery Rate"] = self.recovery_vol_phase[
            time_point, "Liq"
        ]
        var_dict["Membrane Area"] = self.area
        if hasattr(self, "deltaP"):
            var_dict["Pressure Change"] = self.deltaP[time_point]
        # Volume flowrates
        if self.feed_side.properties_in[time_point].is_property_constructed(
            "flow_vol_phase"
        ):
            var_dict[
                f"Volumetric flow rate @ feed inlet"
            ] = self.feed_side.properties_in[time_point].flow_vol_phase["Liq"]
        if self.feed_side.properties_interface[time_point, 0].is_property_constructed(
            "flow_vol_phase"
        ):
            var_dict[
                f"Volumetric flow rate @ inlet interface"
            ] = self.feed_side.properties_interface[time_point, 0].flow_vol_phase["Liq"]
        if self.pore_entrance[time_point, 0].is_property_constructed("flow_vol_phase"):
            var_dict[
                f"Volumetric flow rate @ inlet pore entrance"
            ] = self.pore_entrance[time_point, 0].flow_vol_phase["Liq"]
        if self.pore_exit[time_point, 0].is_property_constructed("flow_vol_phase"):
            var_dict[f"Volumetric flow rate @ inlet pore exit"] = self.pore_exit[
                time_point, 0
            ].flow_vol_phase["Liq"]
        if self.permeate_side[time_point, 0].is_property_constructed("flow_vol_phase"):
            var_dict[f"Volumetric flow rate @ inlet permeate"] = self.permeate_side[
                time_point, 0
            ].flow_vol_phase["Liq"]

        if self.feed_side.properties_out[time_point].is_property_constructed(
            "flow_vol_phase"
        ):
            var_dict[
                f"Volumetric flow rate @ feed outlet"
            ] = self.feed_side.properties_out[time_point].flow_vol_phase["Liq"]
            if self.feed_side.properties_interface[
                time_point, 1
            ].is_property_constructed("flow_vol_phase"):
                var_dict[
                    f"Volumetric flow rate @ outlet interface"
                ] = self.feed_side.properties_interface[time_point, 1].flow_vol_phase[
                    "Liq"
                ]
        if self.pore_entrance[time_point, 1].is_property_constructed("flow_vol_phase"):
            var_dict[
                f"Volumetric flow rate @ outlet pore entrance"
            ] = self.pore_entrance[time_point, 1].flow_vol_phase["Liq"]
        if self.pore_exit[time_point, 1].is_property_constructed("flow_vol_phase"):
            var_dict[f"Volumetric flow rate @ outlet pore exit"] = self.pore_exit[
                time_point, 1
            ].flow_vol_phase["Liq"]
        if self.permeate_side[time_point, 0].is_property_constructed("flow_vol_phase"):
            var_dict[f"Volumetric flow rate @ outlet permeate"] = self.permeate_side[
                time_point, 1
            ].flow_vol_phase["Liq"]
        var_dict[f"Volumetric flow rate @ final permeate"] = self.mixed_permeate[
            time_point
        ].flow_vol_phase["Liq"]

        expr_dict["Average Volumetric Flux (LMH)"] = (
            self.flux_vol_water_avg[time_point] * 3.6e6
        )
        for j in self.config.property_package.component_list:
            expr_dict[f"Average Mole FLux of {j} "] = self.flux_mol_phase_comp_avg[
                time_point, "Liq", j
            ]
            # Molar flowrates
            var_dict[
                f"Molar flow rate of {j} @ feed inlet"
            ] = self.feed_side.properties_in[time_point].flow_mol_phase_comp["Liq", j]
            var_dict[
                f"Molar flow rate of {j} @ feed outlet"
            ] = self.feed_side.properties_out[time_point].flow_mol_phase_comp["Liq", j]
            var_dict[
                f"Molar flow rate of {j} @ membrane-interface inlet"
            ] = self.feed_side.properties_interface[time_point, 0].flow_mol_phase_comp[
                "Liq", j
            ]
            var_dict[
                f"Molar flow rate of {j} @ membrane-interface outlet"
            ] = self.feed_side.properties_interface[time_point, 1].flow_mol_phase_comp[
                "Liq", j
            ]
            var_dict[
                f"Molar flow rate of {j} @ pore entrance, inlet"
            ] = self.pore_entrance[time_point, 0].flow_mol_phase_comp["Liq", j]
            var_dict[
                f"Molar flow rate of {j} @ pore entrance, outlet"
            ] = self.pore_entrance[time_point, 1].flow_mol_phase_comp["Liq", j]
            var_dict[f"Molar flow rate of {j} @ pore exit, inlet"] = self.pore_exit[
                time_point, 0
            ].flow_mol_phase_comp["Liq", j]
            var_dict[f"Molar flow rate of {j} @ pore exit, outlet"] = self.pore_exit[
                time_point, 1
            ].flow_mol_phase_comp["Liq", j]
            var_dict[f"Molar flow rate of {j} @ permeate, inlet"] = self.permeate_side[
                time_point, 0
            ].flow_mol_phase_comp["Liq", j]
            var_dict[f"Molar flow rate of {j} @ permeate, outlet"] = self.permeate_side[
                time_point, 1
            ].flow_mol_phase_comp["Liq", j]
            var_dict[f"Molar flow rate of {j} @ mixed permeate"] = self.mixed_permeate[
                time_point
            ].flow_mol_phase_comp["Liq", j]

            var_dict[
                f"Molar Concentration of {j} @ Feed Inlet"
            ] = self.feed_side.properties_in[time_point].conc_mol_phase_comp["Liq", j]
            var_dict[
                f"Molar Concentration of {j} @ Feed Outlet"
            ] = self.feed_side.properties_out[time_point].conc_mol_phase_comp["Liq", j]
            var_dict[
                f"Molar Concentration of {j} @ Final Permeate"
            ] = self.mixed_permeate[time_point].conc_mol_phase_comp["Liq", j]

            for x in self.io_list:
                if not x:
                    io = "Inlet"
                    prop_feed = self.feed_side.properties_in[0]
                elif x:
                    io = "Outlet"
                    prop_feed = self.feed_side.properties_out[0]

                var_dict[
                    f"Molar Concentration of {j} @ Membrane Interface, {io}"
                ] = self.feed_side.properties_interface[
                    time_point, x
                ].conc_mol_phase_comp[
                    "Liq", j
                ]
                var_dict[
                    f"Molar Concentration of {j} @ Pore Entrance, {io}"
                ] = self.pore_entrance[time_point, x].conc_mol_phase_comp["Liq", j]
                var_dict[
                    f"Molar Concentration of {j} @ Pore Exit, {io}"
                ] = self.pore_exit[time_point, x].conc_mol_phase_comp["Liq", j]
                var_dict[
                    f"Molar Concentration of {j} @ Permeate, {io}"
                ] = self.permeate_side[time_point, x].conc_mol_phase_comp["Liq", j]

        for j in (
            self.config.property_package.solute_set
            | self.config.property_package.ion_set
        ):
            expr_dict[f"Stokes radius of {j}"] = self.feed_side.properties_in[
                time_point
            ].radius_stokes_comp[j]
            expr_dict[f"Stokes:Pore Radius Ratio of {j}"] = self.lambda_comp[
                time_point, j
            ]
            expr_dict[
                f"Diffusive Hindrance Factor of {j}"
            ] = self.hindrance_factor_diffusive_comp[time_point, j]
            expr_dict[
                f"Convective Hindrance Factor of {j}"
            ] = self.hindrance_factor_convective_comp[time_point, j]
            expr_dict[f"Pore Diffusivity of {j}"] = self.diffus_pore_comp[time_point, j]

            expr_dict[
                f"Gibbs Free Energy of Solvation of {j}"
            ] = self.gibbs_solvation_comp[time_point, j]
            expr_dict[
                f"Born Solvation Energy Partitioning Factor of {j}"
            ] = self.partition_factor_born_solvation_comp[time_point, j]
            expr_dict[
                f"Steric Hindrance Partitioning Factor of {j}"
            ] = self.partition_factor_steric_comp[time_point, j]
            expr_dict[
                f"Donnan Partitioning Factor of {j} @ Feed-side Inlet"
            ] = self.partition_factor_donnan_comp_feed[time_point, 0, j]
            expr_dict[
                f"Donnan Partitioning Factor of {j} @ Permeate-side Inlet"
            ] = self.partition_factor_donnan_comp_permeate[time_point, 0, j]
            expr_dict[
                f"Donnan Partitioning Factor of {j} @ Feed-side Outlet"
            ] = self.partition_factor_donnan_comp_feed[time_point, 1, j]
            expr_dict[
                f"Donnan Partitioning Factor of {j} @ Permeate-side Outlet"
            ] = self.partition_factor_donnan_comp_permeate[time_point, 1, j]

            var_dict[
                f"Intrinsic Rejection of {j}"
            ] = self.rejection_intrinsic_phase_comp[time_point, "Liq", j]

            if self.feed_side.properties_in[time_point].is_property_constructed(
                "pressure_osm_phase"
            ):
                var_dict[
                    f"Osmotic Pressure @ Bulk Feed, Inlet (Pa)"
                ] = self.feed_side.properties_in[time_point].pressure_osm_phase["Liq"]
            if self.feed_side.properties_out[time_point].is_property_constructed(
                "pressure_osm_phase"
            ):
                var_dict[
                    f"Osmotic Pressure @ Bulk Feed, Outlet (Pa)"
                ] = self.feed_side.properties_out[time_point].pressure_osm_phase["Liq"]

            for x in self.io_list:
                if not x:
                    io = "Inlet"
                    prop_feed = self.feed_side.properties_in[0]
                elif x:
                    io = "Outlet"
                    prop_feed = self.feed_side.properties_out[0]

                var_dict[
                    f"Osmotic Pressure @ Membrane Interface, {io} (Pa)"
                ] = self.feed_side.properties_interface[
                    time_point, x
                ].pressure_osm_phase[
                    "Liq"
                ]

                var_dict[
                    f"Osmotic Pressure @ Permeate, {io} (Pa)"
                ] = self.permeate_side[time_point, x].pressure_osm_phase["Liq"]
                expr_dict[f"Net Driving Pressure, {io} (Pa)"] = (
                    prop_feed.pressure
                    - self.permeate_side[0, x].pressure
                    - (
                        self.feed_side.properties_interface[0, x].pressure_osm_phase[
                            "Liq"
                        ]
                        - self.permeate_side[0, x].pressure_osm_phase["Liq"]
                    )
                )
                var_dict[
                    f"Electric Potential @ Pore Entrance, {io}"
                ] = self.electric_potential[0, x, "pore_entrance"]
                var_dict[
                    f"Electric Potential @ Pore Exit, {io}"
                ] = self.electric_potential[0, x, "pore_exit"]
                var_dict[
                    f"Electric Potential @ Permeate, {io}"
                ] = self.electric_potential[0, x, "permeate"]
                if hasattr(self, "electric_potential_grad_feed_interface"):
                    var_dict[
                        f"Electric Potential Gradient @ Feed-Membrane Interface, {io}"
                    ] = self.electric_potential_grad_feed_interface[0, x]

        return {"vars": var_dict, "exprs": expr_dict}

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
        This method returns state_args (initial guesses for state variable values) to be passed to
        each stateblock during initialization.
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
        # TODO: enable deltaP guess when pressure drop is added
        if "deltaP" not in initialize_guess:
            initialize_guess["deltaP"] = 0
        if "solvent_recovery" not in initialize_guess:
            initialize_guess["solvent_recovery"] = 0.1
        if "solute_recovery" not in initialize_guess:
            initialize_guess["solute_recovery"] = 0.1
        if "cp_modulus" not in initialize_guess:
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

        if "flow_mol_phase_comp" not in state_args.keys():
            raise ConfigurationError(
                f"{self.__class__.__name__} initialization routine expects "
                "flow_mol_phase_comp as a state variable. Check "
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
            state_args_retentate["flow_mol_phase_comp"][("Liq", j)] *= (
                1 - initialize_guess["solvent_recovery"]
            )
            state_args_permeate["flow_mol_phase_comp"][("Liq", j)] *= initialize_guess[
                "solvent_recovery"
            ]
        for j in (
            self.config.property_package.solute_set
            | self.config.property_package.ion_set
        ):
            state_args_retentate["flow_mol_phase_comp"][("Liq", j)] *= (
                1 - initialize_guess["solute_recovery"]
            )
            state_args_permeate["flow_mol_phase_comp"][("Liq", j)] *= initialize_guess[
                "solute_recovery"
            ]

        state_args_interface_in = deepcopy(state_args)
        state_args_interface_out = deepcopy(state_args_retentate)

        for j in (
            self.config.property_package.solute_set
            | self.config.property_package.ion_set
        ):
            state_args_interface_in["flow_mol_phase_comp"][
                ("Liq", j)
            ] *= initialize_guess["cp_modulus"]
            state_args_interface_out["flow_mol_phase_comp"][
                ("Liq", j)
            ] *= initialize_guess["cp_modulus"]

        return {
            "feed_side": state_args,
            "retentate": state_args_retentate,
            "permeate": state_args_permeate,
            "interface_in": state_args_interface_in,
            "interface_out": state_args_interface_out,
        }

    # stateblock properties need to rescale solute values by a larger factor
    def _rescale_sb_variable(self, var, factor=100):
        if var not in self._sb_scaled_properties:
            sf = iscale.get_scaling_factor(var)
            iscale.set_scaling_factor(var, sf * factor)
            self._sb_scaled_properties.add(var)

    # automatically rescale poorly scaled variables by setting a new scaling factor
    # which multiplies a variable value by the old scaling factor divided by the poorly scaled (resulting) value,
    # bringing the new scaled value to 1. Providing a rescale_factor would just multiply that factor by 1.
    def _automate_rescale_variables(self, rescale_factor=None):
        if rescale_factor is None:
            rescale_factor = 1
        for var, sv in iscale.badly_scaled_var_generator(self):
            if iscale.get_scaling_factor(var) is None:
                print(f"{var} is missing a scaling factor")
                continue
            sf = iscale.get_scaling_factor(var)
            iscale.set_scaling_factor(var, sf / sv * rescale_factor)
            iscale.calculate_scaling_factors(self)

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # setting scaling factors for variables
        for v in self.recovery_vol_phase.values():
            iscale.set_scaling_factor(v, 1)
        if iscale.get_scaling_factor(self.radius_pore) is None:
            iscale.set_scaling_factor(self.radius_pore, 10 / value(self.radius_pore))
        if iscale.get_scaling_factor(self.membrane_thickness_effective) is None:
            iscale.set_scaling_factor(self.membrane_thickness_effective, 1e7)

        # setting scaling factors for variables
        # these variables should have user input, if not there will be a warning
        if iscale.get_scaling_factor(self.area) is None:
            sf = iscale.get_scaling_factor(self.area, default=1, warning=True)
            iscale.set_scaling_factor(self.area, sf)

        for v in self.electric_potential.values():
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, 1e4)
        if hasattr(self, "electric_potential_grad_feed_interface"):
            for v in self.electric_potential_grad_feed_interface.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e8)

        # these variables do not typically require user input,
        # will not override if the user does provide the scaling factor
        for (t, x, p, j), v in self.flux_mol_phase_comp.items():
            if iscale.get_scaling_factor(v) is None:
                comp = self.config.property_package.get_component(j)
                if comp.is_solvent():
                    if x == 0:
                        prop_feed = self.feed_side.properties_in[t]
                    elif x == 1:
                        prop_feed = self.feed_side.properties_out[t]

                    sf = (
                        iscale.get_scaling_factor(prop_feed.dens_mass_phase["Liq"])
                        / iscale.get_scaling_factor(prop_feed.mw_comp[j])
                        * iscale.get_scaling_factor(prop_feed.pressure)
                        * iscale.get_scaling_factor(self.radius_pore) ** 2
                        / iscale.get_scaling_factor(prop_feed.visc_d_phase["Liq"])
                        / iscale.get_scaling_factor(self.membrane_thickness_effective)
                    )
                    iscale.set_scaling_factor(v, sf)

                if comp.is_solute():
                    # Todo: revisit later
                    sf = (
                        iscale.get_scaling_factor(
                            self.flux_mol_phase_comp[t, x, "Liq", "H2O"]
                        )
                        / iscale.get_scaling_factor(
                            self.feed_side.properties_in[t].dens_mass_phase["Liq"]
                        )
                        * iscale.get_scaling_factor(
                            self.feed_side.properties_in[t].mw_comp[j]
                        )
                        * iscale.get_scaling_factor(
                            self.permeate_side[t, x].conc_mol_phase_comp["Liq", j]
                        )
                    )
                    # sf = 1e5
                    iscale.set_scaling_factor(v, sf)

        for v in self.rejection_intrinsic_phase_comp.values():
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, 1e1)

        for j in self.config.property_package.component_list:
            if (
                j
                in self.config.property_package.solute_set
                | self.config.property_package.ion_set
            ):
                iscale.set_scaling_factor(
                    self.feed_side.mass_transfer_term[0, "Liq", j], 1e4
                )
            else:
                iscale.set_scaling_factor(
                    self.feed_side.mass_transfer_term[0, "Liq", j], 1
                )
        if hasattr(self, "Kf_comp"):
            for v in self.Kf_comp.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e6)

        if hasattr(self, "N_Pe_comp"):
            for v in self.N_Pe_comp.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e-4)
        if hasattr(self, "N_Sc_comp"):
            for v in self.N_Sc_comp.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e-3)

        for v in self.velocity.values():
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, 1e1)

        if iscale.get_scaling_factor(self.channel_height) is None:
            iscale.set_scaling_factor(self.channel_height, 1e3)

        if hasattr(self, "length"):
            if iscale.get_scaling_factor(self.length) is None:
                iscale.set_scaling_factor(self.length, 1)

        if hasattr(self, "width"):
            if iscale.get_scaling_factor(self.width) is None:
                iscale.set_scaling_factor(self.width, 1)

        # transforming constraints
        for (t, x, p), con in self.eq_water_flux.items():
            sf = (
                iscale.get_scaling_factor(self.flux_mol_phase_comp[t, x, p, "H2O"])
                * iscale.get_scaling_factor(
                    self.feed_side.properties_in[t].mw_comp["H2O"]
                )
                / iscale.get_scaling_factor(
                    self.feed_side.properties_in[t].dens_mass_solvent
                )
            )
            iscale.constraint_scaling_transform(con, sf / 100)

        # for (t, x, p, j), con in self.eq_solute_solvent_flux.items():
        #     # todo: revisit sf
        #     # sf = (iscale.get_constraint_transform_applied_scaling_factor(self.eq_water_flux[t, x, p])
        #     #       * iscale.get_scaling_factor(self.mixed_permeate[t].conc_mol_phase_comp[p, j]))
        #     iscale.constraint_scaling_transform(con, 1e5)
        #
        # for (t, x, p, j), con in self.eq_solute_flux_concentration_polarization.items():
        #     # todo: revisit sf
        #     # sf = (iscale.get_constraint_transform_applied_scaling_factor(self.eq_water_flux[t, x, p])
        #     #       * iscale.get_scaling_factor(self.mixed_permeate[t].conc_mol_phase_comp[p, j]))
        #     iscale.constraint_scaling_transform(con, 1e2)
        #
        # for con in self.eq_solute_flux_pore_domain.values():
        #     # todo: revisit sf
        #     # sf = (iscale.get_constraint_transform_applied_scaling_factor(self.eq_water_flux[t, x, p])
        #     #       * iscale.get_scaling_factor(self.mixed_permeate[t].conc_mol_phase_comp[p, j]))
        #     iscale.constraint_scaling_transform(con, 1e5)
        #
        for con in self.eq_solute_solvent_flux.values():
            iscale.constraint_scaling_transform(con, 1e1)
        for con in self.eq_solute_flux_pore_domain.values():
            iscale.constraint_scaling_transform(con, 1e1)

        for con in self.eq_solute_flux_concentration_polarization.values():
            iscale.constraint_scaling_transform(con, 1e3)

        for con in self.feed_side.material_balances.values():
            iscale.constraint_scaling_transform(con, 1e-1)

        for con in self.eq_interfacial_partitioning_feed.values():
            iscale.constraint_scaling_transform(con, 1e0)
        for key in self.eq_interfacial_partitioning_feed.keys():
            if key[-1] == "Cl_-":
                iscale.constraint_scaling_transform(
                    self.eq_interfacial_partitioning_feed[key], 1e0
                )
        for con in self.eq_interfacial_partitioning_permeate.values():
            iscale.constraint_scaling_transform(con, 1e-1)
        for con in self.eq_electroneutrality_pore.values():
            iscale.constraint_scaling_transform(con, 1e-4)
        for con in self.eq_electroneutrality_permeate.values():
            iscale.constraint_scaling_transform(con, 1e-5)
        if hasattr(self, "eq_electroneutrality_interface"):
            for con in self.eq_electroneutrality_interface.values():
                iscale.constraint_scaling_transform(con, 1e-3)
        for con in self.eq_permeate_production.values():
            iscale.constraint_scaling_transform(con, 1e-2)
        for con in self.eq_rejection_intrinsic_phase_comp.values():
            iscale.constraint_scaling_transform(con, 1e-1)
        if hasattr(self, "eq_N_Pe_comp"):
            for con in self.eq_N_Pe_comp.values():
                iscale.constraint_scaling_transform(con, 1e-9)

        if hasattr(self, "eq_Kf_comp"):
            for ind, con in self.eq_Kf_comp.items():
                sf = iscale.get_scaling_factor(self.Kf_comp[ind])
                iscale.constraint_scaling_transform(con, sf)

        # Constraints below all scaled by factor of 1 to satisfy test
        # TODO: refine later if improvements in scaling can be made
        if hasattr(self, "eq_N_Sc_comp"):
            for con in self.eq_N_Sc_comp.values():
                iscale.constraint_scaling_transform(con, 1)
        for con in self.eq_area.values():
            iscale.constraint_scaling_transform(con, 1)
        for con in self.eq_permeate_isothermal.values():
            iscale.constraint_scaling_transform(con, 1)
        for con in self.eq_mass_transfer_feed.values():
            iscale.constraint_scaling_transform(con, 1)
        for con in self.eq_recovery_vol_phase.values():
            iscale.constraint_scaling_transform(con, 1)
        for con in self.eq_velocity.values():
            iscale.constraint_scaling_transform(con, 1)
        for con in self.eq_pressure_permeate_io.values():
            iscale.constraint_scaling_transform(con, 1)
        for con in self.eq_pore_isothermal.values():
            iscale.constraint_scaling_transform(con, 1)
        for con in self.eq_permeate_isothermal_mixed.values():
            iscale.constraint_scaling_transform(con, 1)
        for con in self.eq_equal_flow_vol_pore.values():
            iscale.constraint_scaling_transform(con, 1)
        for con in self.eq_equal_flow_vol_permeate.values():
            iscale.constraint_scaling_transform(con, 1)
        for con in self.eq_equal_flow_vol_pore_exit_perm.values():
            iscale.constraint_scaling_transform(con, 1)
        for con in self.feed_side.eq_equal_flow_vol_feed_interface.values():
            iscale.constraint_scaling_transform(con, 1)
        for con in self.feed_side.eq_feed_interface_isothermal.values():
            iscale.constraint_scaling_transform(con, 1)
        for con in self.feed_side.eq_feed_isothermal.values():
            iscale.constraint_scaling_transform(con, 1)
