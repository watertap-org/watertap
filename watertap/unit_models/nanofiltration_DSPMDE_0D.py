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

# Import Pyomo libraries
from pyomo.environ import (Block,
                           Set,
                           Var,
                           Param,
                           Suffix,
                           NonNegativeReals,
                           Reals,
                           Reference,
                           units as pyunits,
                           log,
                           value,
                           Expr_if,
                           Constraint,
                           exp)
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (ControlVolume0DBlock,
                        declare_process_block_class,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        UnitModelBlockData,
                        useDefault,
                        MaterialFlowBasis)
from idaes.core.util import get_solver
from idaes.core.util.math import smooth_min
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
import idaes.core.util.scaling as iscale
from idaes.core.util.constants import Constants

import idaes.logger as idaeslog


_log = idaeslog.getLogger(__name__)


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

    CONFIG.declare("dynamic", ConfigValue(
        domain=In([False]),
        default=False,
        description="Dynamic model flag - must be False",
        doc="""Indicates whether this model will be dynamic or not,
    **default** = False. NF units do not support dynamic
    behavior."""))
    CONFIG.declare("has_holdup", ConfigValue(
        default=False,
        domain=In([False]),
        description="Holdup construction flag - must be False",
        doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False. NF units do not have defined volume, thus
    this must be False."""))
    CONFIG.declare("material_balance_type", ConfigValue(
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
    # CONFIG.declare("energy_balance_type", ConfigValue(
    #     default=EnergyBalanceType.useDefault,
    #     domain=In(EnergyBalanceType),
    #     description="Energy balance construction flag",
    #     doc="""Indicates what type of energy balance should be constructed,
    # **default** - EnergyBalanceType.useDefault.
    # **Valid values:** {
    # **EnergyBalanceType.useDefault - refer to property package for default
    # balance type
    # **EnergyBalanceType.none** - exclude energy balances,
    # **EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
    # **EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
    # **EnergyBalanceType.energyTotal** - single energy balance for material,
    # **EnergyBalanceType.energyPhase** - energy balances for each phase.}"""))
    CONFIG.declare("momentum_balance_type", ConfigValue(
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
    CONFIG.declare("has_pressure_change", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Pressure change term construction flag",
        doc="""Indicates whether terms for pressure change should be
    constructed,
    **default** - False.
    **Valid values:** {
    **True** - include pressure change terms,
    **False** - exclude pressure change terms.}"""))
    CONFIG.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations,
    **default** - useDefault.
    **Valid values:** {
    **useDefault** - use default package from parent model or flowsheet,
    **PhysicalParameterObject** - a PhysicalParameterBlock object.}"""))
    CONFIG.declare("property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
    and used when constructing these,
    **default** - None.
    **Valid values:** {
    see property package for documentation.}"""))


    def _process_config(self):
        if len(self.config.property_package.solvent_set) > 1:
            raise ConfigurationError("NF model only supports one solvent component,"
                                     "the provided property package has specified {} solvent components"
                                     .format(len(self.config.property_package.solvent_set)))

        if len(self.config.property_package.solvent_set) == 0:
            raise ConfigurationError("The NF model was expecting a solvent and did not receive it.")

        if (hasattr(self.config.property_package,'ion_set') and len(self.config.property_package.ion_set) == 0) \
                or (hasattr(self.config.property_package,'solute_set') and len(self.config.property_package.solute_set) == 0):
            raise ConfigurationError("This NF model was expecting ions and did not receive any.")

    def build(self):
        # Call UnitModel.build to setup dynamics
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        units_meta = self.config.property_package.get_metadata().get_derived_units

        self.io_list = io_list = Set(initialize=[0, 1])  # inlet/outlet set

        self._process_config()

        if hasattr(self.config.property_package,'ion_set'):
            solute_set = self.config.property_package.ion_set
        elif hasattr(self.config.property_package,'solute_set'):
            solute_set = self.config.property_package.solute_set

        solvent_set = self.config.property_package.solvent_set
        solvent_solute_set = solvent_set | solute_set
        phase_list = self.config.property_package.phase_list


        # Build control volume for feed side
        self.feed_side = ControlVolume0DBlock(default={
            "dynamic": False,
            "has_holdup": False,
            "property_package": self.config.property_package,
            "property_package_args": self.config.property_package_args})

        self.feed_side.add_state_blocks(
            has_phase_equilibrium=False)

        self.feed_side.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_mass_transfer=True)

        self.feed_side.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change)

        # Make indexed stateblock and separate stateblock for permeate-side and permeate outlet, respectively.
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False  # these blocks are not inlets

        # # Add permeate block
        self.permeate_side = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            io_list,
            doc="Material properties of permeate along permeate channel",
            default=tmp_dict)
        self.mixed_permeate = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of mixed permeate exiting the module",
            default=tmp_dict)

        # Add Ports
        self.add_inlet_port(name='inlet', block=self.feed_side)
        self.add_outlet_port(name='retentate', block=self.feed_side)
        self.add_port(name='permeate', block=self.mixed_permeate)

        # Membrane interface: indexed state block
        self.feed_side.properties_interface = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            io_list,
            doc="Material properties of feed-side membrane interface",
            default=tmp_dict)
        # Pore entrance: indexed state block
        self.pore_entrance = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            io_list,
            doc="Fluid properties within the membrane pore entrance",
            default=tmp_dict)
        # Pore exit: indexed state block
        self.pore_exit = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            io_list,
            doc="Fluid properties within the membrane pore exit",
            default=tmp_dict)

        # References for control volume
        # pressure change
        if (self.config.has_pressure_change is True and
                self.config.momentum_balance_type != 'none'):
            self.deltaP = Reference(self.feed_side.deltaP)
        ###############################################################################################################
        # Parameters
        ###############################################################################################################
        self.tol_electroneutrality = Param(
            initialize=1e-6,
            mutable=True,
            domain=NonNegativeReals,
            units=pyunits.mol/pyunits.m**3,
            doc='Electroneutrality tolerance'
            )
        ###############################################################################################################
        # Variables
        ###############################################################################################################
        #1. Component mole flux, J, DOF=Nj*2 for inlet/outlet
        self.flux_mol_phase_comp = Var(
            self.flowsheet().config.time,
            io_list,
            phase_list,
            solvent_solute_set,
            initialize=lambda b,t,x,p,j : 2.5e-2 if j in solvent_set else 1e-5, #TODO: divide solvent by .02 and solute by .1
            bounds=lambda b,t,x,p,j : (5e-3, 1.5) if j in solvent_set else (1e-7, 1e-2), #TODO: divide solvent by .02 and solute by .1
            units=units_meta('amount')*units_meta('length')**-2*units_meta('time')**-1,
            doc='Component molar flux at inlet and outlet of membrane')

        # 2. Pore radius, rp, DOF = 1
        self.radius_pore = Var(
            initialize=0.5e-9, #TODO: revisit
            domain=NonNegativeReals,
            units=units_meta('length'),
            doc='Membrane pore radius')

        #3. Effective thickness of membrane, Ak, DOF=1
        self.membrane_thickness_effective = Var(
            initialize=1.33e-6, # Value used by Labban et al., 2017
            domain=NonNegativeReals,
            units=units_meta('length'),
            doc='Effective membrane thickness')

        #4. Effective thickness of membrane, Ak, DOF=1
        self.membrane_charge_density = Var(
            self.flowsheet().config.time,
            initialize=-50, # near value used in Roy et al.
            domain=Reals,
            units=pyunits.mol*pyunits.m**-3,
            doc='Membrane charge density')
        self.dielectric_constant_pore = Var(
            self.flowsheet().config.time,
            initialize=42, # near value used in Roy et al.
            bounds=(1, None),
            units=pyunits.dimensionless, # TODO: revisit bounds/domain
            doc='Pore dielectric constant')
        self.electric_potential = Var(
            self.flowsheet().config.time,
            io_list,
            ['pore_entrance','pore_exit', 'permeate'], #TODO: revisit - build in property model w/o constraint?
            initialize=1, #TODO:revisit
            domain=Reals,
            units=pyunits.V,
            doc='Electric potential of pore entrance/exit, and permeate')
        self.electric_potential_grad_feed_interface = Var(
            self.flowsheet().config.time,
            io_list,
            initialize=1, #TODO: revisit
            domain=Reals,
            units= pyunits.V*pyunits.m**-1, # TODO: revisit- Geraldes and Alves give unitless while Roy et al. give V/m
            doc='Electric potential gradient of feed-membrane interface')
        self.Kf_comp = Var(
            self.flowsheet().config.time,
            self.io_list,
            solute_set,
            initialize=5e-5,
            bounds=(1e-6, 1e-3),
            domain=NonNegativeReals,
            units=units_meta('length') * units_meta('time') ** -1,
            doc='Component mass transfer coefficient in feed channel at inlet and outlet')
        self.rejection_phase_comp = Var(
            self.flowsheet().config.time,
            phase_list,
            solute_set,
            initialize=0.9,
            bounds=(-1 + 1e-6, 1 - 1e-6),
            units=pyunits.dimensionless,
            doc='Observed solute rejection')
        self.area = Var(
            initialize=10,
            bounds=(1e-2, 1e3),
            domain=NonNegativeReals,
            units=units_meta('length') ** 2,
            doc='Membrane area')
        # self.length = Var(
        #     initialize=10,
        #     bounds=(0.1, 5e2),
        #     domain=NonNegativeReals,
        #     units=units_meta('length'),
        #     doc='Effective membrane length')
        # self.width = Var(
        #     initialize=1,
        #     bounds=(0.1, 5e2),
        #     domain=NonNegativeReals,
        #     units=units_meta('length'),
        #     doc='Effective feed-channel width')

        def recovery_mol_phase_comp_initialize(b, t, p, j):
            if j in b.config.property_package.solvent_set:
                return 0.8
            elif j in solute_set:
                return 0.1

        def recovery_mol_phase_comp_bounds(b, t, p, j):
            ub = 1 - 1e-6
            if j in b.config.property_package.solvent_set:
                lb = 1e-2
            elif j in solute_set:
                lb = 1e-5
            else:
                lb = 1e-5
            return lb, ub

        self.recovery_mol_phase_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            solvent_solute_set,
            initialize=recovery_mol_phase_comp_initialize,
            bounds=recovery_mol_phase_comp_bounds,
            units=pyunits.dimensionless,
            doc='Mole-based component recovery')
        self.recovery_vol_phase = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            initialize=0.1,
            bounds=(1e-2, 1 - 1e-6),
            units=pyunits.dimensionless,
            doc='Volumetric-based recovery')
        ###############################################################################################################
        # Expressions
        ###############################################################################################################
        # Stokes radius to membrane pore radius ratio (for each solute)
        @self.Expression(self.flowsheet().config.time,
                         solute_set,
                         doc="Ratio of stokes radius to membrane pore radius equation")
        def lambda_comp(b, t, j):
            return smooth_min(1, b.feed_side.properties_in[t].radius_stokes_comp[j]/b.radius_pore)

        @self.Expression(self.flowsheet().config.time,
                         solute_set,
                         doc="Diffusive hindered transport coefficient")
        def hindrance_factor_diffusive_comp(b, t, j):
            eps = 1e-8
            return Expr_if(b.lambda_comp[t,j] > 0.95,
                        0.984 *((1 - b.lambda_comp[t,j])/(b.lambda_comp[t,j])) ** (5/2),
                           (1 + 9. / 8. * b.lambda_comp[t, j] * log(b.lambda_comp[t, j])
                        - 1.56034 * b.lambda_comp[t, j]
                        + 0.528155 * b.lambda_comp[t, j]**2
                        + 1.91521 * b.lambda_comp[t, j]**3
                        - 2.81903 * b.lambda_comp[t, j]**4
                        + 0.270788 * b.lambda_comp[t, j]**5
                        - 1.10115 * b.lambda_comp[t, j]**6
                        - 0.435933 * b.lambda_comp[t, j]**7) /
                        (1 - b.lambda_comp[t, j]+eps) ** 2,
                    )

        @self.Expression(self.flowsheet().config.time,
                         solute_set,
                         doc="Pore diffusion coefficient")
        def diffus_pore_comp(b, t, j):
            return b.hindrance_factor_diffusive_comp[t, j] * b.feed_side.properties_in[t].diffus_phase_comp['Liq', j]

        @self.Expression(self.flowsheet().config.time,
                         solute_set,
                         doc="Convective hindered transport coefficient")
        def hindrance_factor_convective_comp(b, t, j):
            return ((1 + 3.867 * b.lambda_comp[t, j]
                     - 1.907 * b.lambda_comp[t, j] ** 2
                     - 0.834 * b.lambda_comp[t, j] ** 3)
                    / (1 + 1.867 * b.lambda_comp[t, j] - 0.741 * b.lambda_comp[t, j]))

        @self.Expression(self.flowsheet().config.time,
                         solute_set,
                         doc="Steric partitioning factor")
        def partition_factor_steric_comp(b, t, j):
            return (1 - b.lambda_comp[t, j]) ** 2

        @self.Expression(self.flowsheet().config.time,
                         solute_set,
                         doc="Gibbs free energy of solvation for each ion")
        def gibbs_solvation_comp(b, t, j):
            return (b.feed_side.properties_in[t].charge_comp[j] ** 2
                    * Constants.elemental_charge ** 2
                    / (8 * Constants.pi
                       * Constants.vacuum_electric_permittivity
                       * b.feed_side.properties_in[t].radius_stokes_comp[j])
                    * (1 / b.feed_side.properties_in[t].dielectric_constant - 1 / b.dielectric_constant_pore[t]))

        @self.Expression(self.flowsheet().config.time,
                         solute_set,
                         doc="Born solvation contribution to partitioning")
        def partition_factor_born_solvation_comp(b, t, j):
            return (-b.gibbs_solvation_comp[t, j]
                    / (Constants.boltzmann_constant * b.feed_side.properties_in[t].temperature))

        @self.Expression(self.flowsheet().config.time,
                         io_list,
                         solute_set,
                         doc="Donnan exclusion contribution to partitioning on feed side")
        def partition_factor_donnan_comp_feed(b, t, x, j):
            return (exp(-b.feed_side.properties_in[t].charge_comp[j] * Constants.faraday_constant
                    / (Constants.gas_constant * b.pore_entrance[t, x].temperature)
                    * b.electric_potential[t, x, 'pore_entrance']))

        @self.Expression(self.flowsheet().config.time,
                         io_list,
                         solute_set,
                         doc="Donnan exclusion contribution to partitioning on permeate side")
        def partition_factor_donnan_comp_permeate(b, t, x, j):
            return (exp(-b.feed_side.properties_in[t].charge_comp[j] * Constants.faraday_constant
                    / (Constants.gas_constant * b.pore_exit[t, x].temperature) #todo: switch to permeate side temp
                    * (b.electric_potential[t, x, 'pore_exit'] - b.electric_potential[t, x, 'permeate'])))

        # Volumetric Water Flux at inlet and outlet ------------------------------------#
        @self.Expression(self.flowsheet().config.time,
                         io_list,
                         doc="Volumetric water flux at inlet and outlet")
        def flux_vol_water(b, t, x):
            prop = b.feed_side.properties_in[t]
            return b.flux_mol_phase_comp[t, x, 'Liq', 'H2O'] * prop.mw_comp['H2O'] / prop.dens_mass_comp['H2O']

        # Average Volumetric Water Flux ------------------------------------#
        @self.Expression(self.flowsheet().config.time,
                         doc="Average volumetric water flux")
        def flux_vol_water_avg(b, t):
            return sum(b.flux_vol_water[t, x] for x in io_list) * 0.5

        # Average mole flux of each component ------------------------------------#
        @self.Expression(self.flowsheet().config.time,
                         phase_list,
                         solvent_solute_set,
                         doc="Average molar component flux")
        def flux_mol_phase_comp_avg(b, t, p, j):
            return sum(b.flux_mol_phase_comp[t, x, p, j] for x in io_list) * 0.5


        ################################################################################################################
        # Constraints
        ################################################################################################################
        # Membrane area
        # @self.Constraint(doc="Membrane area")
        # def eq_area(b):
        #     return b.area == b.length * b.width

        # 1. Feed-solution/membrane equilibrium, DOF= Nj * 2 for inlet/outlet
        @self.Constraint(self.flowsheet().config.time,
                         io_list,
                         phase_list,
                         solute_set,
                         doc="Interfacial partitioning at feed side of membrane")
        def interfacial_partitioning_feed_eq(b, t, x, p, j):
            return (b.pore_entrance[t, x].act_coeff_phase_comp[p, j] * b.pore_entrance[t, x].conc_mol_phase_comp[p, j] ==
                    b.feed_side.properties_interface[t, x].act_coeff_phase_comp[p, j]
                    * b.feed_side.properties_interface[t, x].conc_mol_phase_comp[p, j]
                    * b.partition_factor_steric_comp[t, j]
                    * b.partition_factor_born_solvation_comp[t, j]
                    * b.partition_factor_donnan_comp_feed[t, x, j])

        # 2. Permeate solution/membrane equilibrium, DOF= Nj * 2 for inlet/outlet
        @self.Constraint(self.flowsheet().config.time,
                         io_list,
                         phase_list,
                         solute_set,
                         doc="Interfacial partitioning at permeate side of membrane")
        def interfacial_partitioning_permeate_eq(b, t, x, p, j):
            return (b.pore_exit[t, x].act_coeff_phase_comp[p, j] * b.pore_exit[t, x].conc_mol_phase_comp[p, j] ==
                    b.permeate_side[t, x].act_coeff_phase_comp[p, j]
                    * b.permeate_side[t, x].conc_mol_phase_comp[p, j]
                    * b.partition_factor_steric_comp[t, j]
                    * b.partition_factor_born_solvation_comp[t, j]
                    * b.partition_factor_donnan_comp_permeate[t, x, j])

        # 3. Feed-solution/membrane electroneutrality, DOF=1 *2 for inlet/outlet: DOF= 2
        @self.Constraint(self.flowsheet().config.time,
                         io_list,
                         phase_list,
                         doc="Electroneutrality at feed-side membrane interface")
        def electroneutrality_interface_eq(b, t, x, p):
            return (abs(sum(b.feed_side.properties_interface[t, x].conc_mol_phase_comp[p, j] *
                        b.feed_side.properties_interface[t, x].charge_comp[j] for j in solute_set))
                    == b.tol_electroneutrality)
            #todo: (1) consider adding tol_electroneutrality as property pack param
            # to match assert_electroneutrality method in dspmde prop pack
            # (2) probe whether inequality constraint with tolerance improves model stability and
            # returns optimal solution or if an equality constraint leads to more consistent solution

        # 4. Charge balance inside the membrane, DOF=N nodes across membrane thickness *2 for inlet/outlet: N=2, DOF=4
        @self.Constraint(self.flowsheet().config.time,
                         io_list,
                         ['pore_entrance', 'pore_exit'],
                         phase_list,
                         doc="Electroneutrality within membrane pore")
        def electroneutrality_pore_eq(b, t, x, y, p):
            if y == 'pore_entrance':
                pore_loc = b.pore_entrance[t, x]
            elif y == 'pore_exit':
                pore_loc = b.pore_exit[t, x]
            return (abs(sum(pore_loc.conc_mol_phase_comp[p, j]
                            * pore_loc.charge_comp[j] for j in solute_set)
                        + b.membrane_charge_density[t]) == b.tol_electroneutrality)

        # 4. Permeate electroneutrality, DOF=1 *2 for inlet/outlet:  DOF=2
        @self.Constraint(self.flowsheet().config.time,
                         io_list,
                         phase_list,
                         doc="Electroneutrality in permeate")
        def electroneutrality_permeate_eq(b, t, x, p):
            return (abs(sum(b.permeate_side[t, x].conc_mol_phase_comp[p, j] *
                            b.permeate_side[t, x].charge_comp[j] for j in solute_set)) == b.tol_electroneutrality)

        # 5. Water flux via Hagen-Poiseuille relationship, DOF= 1 * 2 for inlet/outlet: DOF= 2
        @self.Constraint(self.flowsheet().config.time,
                         io_list,
                         phase_list,
                         doc="Hagen-Poiseuille relationship for water flux across membrane")
        def water_flux_eq(b, t, x, p):
            if x == 0:
                prop_feed = b.feed_side.properties_in[t]
            elif x == 1:
                prop_feed = b.feed_side.properties_out[t]
            prop_perm = b.permeate_side[t, x]
            prop_feed_inter = b.feed_side.properties_interface[t, x]
            return (b.flux_vol_water[t, x] ==
                    (prop_feed.pressure - prop_perm.pressure -
                     (prop_feed_inter.pressure_osm - prop_perm.pressure_osm))
                    * (b.radius_pore ** 2)
                    / (8 * prop_feed.visc_d_phase[p] * b.membrane_thickness_effective)
                    )

        # 6. Unhindered mass transfer; Js,i=Jw*cp,i; DOF= Nj * 2 for inlet/outlet
        @self.Constraint(self.flowsheet().config.time,
                         io_list,
                         phase_list,
                         solute_set,
                         doc="Solute flux as function of solvent flux")
        def solute_solvent_flux_eq(b, t, x, p, j):
            return (b.flux_mol_phase_comp[t, x, p, j] ==
                    b.flux_vol_water[t, x] * b.permeate_side[t, x].conc_mol_phase_comp[p, j])

        # 7. Extended Nernst Planck equation, DOF= Nj * 2 for inlet/outlet
        @self.Constraint(self.flowsheet().config.time,
                         io_list,
                         phase_list,
                         solute_set,
                         doc="Solute flux within pore domain")
        def solute_flux_pore_domain_eq(b, t, x, p, j):
            return (b.flux_mol_phase_comp[t, x, p, j] ==
                    - b.diffus_pore_comp[t, j]
                    * (b.pore_exit[t, x].conc_mol_phase_comp[p, j] - b.pore_entrance[t, x].conc_mol_phase_comp[p, j])
                    / b.membrane_thickness_effective
                    - 0.5 * b.feed_side.properties_in[t].charge_comp[j]
                    * (b.pore_entrance[t, x].conc_mol_phase_comp[p, j] + b.pore_exit[t, x].conc_mol_phase_comp[p, j])
                    * b.diffus_pore_comp[t, j]
                    * Constants.faraday_constant / (Constants.gas_constant * b.feed_side.properties_in[t].temperature)
                    * (b.electric_potential[t, x, 'pore_exit'] - b.electric_potential[t, x, 'pore_entrance'])
                    / b.membrane_thickness_effective
                    + 0.5 * b.hindrance_factor_convective_comp[t, j]
                    * (b.pore_entrance[t, x].conc_mol_phase_comp[p, j] + b.pore_exit[t, x].conc_mol_phase_comp[p, j])
                    * b.flux_vol_water[t, x])

        # 8. Feed-solution/membrane mass transfer resistance, DOF= Nj * 2 for inlet/outlet
        @self.Constraint(self.flowsheet().config.time,
                         io_list,
                         phase_list,
                         solute_set,
                         doc='Feed-interface mass transfer resistance accounting for concentration polarization')
        def solute_flux_concentration_polarization_eq(b, t, x, p, j):
            if x == 0:
                bulk = b.feed_side.properties_in[t]
            elif x:
                bulk = b.feed_side.properties_out[t]
            interface = b.feed_side.properties_interface[t, x]
            return (b.flux_mol_phase_comp[t, x, p, j] ==
                    - b.Kf_comp[t, x, j]
                    * (interface.conc_mol_phase_comp[p, j]
                       - bulk.conc_mol_phase_comp[p, j])
                    + b.flux_vol_water[t, x]
                    * interface.conc_mol_phase_comp[p, j]
                    - interface.charge_comp[j]
                    * interface.conc_mol_phase_comp[p, j]
                    * interface.diffus_phase_comp[p, j]
                    * Constants.faraday_constant
                    / Constants.gas_constant
                    / interface.temperature
                    * b.electric_potential_grad_feed_interface[t, x])

        # 9. Isothermal conditions at permeate inlet/outlet, DOF= 1*2 for inlet/outlet
        @self.Constraint(self.flowsheet().config.time,
                         io_list,
                         doc="Isothermal assumption for permeate")
        def eq_permeate_isothermal(b, t, x):
            return b.feed_side.properties_in[t].temperature == \
                   b.permeate_side[t, x].temperature

        # 10. Isothermal conditions at feed/membrane interface, DOF= 1*2 for inlet/outlet
        @self.feed_side.Constraint(self.flowsheet().config.time,
                                   io_list,
                                   doc="Isothermal assumption for feed-membrane interface")
        def eq_feed_interface_isothermal(b, t, x):
            return b.properties_in[t].temperature == \
                   b.properties_interface[t, x].temperature

        # Todo: since mixed permeate temp variable is unused, this constraint is unnecessary- confirm+remove
        # @self.Constraint(self.flowsheet().config.time,
        #                  doc="Isothermal assumption for mixed permeate")
        # def eq_mixed_permeate_isothermal(b, t):
        #     return b.feed_side.properties_in[t].temperature == \
        #            b.mixed_permeate[t].temperature

        # 11. Isobaric conditions at permeate side, DOF= 1*2 for inlet/outlet
        # TOdo: mixed permeate pressure is currently unused variable, but fixing its value satisfies this constraint
        @self.Constraint(self.flowsheet().config.time,
                         io_list,
                         doc="Isobaric permeate pressure")
        def eq_pressure_permeate_io(b, t, x):
            return b.permeate_side[t, x].pressure == b.mixed_permeate[t].pressure

        # 12. Javg * area = -dMf
        @self.Constraint(self.flowsheet().config.time,
                         phase_list,
                         solvent_solute_set,
                         doc="Component mass transfer from feed")
        def eq_mass_transfer_feed(b, t, p, j):
            if b.feed_side.properties_in[0].get_material_flow_basis() == MaterialFlowBasis.molar:
                return b.flux_mol_phase_comp_avg[t, p, j] * b.area == -b.feed_side.mass_transfer_term[t, p, j]

        # 13. Mass transfer equal to permeate flow terms; mole_flow,perm final = -dMf = Javg * area
        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.phase_list,
                         solvent_solute_set,
                         doc="Permeate production/average mass transfer constraint")
        def eq_permeate_production(b, t, p, j):
            return (b.mixed_permeate[t].get_material_flow_terms(p, j)
                    == b.flux_mol_phase_comp_avg[t, p, j] * b.area)

        # 14. Mole fraction at permeate inlet/outlet, DOF= Nj* 2 for inlet/out
        # @self.Constraint(self.flowsheet().config.time,
        #                  io_list,
        #                  solute_set,
        #                  doc="Permeate mass fraction at inlet/outlet")
        # def eq_mole_frac_permeate(b, t, x, j):
        #     return (b.permeate_side[t, x].mole_frac_phase_comp['Liq', j]
        #             * sum(b.flux_mol_phase_comp[t, x, 'Liq', jj]
        #                   for jj in solvent_solute_set)
        #             == b.flux_mol_phase_comp[t, x, 'Liq', j])

        # 14. Experimental constraint: Electroneutrality of final permeate
        @self.Constraint(self.flowsheet().config.time,
                         phase_list,
                         doc="Electroneutrality in mixed permeate")
        def electroneutrality_mixed_permeate_eq(b, t, p):
            return (abs(sum(b.mixed_permeate[t].conc_mol_phase_comp[p, j] *
                            b.mixed_permeate[t].charge_comp[j] for j in solute_set)) == b.tol_electroneutrality)

        # Experimental constraint: feed electroneutrality
        @self.Constraint(self.flowsheet().config.time,
                         io_list,
                         phase_list,
                         doc="Electroneutrality at feed")
        def electroneutrality_feed_eq(b, t, x, p):
            if x == 0:
                prop = b.feed_side.properties_in[t]
            elif x:
                prop = b.feed_side.properties_out[t]
            return (abs(sum(prop.conc_mol_phase_comp[p, j] *
                            prop.charge_comp[j] for j in solute_set))
                    == b.tol_electroneutrality)

        # Experimental constraint
        @self.Constraint(self.flowsheet().config.time,
                         io_list,
                         doc="Volumetric flow at interface of inlet and outlet")
        def eq_equal_flow_vol_interface(b, t, x):
            if x == 0:
                bulk = b.feed_side.properties_in[t]
            elif x:
                bulk = b.feed_side.properties_out[t]
            interface = b.feed_side.properties_interface[t, x]
            return interface.flow_vol_phase['Liq'] ==\
                   bulk.flow_vol_phase['Liq']

        # Experimental constraint
        @self.Constraint(self.flowsheet().config.time,
                         io_list,
                         doc="Volumetric flow at pore exit and permeate of inlet and outlet")
        def eq_equal_flow_vol_pore_permeate(b, t, x):
            pore = b.pore_exit[t, x]
            permeate = b.permeate_side[t, x]
            return pore.flow_vol_phase['Liq'] ==\
                   permeate.flow_vol_phase['Liq']

        # 15. Mole component recovery rate
        @self.Constraint(self.flowsheet().config.time,
                         solvent_solute_set)
        def eq_recovery_mol_phase_comp(b, t, j):
            return (b.recovery_mol_phase_comp[t, 'Liq', j] ==
                    b.mixed_permeate[t].flow_mol_phase_comp['Liq', j] /
                    b.feed_side.properties_in[t].flow_mol_phase_comp['Liq', j])

        # 16. Volumetric recovery rate
        @self.Constraint(self.flowsheet().config.time,
                         phase_list)
        def eq_recovery_vol_phase(b, t, p):
            return (b.recovery_vol_phase[t, p]
                    * b.feed_side.properties_in[t].flow_vol_phase[p]
                    == b.mixed_permeate[t].flow_vol_phase[p])

        # 17. Observed rejection rate
        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.phase_list,
                         solute_set,
                         doc="Observed solute rejection")
        def eq_rejection_phase_comp(b, t, p, j):
            return (b.mixed_permeate[t].conc_mol_phase_comp['Liq', j]
                    == b.feed_side.properties_in[t].conc_mol_phase_comp['Liq', j]
                    * (1 - b.rejection_phase_comp[t, p, j]))

        # TODO: seems stale since temperature unused at pore entrance/exit- confirm+remove;
        #  1/17/22: after including temp variables for pore in interfacial equilib eqns, this is relevant
        @self.Constraint(self.flowsheet().config.time,
                         io_list,
                         ['pore_entrance', 'pore_exit'],
                         doc="Isothermal assumption for pore inlet/outlet")
        def eq_pore_isothermal(b, t, x, y):
            if y == 'pore_entrance':
                prop = b.pore_entrance[t, x]
            elif y == 'pore_exit':
                prop = b.pore_exit[t, x]
            return b.feed_side.properties_in[t].temperature == \
                   prop.temperature


    def initialize(
            blk,
            state_args=None,
            outlvl=idaeslog.NOTSET,
            solver=None,
            optarg=None):

        """
        General wrapper for pressure changer initialization routines

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None)
            solver : str indicating which solver to use during
                     initialization (default = None)

        Returns: None
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")
        # Set solver options
        if optarg is None:
            optarg = {'bound_push': 1e-8}

        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Initialize holdup block
        flags = blk.feed_side.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )
        init_log.info_high("Initialization Step 1 Complete.")
        # ---------------------------------------------------------------------
        # Initialize permeate
        # Set state_args from inlet state
        if state_args is None:
            state_args = {}
            state_dict = blk.feed_side.properties_in[
                blk.flowsheet().config.time.first()].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        blk.mixed_permeate.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )
        init_log.info_high("Initialization Step 2 Complete.")

        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high(
            "Initialization Step 3 {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.feed_side.release_state(flags, outlvl + 1)
        init_log.info(
            "Initialization Complete: {}".format(idaeslog.condition(res))
        )

    def _get_performance_contents(self, time_point=0):
        pass
        # for k in ('ion_set', 'solute_set'):
        #     if hasattr(self.config.property_package, k):
        #         solute_set = getattr(self.config.property_package, k)
        #         break
        # var_dict = {}
        # expr_dict = {}
        # var_dict["Volumetric Recovery Rate"] = self.recovery_vol_phase[time_point, 'Liq']
        # var_dict["Solvent Mass Recovery Rate"] = self.recovery_mass_phase_comp[time_point, 'Liq', 'H2O']
        # var_dict["Membrane Area"] = self.area
        # if hasattr(self, "deltaP"):
        #     var_dict["Pressure Change"] = self.deltaP[time_point]
        # if self.feed_side.properties_in[time_point].is_property_constructed('flow_vol'):
        #     if self.feed_side.properties_in[time_point].flow_vol.is_variable_type():
        #         obj_dict = var_dict
        #     elif self.feed_side.properties_in[time_point].flow_vol.is_named_expression_type():
        #         obj_dict = expr_dict
        #     else:
        #         raise Exception(f"{self.feed_side.properties_in[time_point].flow_vol} isn't a variable nor expression")
        #     obj_dict['Volumetric Flowrate @Inlet'] = self.feed_side.properties_in[time_point].flow_vol
        # if self.feed_side.properties_out[time_point].is_property_constructed('flow_vol'):
        #     if self.feed_side.properties_out[time_point].flow_vol.is_variable_type():
        #         obj_dict = var_dict
        #     elif self.feed_side.properties_out[time_point].flow_vol.is_named_expression_type():
        #         obj_dict = expr_dict
        #     else:
        #         raise Exception(f"{self.feed_side.properties_in[time_point].flow_vol} isn't a variable nor expression")
        #     obj_dict['Volumetric Flowrate @Outlet'] = self.feed_side.properties_out[time_point].flow_vol
        # var_dict['Solvent Volumetric Flux']= self.flux_vol_solvent[time_point, 'H2O']
        # for j in solute_set:
        #     var_dict[f'{j} Rejection'] = self.rejection_phase_comp[time_point, 'Liq', j]
        #     if self.feed_side.properties_in[time_point].conc_mol_phase_comp['Liq', j].is_expression_type():
        #         obj_dict = expr_dict
        #     elif self.feed_side.properties_in[time_point].conc_mol_phase_comp['Liq', j].is_variable_type():
        #         obj_dict = var_dict
        #     obj_dict[f'{j} Molar Concentration @Inlet'] = self.feed_side.properties_in[time_point].conc_mol_phase_comp['Liq', j]
        #     obj_dict[f'{j} Molar Concentration @Outlet'] = self.feed_side.properties_out[time_point].conc_mol_phase_comp['Liq', j]
        #     obj_dict[f'{j} Molar Concentration @Permeate'] = self.properties_permeate[time_point].conc_mol_phase_comp['Liq', j]
        #
        # return {"vars": var_dict, "exprs": expr_dict}

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Feed Inlet": self.inlet,
                "Feed Outlet": self.retentate,
                "Permeate Outlet": self.permeate,
            },
            time_point=time_point,
        )

    def get_costing(self, module=None, **kwargs):
        self.costing = Block()
        module.Nanofiltration_costing(self.costing, **kwargs)

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        for k in ('ion_set', 'solute_set'):
            if hasattr(self.config.property_package, k):
                solute_set = getattr(self.config.property_package, k)
                break
        if iscale.get_scaling_factor(self.radius_pore) is None:
            sf = iscale.get_scaling_factor(self.radius_pore, default=1e10, warning=True)
            iscale.set_scaling_factor(self.radius_pore, sf)

        # # TODO: require users to set scaling factor for area or calculate it based on mass transfer and flux
        # iscale.set_scaling_factor(self.area, 1e-1)
        #
        # # setting scaling factors for variables
        # # these variables should have user input, if not there will be a warning
        # if iscale.get_scaling_factor(self.area) is None:
        #     sf = iscale.get_scaling_factor(self.area, default=1, warning=True)
        #     iscale.set_scaling_factor(self.area, sf)
        #
        # # these variables do not typically require user input,
        # # will not override if the user does provide the scaling factor
        # # TODO: this default scaling assumes SI units rather than being based on the property package
        # if iscale.get_scaling_factor(self.dens_solvent) is None:
        #     iscale.set_scaling_factor(self.dens_solvent, 1e-3)
        #
        # for t, v in self.flux_vol_solvent.items():
        #     if iscale.get_scaling_factor(v) is None:
        #         iscale.set_scaling_factor(v, 1e6)
        #
        # for (t, p, j), v in self.rejection_phase_comp.items():
        #     if iscale.get_scaling_factor(v) is None:
        #         iscale.set_scaling_factor(v, 1e1)
        #
        # for (t, p, j), v in self.mass_transfer_phase_comp.items():
        #     if iscale.get_scaling_factor(v) is None:
        #         sf = 10 * iscale.get_scaling_factor(self.feed_side.properties_in[t].get_material_flow_terms(p, j))
        #         iscale.set_scaling_factor(v, sf)
        # if iscale.get_scaling_factor(self.recovery_vol_phase) is None:
        #     iscale.set_scaling_factor(self.recovery_vol_phase, 1)
        #
        # for (t, p, j), v in self.recovery_mass_phase_comp.items():
        #     if j in self.config.property_package.solvent_set:
        #         sf = 1
        #     elif j in solute_set:
        #         sf = 10
        #     if iscale.get_scaling_factor(v) is None:
        #         iscale.set_scaling_factor(v, sf)
        #
        # # transforming constraints
        # for ind, c in self.feed_side.eq_isothermal.items():
        #     sf = iscale.get_scaling_factor(self.feed_side.properties_in[0].temperature)
        #     iscale.constraint_scaling_transform(c, sf)
        #
        # for ind, c in self.eq_mass_transfer_term.items():
        #     sf = iscale.get_scaling_factor(self.mass_transfer_phase_comp[ind])
        #     iscale.constraint_scaling_transform(c, sf)
        #
        # for ind, c in self.eq_solvent_transfer.items():
        #     sf = iscale.get_scaling_factor(self.mass_transfer_phase_comp[ind])
        #     iscale.constraint_scaling_transform(c, sf)
        #
        # for ind, c in self.eq_permeate_production.items():
        #     sf = iscale.get_scaling_factor(self.mass_transfer_phase_comp[ind])
        #     iscale.constraint_scaling_transform(c, sf)
        #
        # for ind, c in self.eq_rejection_phase_comp.items():
        #     sf = iscale.get_scaling_factor(self.rejection_phase_comp[ind])
        #     iscale.constraint_scaling_transform(c, sf)
        #
        # for t, c in self.eq_permeate_isothermal.items():
        #     sf = iscale.get_scaling_factor(self.feed_side.properties_in[t].temperature)
        #     iscale.constraint_scaling_transform(c, sf)
        # for t, c in self.eq_recovery_vol_phase.items():
        #     sf = iscale.get_scaling_factor(self.recovery_vol_phase[t, 'Liq'])
        #     iscale.constraint_scaling_transform(c, sf)
        #
        # for (t, j), c in self.eq_recovery_mass_phase_comp.items():
        #     sf = iscale.get_scaling_factor(self.recovery_mass_phase_comp[t, 'Liq', j])
        #     iscale.constraint_scaling_transform(c, sf)
