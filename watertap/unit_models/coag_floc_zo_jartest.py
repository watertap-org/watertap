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
                           Reference,
                           units as pyunits)
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
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

__author__ = "Austin Ladshaw"

_log = idaeslog.getLogger(__name__)

# Name of the unit model
@declare_process_block_class("CoagulationFlocculationZO_JarTestModel")
class CoagulationFlocculationZO_JarTestModelData(UnitModelBlockData):
    """
    Zero order Coagulation-Flocculation model based on Jar Tests
    """
    # CONFIG are options for the unit model
    CONFIG = ConfigBlock()

    CONFIG.declare("dynamic", ConfigValue(
        domain=In([False]),
        default=False,
        description="Dynamic model flag - must be False",
        doc="""Indicates whether this model will be dynamic or not,
    **default** = False. The filtration unit does not support dynamic
    behavior, thus this must be False."""))

    CONFIG.declare("has_holdup", ConfigValue(
        default=False,
        domain=In([False]),
        description="Holdup construction flag - must be False",
        doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False. The filtration unit does not have defined volume, thus
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

    CONFIG.declare("energy_balance_type", ConfigValue(
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

    CONFIG.declare("chemical_additives", ConfigValue(
        default={},
        domain=dict,
        description="""Dictionary of chemical additives used in coagulation process,
        along with their molecular weights, the moles of salt produced per mole of
        chemical added, and the molecular weights of the salt produced by the chemical
        additive with the format of: {'chem': (mw_chem, moles_salt_added, mw_salt)} """))


    def build(self):
        # build always starts by calling super().build()
        # This triggers a lot of boilerplate in the background for you
        super().build()

        # this creates blank scaling factors, which are populated later
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # Next, get the base units of measurement from the property definition
        units_meta = self.config.property_package.get_metadata().get_derived_units

        # check the optional config arg 'chemical_additives'
        common_msg = "The 'chemical_additives' dict MUST contain a tuple of 3 numbers for each key " + \
                     "in the following format: " + \
                     "\n {'chem': (mw_chem, moles_salt_added_per_mole_chem, mw_salt)}.\n"
        for j in self.config.chemical_additives:
            if type(self.config.chemical_additives[j]) != tuple:
                raise ConfigurationError(common_msg + "\n Did not provide a 'tuple'\n")
            if len(self.config.chemical_additives[j]) != 3:
                raise ConfigurationError(common_msg + "\n Did not provide 3 args in 'tuple'\n")
            for item in self.config.chemical_additives[j]:
                if not isinstance(item, (int,float)):
                    raise ConfigurationError(common_msg + "\n Provided arg is not a number\n")

        # Add unit variables
        # Linear relationship between TSS (mg/L) and Turbidity (NTU)
        #           TSS (mg/L) = Turbidity (NTU) * slope + intercept
        #   Default values come from the following paper:
        #       H. Rugner, M. Schwientek,B. Beckingham, B. Kuch, P. Grathwohl,
        #       Environ. Earth Sci. 69 (2013) 373-380. DOI: 10.1007/s12665-013-2307-1
        self.slope = Var(
            self.flowsheet().config.time,
            initialize=1.86,
            bounds=(1e-8, 10),
            domain=NonNegativeReals,
            units=pyunits.mg/pyunits.L,
            doc='Slope relation between TSS (mg/L) and Turbidity (NTU)')

        self.intercept = Var(
            self.flowsheet().config.time,
            initialize=0,
            bounds=(0, 10),
            domain=NonNegativeReals,
            units=pyunits.mg/pyunits.L,
            doc='Intercept relation between TSS (mg/L) and Turbidity (NTU)')

        self.initial_turbidity_ntu = Var(
            self.flowsheet().config.time,
            initialize=50,
            bounds=(0, 1000),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc='Initial measured Turbidity (NTU) from Jar Test')

        self.final_turbidity_ntu = Var(
            self.flowsheet().config.time,
            initialize=1,
            bounds=(0, 1000),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc='Final measured Turbidity (NTU) from Jar Test')

        self.chemical_doses = Var(
            self.flowsheet().config.time,
            self.config.chemical_additives.keys(),
            initialize=0,
            bounds=(0, 100),
            domain=NonNegativeReals,
            units=pyunits.mg/pyunits.L,
            doc='Dosages of the set of chemical additives')


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

        self.feed_side.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_enthalpy_transfer=False)

        self.feed_side.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=False)

        # Add ports
        self.add_inlet_port(name='inlet', block=self.feed_side)
        self.add_outlet_port(name='outlet', block=self.feed_side)

        # Check _phase_component_set for required items
        if ('Liq', 'TSS') not in self.config.property_package._phase_component_set:
            raise ConfigurationError(
                "Coagulation-Flocculation model MUST contain ('Liq','TSS') as a component, but "
                "the property package has only specified the following components {}"
                    .format([p for p in self.config.property_package._phase_component_set]))
        if ('Liq', 'TDS') not in self.config.property_package._phase_component_set:
            raise ConfigurationError(
                "Coagulation-Flocculation model MUST contain ('Liq','TDS') as a component, but "
                "the property package has only specified the following components {}"
                    .format([p for p in self.config.property_package._phase_component_set]))
        if ('Liq', 'Sludge') not in self.config.property_package._phase_component_set:
            raise ConfigurationError(
                "Coagulation-Flocculation model MUST contain ('Liq','Sludge') as a component, but "
                "the property package has only specified the following components {}"
                    .format([p for p in self.config.property_package._phase_component_set]))

        # -------- Add constraints ---------

        # Constraint for tss loss rate based on measured final turbidity
        self.tss_loss_rate = Var(
            self.flowsheet().config.time,
            initialize=1,
            bounds=(0, 100),
            domain=NonNegativeReals,
            units=units_meta('mass')*units_meta('time')**-1,
            doc='Mass per time loss rate of TSS based on the measured final turbidity')

        @self.Constraint(self.flowsheet().config.time,
                         doc="Constraint for the loss rate of TSS to be used in mass_transfer_term")
        def eq_tss_loss_rate(self, t):
            tss_out = pyunits.convert(self.slope[t]*self.final_turbidity_ntu[t] + self.intercept[t],
                                    to_units=units_meta('mass')*units_meta('length')**-3)
            input_rate = self.feed_side.properties_in[t].flow_mass_phase_comp['Liq','TSS']
            exit_rate = self.feed_side.properties_out[t].flow_vol_phase['Liq']*tss_out

            return (self.tss_loss_rate[t] == input_rate - exit_rate)

        # Add constraints for mass transfer terms
        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Mass transfer term")
        def eq_mass_transfer_term(self, t, p, j):
            if (p, j) == ('Liq', 'TSS'):
                return self.feed_side.mass_transfer_term[t, p, j] == -self.tss_loss_rate[t]
            elif (p, j) == ('Liq', 'Sludge'):
                return self.feed_side.mass_transfer_term[t, p, j] == self.tss_loss_rate[t]
            else:
                return self.feed_side.mass_transfer_term[t, p, j] == 0.0

        self.eq_tss_loss_rate.pprint()

        self.eq_mass_transfer_term.pprint()

        self.feed_side.material_balances.pprint()
        self.feed_side.enthalpy_balances.pprint()
        self.feed_side.pressure_balance.pprint()

        # enthalpy transfer
        #@self.Constraint(self.flowsheet().config.time,
        #                self.config.property_package.phase_list,
        #                 doc="Enthalpy transfer from inlet to outlet")
        #def eq_connect_enthalpy_transfer(b, t, p):
        #    return (0.0 == b.feed_side.enthalpy_transfer[t])




    # # TODO: Add user access functions for setting up the problem

    # # TODO: Add initialize method


    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # # TODO:  Fill out scaling
