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
from pyomo.environ import (
    Block,
    Set,
    Var,
    Param,
    Suffix,
    NonNegativeReals,
    Reference,
    Constraint,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog


_log = idaeslog.getLogger(__name__)

# when using this file the name "Filtration" is what is imported
@declare_process_block_class("Filtration")
class FiltrationData(UnitModelBlockData):
    """
    Zero order filtration model
    """

    # CONFIG are options for the unit model, this simple model only has the mandatory config options
    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False. The filtration unit does not support dynamic
    behavior, thus this must be False.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False. The filtration unit does not have defined volume, thus
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
    and used when constructing these,
    **default** - None.
    **Valid values:** {
    see property package for documentation.}""",
        ),
    )

    def build(self):
        # build always starts by calling super().build()
        # This triggers a lot of boilerplate in the background for you
        super().build()

        # this creates blank scaling factors, which are populated later
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # Next, get the base units of measurement from the property definition
        units_meta = self.config.property_package.get_metadata().get_derived_units

        # Add unit variables
        self.recovery_mass_phase_comp = Var(
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=0.5,
            bounds=(1e-8, 1),
            units=pyunits.dimensionless,
        )
        self.removal_fraction_mass_phase_comp = Var(
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=0.5,
            bounds=(1e-8, 1),
            units=pyunits.dimensionless,
        )

        # Add state blocks for inlet, outlet, and waste
        # These include the state variables and any other properties on demand
        # Add inlet block
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True  # inlet block is an inlet
        self.properties_in = self.config.property_package.state_block_class(
            self.flowsheet().config.time, doc="Material properties of inlet", **tmp_dict
        )
        # Add outlet and waste block
        tmp_dict["defined_state"] = False  # outlet and waste block is not an inlet
        self.properties_out = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of outlet",
            **tmp_dict
        )
        self.properties_waste = self.config.property_package.state_block_class(
            self.flowsheet().config.time, doc="Material properties of waste", **tmp_dict
        )

        # Add ports - oftentimes users interact with these rather than the state blocks
        self.add_port(name="inlet", block=self.properties_in)
        self.add_port(name="outlet", block=self.properties_out)
        self.add_port(name="waste", block=self.properties_waste)

        # Add constraints
        # Usually unit models use a control volume to do the mass, energy, and momentum
        # balances, however, in this example I'm just going to be explicit and write them myself
        @self.Constraint(
            self.config.property_package.component_list,
            doc="Mass balance for components",
        )
        def eq_mass_balance(b, j):
            return (
                b.properties_in[0].flow_mass_phase_comp["Liq", j]
                == b.properties_out[0].flow_mass_phase_comp["Liq", j]
                + b.properties_waste[0].flow_mass_phase_comp["Liq", j]
            )

        @self.Constraint(doc="Isothermal assumption 1")
        def eq_isothermal_1(b):
            return b.properties_in[0].temperature == b.properties_out[0].temperature

        @self.Constraint(doc="Isothermal assumption 2")
        def eq_isothermal_2(b):
            return b.properties_in[0].temperature == b.properties_waste[0].temperature

        @self.Constraint(doc="Isobaric assumption")
        def eq_isobaric(b):
            return b.properties_in[0].pressure == b.properties_out[0].pressure

        # waste pressure should be set by the user

        # Next, add constraints linking variables
        @self.Constraint(
            self.config.property_package.component_list,
            doc="Relationship between removal and recovery",
        )
        def eq_removal_recovery(b, j):
            return (
                b.removal_fraction_mass_phase_comp["Liq", j]
                == 1 - b.recovery_mass_phase_comp["Liq", j]
            )

        @self.Constraint(
            self.config.property_package.component_list,
            doc="Component removal equation",
        )
        def eq_removal_balance(b, j):
            return (
                b.properties_in[0].flow_mass_phase_comp["Liq", j]
                * b.removal_fraction_mass_phase_comp["Liq", j]
                == b.properties_waste[0].flow_mass_phase_comp["Liq", j]
            )

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        iscale.set_scaling_factor(self.recovery_mass_phase_comp, 1)
        iscale.set_scaling_factor(self.removal_fraction_mass_phase_comp, 1)

        # transforming constraints
        for j, c in self.eq_mass_balance.items():
            sf = iscale.get_scaling_factor(
                self.properties_in[0].flow_mass_phase_comp["Liq", j]
            )
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_isothermal_1.items():
            sf = iscale.get_scaling_factor(self.properties_in[0].temperature)
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_isothermal_2.items():
            sf = iscale.get_scaling_factor(self.properties_in[0].temperature)
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_isobaric.items():
            sf = iscale.get_scaling_factor(self.properties_in[0].pressure)
            iscale.constraint_scaling_transform(c, sf)

        # other constraints don't need to be transformed
