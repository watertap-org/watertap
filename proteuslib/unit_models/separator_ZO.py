##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################

# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.environ import Block, Constraint, Var, units as pyunits
from pyomo.network import Port

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        UnitModelBlockData,
                        useDefault)
from idaes.core.util.config import is_physical_parameter_block

@declare_process_block_class("ZeroOrder")
class ZeroOrderData(UnitModelBlockData):
    """
    This class describes the rules for a zero-order model
    """
    # The Config Block is used tpo process arguments from when the model is
    # instantiated. In IDAES, this serves two purposes:
    #     1. Allows us to separate physical properties from unit models
    #     2. Lets us give users options for configuring complex units
    # For WaterTAP3, this will mainly be boilerplate to keep things consistent
    # with ProteusLib and IDAES.
    # The dynamic and has_holdup options are expected arguments which must exist
    # The property package arguments let us define different sets of contaminants
    # without needing to write a new model.
    CONFIG = ConfigBlock()
    CONFIG.declare("dynamic", ConfigValue(
        domain=In([False]),
        default=False,
        description="Dynamic model flag - must be False",
        doc="""Indicates whether this model will be dynamic or not,
**default** = False. Equilibrium Reactors do not support dynamic behavior."""))
    CONFIG.declare("has_holdup", ConfigValue(
        default=False,
        domain=In([False]),
        description="Holdup construction flag - must be False",
        doc="""Indicates whether holdup terms should be constructed or not.
**default** - False. Equilibrium reactors do not have defined volume, thus
this must be False."""))
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

    def build(self):
        """
        The build method is the core of the unit model, and contains the rules
        for building the Vars and Constraints that make up the unit model.
        """
        # build always starts by calling super().build()
        # This triggers a lot of boilerplate in the background for you
        super(ZeroOrderData, self).build()

        # Next, get the base units of measurement from the property definition
        units_meta = self.config.property_package.get_metadata().get_derived_units

        # Also need to get time domain
        # This will not be used for WaterTAP3, but will be needed to integrate
        # with ProteusLib dynamic models
        time = self.flowsheet().config.time

        # Add state blocks for inlet, outlet, and waste
        # These include the state variables and any other properties on demand
        # Add inlet block
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True  # inlet block is an inlet
        self.inlet_state = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of inlet",
            default=tmp_dict)
        # Add outlet and waste block
        tmp_dict["defined_state"] = False  # outlet and waste block is not an inlet
        self.outlet_state = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of outlet",
            default=tmp_dict)
        self.waste_state = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of waste",
            default=tmp_dict)

        # Next, add additional variables for unit performance
        self.water_recovery = Var(time,
                                  initialize=0.8,
                                  bounds=(0, 1),
                                  units=pyunits.dimensionless,
                                  doc="Water recovery fraction")
        self.removal_fraction = Var(time,
                                    self.config.property_package.component_list,
                                    initialize=0.75,
                                    bounds=(0, 1),
                                    units=pyunits.dimensionless,
                                    doc="Component removal fraction")
        self.deltaP_outlet = Var(time,
                                 initialize=1e4,
                                 units=units_meta("pressure"),
                                 doc="Pressure change between inlet and outlet")
        self.deltaP_waste = Var(time,
                                initialize=1e4,
                                units=units_meta("pressure"),
                                doc="Pressure change between inlet and waste")

        # Next, add constraints linking these variables
        @self.Constraint(time, doc="Water recovery to removal fraction")
        def recovery_to_removal(b, t):
            return b.removal_fraction[t, 'H2O'] == 1 - b.water_recovery[t]

        @self.Constraint(time,
                         self.config.property_package.component_list,
                         doc="Component mass balance")
        def component_mass_balance(b, t, j):
            return (b.inlet_state[t].flow_mass_comp[j] ==
                    b.outlet_state[t].flow_mass_comp[j] + b.waste_state[t].flow_mass_comp[j])

        if self.inlet_state[0].is_property_constructed('temperature'):
            @self.Constraint(time, doc="Outlet temperature equation")
            def outlet_temperature_constraint(b, t):
                return b.outlet_state[t].temperature == b.inlet_state[t].temperature

            @self.Constraint(time, doc="Waste temperature equation")
            def waste_temperature_constraint(b, t):
                return b.waste_state[t].temperature == b.inlet_state[t].temperature

        if self.outlet_state[0].is_property_constructed('pressure'):
            @self.Constraint(time, doc="Outlet pressure equation")
            def outlet_pressure_constraint(b, t):
                return (b.inlet_state[t].pressure ==
                        b.outlet_state[t].pressure + b.deltaP_outlet[t])

            @self.Constraint(time, doc="Waste pressure equation")
            def waste_pressure_constraint(b, t):
                return (b.inlet_state[t].pressure ==
                        b.waste_state[t].pressure + b.deltaP_waste[t])

        # Finally, add removal (and recovery) equations
        @self.Constraint(time,
                         self.config.property_package.component_list,
                         doc="Component removal equation")
        def component_removal_equation(b, t, j):
            return (b.removal_fraction[t, j] * b.inlet_state[t].flow_mass_comp[j] ==
                    b.waste_state[t].flow_mass_comp[j])

        # The last step is to create Ports representing the three streams
        self.add_inlet_port(name='inlet', block=self.inlet_state)
        self.add_outlet_port(name='outlet', block=self.outlet_state)
        self.add_port(name='waste', block=self.waste_state)

    # def initialization(self, *args, **kwargs):
    #     # All IDAES models are expected ot have an initialization routine
    #     # We will need to add one here and it will be fairly simple,
    #     # but I will skip it for now
    #     pass
    #
    # def get_costing(self, module=financials, cost_method="wt", year=None):
    #     """
    #     We need a get_costing method here to provide a point to call the
    #     costing methods, but we call out to an external costing module
    #     for the actual calculations. This lets us easily swap in different
    #     methods if needed.
    #
    #     Within IDAES, the year argument is used to set the initial value for
    #     the cost index when we build the model.
    #     """
    #     # First, check to see if global costing module is in place
    #     # Construct it if not present and pass year argument
    #     if not hasattr(self.flowsheet(), "costing"):
    #         self.flowsheet().get_costing(module=module, year=year)
    #
    #     # Next, add a sub-Block to the unit model to hold the cost calculations
    #     # This is to let us separate costs from model equations when solving
    #     self.costing = Block()
    #     # Then call the appropriate costing function out of the costing module
    #     # The first argument is the Block in which to build the equations
    #     # Can pass additional arguments a needed
    #     module.nf_costing(self.costing, cost_method=cost_method)
