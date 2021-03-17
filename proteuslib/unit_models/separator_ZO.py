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
from pyomo.environ import Block, Constraint, Var, Suffix, NonNegativeReals, Reals, \
    SolverFactory, units as pyunits
from pyomo.network import Port

# Import IDAES cores
import idaes.logger as idaeslog
from idaes.core import (declare_process_block_class,
                        UnitModelBlockData,
                        useDefault)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.testing import get_default_solver
import idaes.core.util.scaling as iscale

_log = idaeslog.getLogger(__name__)

@declare_process_block_class("SeparatorZO")
class SeparatorZOData(UnitModelBlockData):
    """
    Standard Zero Order Separator Unit Model Class:
    - steady state only
    """
    CONFIG = ConfigBlock()

    CONFIG.declare("dynamic", ConfigValue(
        domain=In([False]),
        default=False,
        description="Dynamic model flag - must be False",
        doc="""Indicates whether this model will be dynamic or not,
    **default** = False. Zero Order Separators do not support dynamic behavior."""))
    CONFIG.declare("has_holdup", ConfigValue(
        default=False,
        domain=In([False]),
        description="Holdup construction flag - must be False",
        doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False. Zero Order Separators do not have defined volume, thus
    this must be False."""))
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

    def build(self):
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        units_meta = self.config.property_package.get_metadata().get_derived_units

        # Add state blocks for inlet, outlet, and waste
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

        # Add ports
        self.add_inlet_port(name='inlet', block=self.inlet_state)
        self.add_outlet_port(name='outlet', block=self.outlet_state)
        self.add_port(name='waste', block=self.waste_state)

        # Add additional variables
        self.recovery_frac_phase_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=0.8,
            bounds=(0, 1),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Water recovery fraction")
        self.removal_frac_phase_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=0.75,
            bounds=(0, 1),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Component removal fraction")

        if (self.config.has_pressure_change
                and self.outlet_state[0].is_property_constructed('pressure')):
            self.deltaP_outlet = Var(
                self.flowsheet().config.time,
                initialize=1e4,
                bounds=(1e-8, 1e8),
                domain=Reals,
                units=units_meta("pressure"),
                doc="Pressure change between inlet and outlet")
            self.deltaP_waste = Var(
                self.flowsheet().config.time,
                initialize=1e4,
                bounds=(1e-8, 1e8),
                domain=Reals,
                units=units_meta("pressure"),
                doc="Pressure change between inlet and waste")

        # Add performance equations
        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Component mass balance")
        def eq_component_mass_balance(b, t, p, j):
            return (b.inlet_state[t].get_material_flow_terms(p, j) ==
                    b.outlet_state[t].get_material_flow_terms(p, j)
                    + b.waste_state[t].get_material_flow_terms(p, j))

        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Component removal")
        def eq_component_removal(b, t, p, j):
            return (b.removal_frac_phase_comp[t, p, j]
                    * b.inlet_state[t].get_material_flow_terms(p, j) ==
                    b.waste_state[t].get_material_flow_terms(p, j))

        @self.Constraint(self.flowsheet().config.time,
                         self.config.property_package.phase_list,
                         self.config.property_package.component_list,
                         doc="Removal fraction to recovery fraction")
        def eq_removal_to_recovery(b, t, p, j):
            return (b.removal_frac_phase_comp[t, p, j] ==
                    1 - b.recovery_frac_phase_comp[t, p, j])

        if self.inlet_state[0].is_property_constructed('temperature'):
            @self.Constraint(self.flowsheet().config.time,
                             doc="Isothermal outlet temperature")
            def eq_outlet_temperature(b, t):
                return b.outlet_state[t].temperature == b.inlet_state[t].temperature

            @self.Constraint(self.flowsheet().config.time,
                             doc="Isothermal waste temperature")
            def eq_waste_temperature(b, t):
                return b.waste_state[t].temperature == b.inlet_state[t].temperature

        if self.inlet_state[0].is_property_constructed('pressure'):
            if self.config.has_pressure_change:
                @self.Constraint(self.flowsheet().config.time,
                                 doc="Outlet pressure equation")
                def eq_outlet_pressure(b, t):
                    return (b.inlet_state[t].pressure ==
                            b.outlet_state[t].pressure + b.deltaP_outlet[t])

                @self.Constraint(self.flowsheet().config.time,
                                 doc="Waste pressure equation")
                def eq_waste_pressure(b, t):
                    return (b.inlet_state[t].pressure ==
                            b.waste_state[t].pressure + b.deltaP_waste[t])
            else:
                @self.Constraint(self.flowsheet().config.time,
                                 doc="Isobaric outlet pressure")
                def eq_outlet_pressure(b, t):
                    return b.inlet_state[t].pressure == b.outlet_state[t].pressure

                @self.Constraint(self.flowsheet().config.time,
                                 doc="Isobaric waste pressure")
                def eq_waste_pressure(b, t):
                    return b.inlet_state[t].pressure == b.waste_state[t].pressure

    def initialize(
            self,
            state_args=None,
            routine=None,
            outlvl=idaeslog.NOTSET,
            solver=None,
            optarg={"tol": 1e-6}):
        """
        General wrapper for pressure changer initialization routines
        Keyword Arguments:
            routine : str stating which initialization routine to execute
                        * None - currently no specialized routine for RO unit
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default={'tol': 1e-6})
            solver : str indicating whcih solver to use during
                     initialization (default = 'ipopt')
        Returns:
            None
        """
        # Get loggers
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="properties")

        # Set solver and options
        # TODO: clean up once IDAES new API for initialize solvers is released
        if isinstance(solver, str):
            opt = SolverFactory(solver)
            opt.options = optarg
        else:
            if solver is None:
                opt = get_default_solver()
            else:
                opt = solver
                opt.options = optarg

        # Solve simple unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high(
            "Initialized zero order separator.".format(idaeslog.condition(res)))

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # scale variables
        for (t, p, j), v in self.recovery_frac_phase_comp.items():
            recovery_var = v
            removal_var = self.removal_frac_phase_comp[t, p, j]
            if (iscale.get_scaling_factor(recovery_var) is None
                    and iscale.get_scaling_factor(removal_var) is None):
                sf = iscale.get_scaling_factor(removal_var, default=1, warning=True)
                iscale.set_scaling_factor(recovery_var, sf)
                iscale.set_scaling_factor(removal_var, sf)
            elif iscale.get_scaling_factor(recovery_var) is None:
                sf = 1 / iscale.get_scaling_factor(removal_var)
                iscale.set_scaling_factor(recovery_var, sf)
            elif iscale.get_scaling_factor(removal_var) is None:
                sf = 1 / iscale.get_scaling_factor(recovery_var)
                iscale.set_scaling_factor(removal_var, sf)

        for t, v in self.deltaP_outlet.items():
            if iscale.get_scaling_factor(v) is None:
                sf = iscale.get_scaling_factor(self.inlet_state[t].pressure)
                sf = sf / 10
                iscale.set_scaling_factor(v, sf)

        for t, v in self.deltaP_waste.items():
            if iscale.get_scaling_factor(v) is None:
                sf = iscale.get_scaling_factor(self.inlet_state[t].pressure)
                sf = sf / 10
                iscale.set_scaling_factor(v, sf)

        # transform constraints
        for (t, p, j), c in self.eq_component_mass_balance.items():
            sf = iscale.get_scaling_factor(self.inlet_state[t].get_material_flow_terms(p, j))
            iscale.constraint_scaling_transform(c, sf)

        for (t, p, j), c in self.eq_component_removal.items():
            sf = iscale.get_scaling_factor(self.inlet_state[t].get_material_flow_terms(p, j))
            iscale.constraint_scaling_transform(c, sf)

        for (t, p, j), c in self.eq_removal_to_recovery.items():
            sf = iscale.get_scaling_factor(self.removal_frac_phase_comp[t, p, j])
            iscale.constraint_scaling_transform(c, sf)

        for t, c in self.eq_outlet_temperature.items():
            sf = iscale.get_scaling_factor(self.inlet_state[t].temperature)
            iscale.constraint_scaling_transform(c, sf)

        for t, c in self.eq_waste_temperature.items():
            sf = iscale.get_scaling_factor(self.inlet_state[t].temperature)
            iscale.constraint_scaling_transform(c, sf)

        for t, c in self.eq_outlet_pressure.items():
            sf = iscale.get_scaling_factor(self.inlet_state[t].pressure)
            iscale.constraint_scaling_transform(c, sf)

        for t, c in self.eq_waste_pressure.items():
            sf = iscale.get_scaling_factor(self.inlet_state[t].pressure)
            iscale.constraint_scaling_transform(c, sf)
