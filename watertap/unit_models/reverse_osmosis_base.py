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

from copy import deepcopy
from pyomo.common.collections import ComponentSet
from pyomo.common.config import Bool, ConfigValue
from pyomo.environ import (
    NonNegativeReals,
    Param,
    Suffix,
    Var,
    check_optimal_termination,
    exp,
    units as pyunits,
)
from idaes.core import UnitModelBlockData
from watertap.core.solvers import get_solver
from idaes.core.util import scaling as iscale
from idaes.core.util.exceptions import ConfigurationError, InitializationError
from idaes.core.util.misc import add_object_reference
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.tables import create_stream_table_dataframe
import idaes.logger as idaeslog

from watertap.core import InitializationMixin
from watertap.core.membrane_channel_base import (
    validate_membrane_config_args,
    ConcentrationPolarizationType,
    TransportModel,
    ModuleType,
)
from watertap.core.util.initialization import interval_initializer
from watertap.costing.unit_models.reverse_osmosis import cost_reverse_osmosis


def _add_has_full_reporting(config_obj):
    config_obj.declare(
        "has_full_reporting",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Level of reporting results",
            doc="""Level of reporting results.
            **default** - False.
            **Valid values:** {
            **False** - include minimal reporting of results,
            **True** - report additional properties of interest that aren't constructed by
            the unit model by default. Also, report averaged expression values""",
        ),
    )


class ReverseOsmosisBaseData(InitializationMixin, UnitModelBlockData):
    """
    Reverse Osmosis base class
    """

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

        self._add_feed_side_membrane_channel_and_geometry()

        self.feed_side.add_state_blocks(has_phase_equilibrium=False)

        self.feed_side.add_material_balances(
            balance_type=self.config.material_balance_type, has_mass_transfer=True
        )

        self.feed_side.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            pressure_change_type=self.config.pressure_change_type,
            has_pressure_change=self.config.has_pressure_change,
        )

        self.feed_side.add_control_volume_isothermal_conditions()
        self.feed_side.add_interface_isothermal_conditions()

        self.feed_side.add_extensive_flow_to_interface()

        self.feed_side.add_concentration_polarization(
            concentration_polarization_type=self.config.concentration_polarization_type,
            mass_transfer_coefficient=self.config.mass_transfer_coefficient,
        )

        self.feed_side.apply_transformation()

        self.feed_side.add_expressions()

        add_object_reference(self, "length_domain", self.feed_side.length_domain)
        add_object_reference(
            self, "difference_elements", self.feed_side.difference_elements
        )
        add_object_reference(self, "first_element", self.feed_side.first_element)
        add_object_reference(self, "nfe", self.feed_side.nfe)

        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["defined_state"] = False  # these blocks are not inlets or outlets
        self.permeate_side = self.config.property_package.build_state_block(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Material properties of permeate along permeate channel",
            **tmp_dict,
        )
        self.mixed_permeate = self.config.property_package.build_state_block(
            self.flowsheet().config.time,
            doc="Material properties of mixed permeate exiting the module",
            **tmp_dict,
        )

        # Add Ports
        self.add_inlet_port(name="inlet", block=self.feed_side)
        self.add_outlet_port(name="retentate", block=self.feed_side)
        self.add_port(name="permeate", block=self.mixed_permeate)

        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Isobaric assumption for permeate out",
        )

        # Add other equations
        def eq_permeate_outlet_isobaric(b, t, x):
            return b.permeate_side[t, x].pressure == b.mixed_permeate[t].pressure

        # # ==========================================================================
        # Feed and permeate-side isothermal conditions
        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Isothermal assumption for permeate",
        )
        def eq_permeate_isothermal(b, t, x):
            return (
                b.feed_side.properties[t, x].temperature
                == b.permeate_side[t, x].temperature
            )

        # ==========================================================================
        # isothermal conditions at permeate outlet
        @self.Constraint(
            self.flowsheet().config.time, doc="Isothermal assumption for permeate out"
        )
        def eq_permeate_outlet_isothermal(b, t):
            return (
                b.feed_side.properties[t, b.length_domain.last()].temperature
                == b.mixed_permeate[t].temperature
            )

        self.recovery_vol_phase = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            initialize=0.4,
            units=pyunits.dimensionless,
            doc="Volumetric recovery rate",
        )

        # ==========================================================================
        # Volumetric Recovery rate
        @self.Constraint(self.flowsheet().config.time)
        def eq_recovery_vol_phase(b, t):
            return (
                b.recovery_vol_phase[t, "Liq"]
                == b.mixed_permeate[t].flow_vol_phase["Liq"]
                / b.feed_side.properties[t, self.first_element].flow_vol_phase["Liq"]
            )

        solvent_set = self.config.property_package.solvent_set
        solute_set = self.config.property_package.solute_set

        self.recovery_mass_phase_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=lambda b, t, p, j: 0.4037 if j in solvent_set else 0.0033,
            bounds=lambda b, t, p, j: (
                (1e-2, 1 - 1e-6) if j in solvent_set else (1e-5, 1 - 1e-6)
            ),
            units=pyunits.dimensionless,
            doc="Mass-based component recovery",
        )

        self.rejection_phase_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            solute_set,
            initialize=0.9,
            bounds=(1e-2, 1 - 1e-6),
            units=pyunits.dimensionless,
            doc="Observed solute rejection",
        )

        # ==========================================================================
        # Mass-based Component Recovery rate
        @self.Constraint(
            self.flowsheet().config.time, self.config.property_package.component_list
        )
        def eq_recovery_mass_phase_comp(b, t, j):
            return (
                b.recovery_mass_phase_comp[t, "Liq", j]
                * b.feed_side.properties[t, b.first_element].flow_mass_phase_comp[
                    "Liq", j
                ]
                == b.mixed_permeate[t].flow_mass_phase_comp["Liq", j]
            )

        # rejection
        @self.Constraint(self.flowsheet().config.time, solute_set)
        def eq_rejection_phase_comp(b, t, j):
            return b.rejection_phase_comp[t, "Liq", j] == 1 - (
                b.mixed_permeate[t].conc_mass_phase_comp["Liq", j]
                / b.feed_side.properties[t, self.first_element].conc_mass_phase_comp[
                    "Liq", j
                ]
            )

        self._add_flux_balance()
        if self.config.has_pressure_change:
            self._add_deltaP()
        self._add_mass_transfer()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

    def _add_deltaP(self):
        raise NotImplementedError()

    def _add_mass_transfer(self):
        raise NotImplementedError()

    def _add_length_and_width(self):
        units_meta = self.config.property_package.get_metadata().get_derived_units
        self.length = Var(
            initialize=10,
            bounds=(0.1, 5e2),
            domain=NonNegativeReals,
            units=units_meta("length"),
            doc="Effective membrane length",
        )
        self.width = Var(
            initialize=1,
            bounds=(1e-1, 1e3),
            domain=NonNegativeReals,
            units=units_meta("length"),
            doc="Membrane width",
        )

    def _add_area(self, include_constraint=True):
        units_meta = self.config.property_package.get_metadata().get_derived_units

        self.area = Var(
            initialize=10,
            bounds=(1e-1, 1e5),
            domain=NonNegativeReals,
            units=units_meta("length") ** 2,
            doc="Total Membrane area",
        )

        if include_constraint:
            if self.config.module_type == ModuleType.flat_sheet:
                # Membrane area equation for flat plate membranes
                @self.Constraint(doc="Total Membrane area")
                def eq_area(b):
                    return b.area == b.length * b.width

            elif self.config.module_type == ModuleType.spiral_wound:
                # Membrane area equation
                @self.Constraint(doc="Total Membrane area")
                def eq_area(b):
                    return b.area == b.length * 2 * b.width

            else:
                raise ConfigurationError(
                    "Unsupported membrane module type: {}".format(
                        self.config.module_type
                    )
                )

    def _add_flux_balance(self):

        solvent_set = self.config.property_package.solvent_set
        solute_set = self.config.property_package.solute_set

        units_meta = self.config.property_package.get_metadata().get_derived_units

        self.A_comp = Var(
            self.flowsheet().config.time,
            solvent_set,
            initialize=1e-12,
            bounds=(1e-18, 1e-6),
            domain=NonNegativeReals,
            units=units_meta("length")
            * units_meta("pressure") ** -1
            * units_meta("time") ** -1,
            doc="Solvent permeability coeff.",
        )

        self.B_comp = Var(
            self.flowsheet().config.time,
            solute_set,
            initialize=1e-8,
            bounds=(1e-11, 1e-5),
            domain=NonNegativeReals,
            units=units_meta("length") * units_meta("time") ** -1,
            doc="Solute permeability coeff.",
        )

        # TODO: add water density to NaCl prop model and remove here (or use IDAES version)
        self.dens_solvent = Param(
            initialize=1000,
            units=units_meta("mass") * units_meta("length") ** -3,
            doc="Pure water density",
        )

        if self.config.transport_model == TransportModel.SKK:
            self.reflect_coeff = Var(
                initialize=0.9,
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Reflection coefficient of the membrane",
            )

            self.alpha = Var(
                initialize=1e8,
                domain=NonNegativeReals,
                units=units_meta("time") * units_meta("length") ** -1,
                doc="Alpha coefficient of the membrane",
            )

        self.flux_mass_phase_comp = Var(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=lambda b, t, x, p, j: 5e-4 if j in solvent_set else 1e-6,
            bounds=lambda b, t, x, p, j: (
                (1e-4, 3e-2) if j in solvent_set else (1e-8, 1e-3)
            ),
            units=units_meta("mass")
            * units_meta("length") ** -2
            * units_meta("time") ** -1,
            doc="Mass flux across membrane at inlet and outlet",
        )

        if self.config.transport_model == TransportModel.SD:

            @self.Constraint(
                self.flowsheet().config.time,
                self.difference_elements,
                self.config.property_package.phase_list,
                self.config.property_package.component_list,
                doc="Solvent and solute mass flux using SD model",
            )
            def eq_flux_mass(b, t, x, p, j):
                prop_feed = b.feed_side.properties[t, x]
                prop_perm = b.permeate_side[t, x]
                interface = b.feed_side.properties_interface[t, x]
                comp = self.config.property_package.get_component(j)
                if comp.is_solvent():
                    return b.flux_mass_phase_comp[t, x, p, j] == b.A_comp[
                        t, j
                    ] * b.dens_solvent * (
                        (prop_feed.pressure - prop_perm.pressure)
                        - (
                            interface.pressure_osm_phase[p]
                            - prop_perm.pressure_osm_phase[p]
                        )
                    )
                elif comp.is_solute():
                    return b.flux_mass_phase_comp[t, x, p, j] == b.B_comp[t, j] * (
                        interface.conc_mass_phase_comp[p, j]
                        - prop_perm.conc_mass_phase_comp[p, j]
                    )

        elif self.config.transport_model == TransportModel.SKK:

            @self.Constraint(
                self.flowsheet().config.time, solute_set, doc="SKK alpha coeff."
            )
            def eq_alpha(b, t, j):
                return b.alpha == (1 - b.reflect_coeff) / b.B_comp[t, j]

            @self.Constraint(
                self.flowsheet().config.time,
                self.difference_elements,
                self.config.property_package.phase_list,
                self.config.property_package.component_list,
                doc="Solvent and solute mass flux using SKK model",
            )
            def eq_flux_mass(b, t, x, p, j):
                prop_feed = b.feed_side.properties[t, x]
                prop_perm = b.permeate_side[t, x]
                interface = b.feed_side.properties_interface[t, x]
                comp = self.config.property_package.get_component(j)
                if comp.is_solvent():
                    return b.flux_mass_phase_comp[t, x, p, j] == b.A_comp[
                        t, j
                    ] * b.dens_solvent * (
                        (prop_feed.pressure - prop_perm.pressure)
                        - b.reflect_coeff
                        * (
                            interface.pressure_osm_phase[p]
                            - prop_perm.pressure_osm_phase[p]
                        )
                    )
                elif comp.is_solute():
                    return b.flux_mass_phase_comp[t, x, p, j] == b.B_comp[t, j] * (
                        interface.conc_mass_phase_comp[p, j]
                        - prop_perm.conc_mass_phase_comp[p, j]
                    ) + (1 - b.reflect_coeff) * (
                        (
                            (b.flux_mass_phase_comp[t, x, p, "H2O"] / b.dens_solvent)
                            * interface.conc_mass_phase_comp[p, j]
                        )
                    )

        else:
            raise ConfigurationError(
                "Unsupported transport model: {}".format(self.config.transport_model)
            )

        @self.Expression(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Average flux expression",
        )
        def flux_mass_phase_comp_avg(b, t, p, j):
            return (
                sum(
                    b.flux_mass_phase_comp[t, x, p, j] for x in self.difference_elements
                )
                / self.nfe
            )

        if (
            self.config.concentration_polarization_type
            == ConcentrationPolarizationType.calculated
        ):

            @self.Constraint(
                self.flowsheet().config.time,
                self.difference_elements,
                solute_set,
                doc="Concentration polarization",
            )
            def eq_concentration_polarization(b, t, x, j):
                jw = b.flux_mass_phase_comp[t, x, "Liq", "H2O"] / self.dens_solvent
                js = b.flux_mass_phase_comp[t, x, "Liq", j]
                return b.feed_side.properties_interface[t, x].conc_mass_phase_comp[
                    "Liq", j
                ] == b.feed_side.properties[t, x].conc_mass_phase_comp["Liq", j] * exp(
                    jw / self.feed_side.K[t, x, j]
                ) - js / jw * (
                    exp(jw / self.feed_side.K[t, x, j]) - 1
                )

        return self.eq_flux_mass

    def _get_state_args_permeate(self, initialize_guess, state_args):
        """
        Arguments:
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

        source = self.feed_side.properties[
            self.flowsheet().config.time.first(), self.first_element
        ]
        mixed_permeate_properties = self.mixed_permeate[
            self.flowsheet().config.time.first()
        ]

        # assumptions
        if initialize_guess is None:
            initialize_guess = {}
        if "solvent_recovery" not in initialize_guess:
            initialize_guess["solvent_recovery"] = 0.5
        if "solute_recovery" not in initialize_guess:
            initialize_guess["solute_recovery"] = 0.01

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
        state_args_permeate = deepcopy(state_args)

        state_args_permeate["pressure"] = mixed_permeate_properties.pressure.value
        for j in self.config.property_package.solvent_set:
            state_args_permeate["flow_mass_phase_comp"][("Liq", j)] *= initialize_guess[
                "solvent_recovery"
            ]
        for j in self.config.property_package.solute_set:
            state_args_permeate["flow_mass_phase_comp"][("Liq", j)] *= initialize_guess[
                "solute_recovery"
            ]

        return state_args_permeate

    def initialize_build(
        self,
        initialize_guess=None,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        General wrapper for RO initialization routines

        Keyword Arguments:

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
                         property package) (default = None).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None)
            solver : solver object or string indicating which solver to use during
                     initialization, if None provided the default solver will be used
                     (default = None)
        Returns:
            None
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        source_flags = self.feed_side.initialize(
            state_args=state_args,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            initialize_guess=initialize_guess,
        )

        init_log.info_high("Initialization Step 1a (feed side) Complete")

        state_args_permeate = self._get_state_args_permeate(
            initialize_guess, state_args
        )

        self.permeate_side.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_permeate,
        )

        self.mixed_permeate.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_permeate,
        )

        init_log.info_high("Initialization Step 1b (permeate side) Complete")

        if degrees_of_freedom(self) != 0:
            # TODO: should we have a separate error for DoF?
            raise Exception(
                f"{self.name} degrees of freedom were not 0 at the beginning "
                f"of initialization. DoF = {degrees_of_freedom(self)}"
            )

        # pre-solve using interval arithmetic
        interval_initializer(self)

        # Create solver
        opt = get_solver(solver, optarg)

        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
            # occasionally it might be worth retrying a solve
            if not check_optimal_termination(res):
                init_log.warning(
                    f"Trouble solving unit model {self.name}, trying one more time"
                )
                res = opt.solve(self, tee=slc.tee)
        init_log.info_high(f"Initialization Step 2 {idaeslog.condition(res)}")
        # release inlet state, in case this error is caught

        self.feed_side.release_state(source_flags, outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize")

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Feed Inlet": self.inlet,
                "Feed Outlet": self.retentate,
                "Permeate Outlet": self.permeate,
            },
            time_point=time_point,
        )

    def _get_performance_contents(self, time_point=0):
        x_in = self.first_element
        x_interface_in = self.difference_elements.first()
        x_out = self.length_domain.last()
        feed_inlet = self.feed_side.properties[time_point, x_in]
        feed_outlet = self.feed_side.properties[time_point, x_out]
        interface_inlet = self.feed_side.properties_interface[
            time_point, x_interface_in
        ]
        interface_outlet = self.feed_side.properties_interface[time_point, x_out]
        permeate = self.mixed_permeate[time_point]
        var_dict = {}
        expr_dict = {}
        var_dict["Volumetric Recovery Rate"] = self.recovery_vol_phase[
            time_point, "Liq"
        ]
        var_dict["Solvent Mass Recovery Rate"] = self.recovery_mass_phase_comp[
            time_point, "Liq", "H2O"
        ]
        var_dict["Membrane Area"] = self.area
        if hasattr(self, "length") and self.config.has_full_reporting:
            var_dict["Membrane Length"] = self.length
        if hasattr(self, "width") and self.config.has_full_reporting:
            var_dict["Membrane Width"] = self.width
        if hasattr(self, "deltaP") and self.config.has_full_reporting:
            var_dict["Pressure Change"] = self.deltaP[time_point]
        if hasattr(self.feed_side, "N_Re") and self.config.has_full_reporting:
            var_dict["Reynolds Number @Inlet"] = self.feed_side.N_Re[time_point, x_in]
            var_dict["Reynolds Number @Outlet"] = self.feed_side.N_Re[time_point, x_out]
        if hasattr(self.feed_side, "velocity") and self.config.has_full_reporting:
            var_dict["Velocity @Inlet"] = self.feed_side.velocity[time_point, x_in]
            var_dict["Velocity @Outlet"] = self.feed_side.velocity[time_point, x_out]
        for j in self.config.property_package.solute_set:
            if (
                interface_inlet.is_property_constructed("conc_mass_phase_comp")
                and self.config.has_full_reporting
            ):
                var_dict[f"{j} Concentration @Inlet,Membrane-Interface "] = (
                    interface_inlet.conc_mass_phase_comp["Liq", j]
                )
            if (
                interface_outlet.is_property_constructed("conc_mass_phase_comp")
                and self.config.has_full_reporting
            ):
                var_dict[f"{j} Concentration @Outlet,Membrane-Interface "] = (
                    interface_outlet.conc_mass_phase_comp["Liq", j]
                )
            if (
                feed_inlet.is_property_constructed("conc_mass_phase_comp")
                and self.config.has_full_reporting
            ):
                var_dict[f"{j} Concentration @Inlet,Bulk"] = (
                    feed_inlet.conc_mass_phase_comp["Liq", j]
                )
            if (
                feed_outlet.is_property_constructed("conc_mass_phase_comp")
                and self.config.has_full_reporting
            ):
                var_dict[f"{j} Concentration @Outlet,Bulk"] = (
                    feed_outlet.conc_mass_phase_comp["Liq", j]
                )
            if (
                permeate.is_property_constructed("conc_mass_phase_comp")
                and self.config.has_full_reporting
            ):
                var_dict[f"{j} Permeate Concentration"] = permeate.conc_mass_phase_comp[
                    "Liq", j
                ]
        if (
            interface_outlet.is_property_constructed("pressure_osm_phase")
            and self.config.has_full_reporting
        ):
            var_dict["Osmotic Pressure @Outlet,Membrane-Interface "] = (
                interface_outlet.pressure_osm_phase["Liq"]
            )
        if (
            feed_outlet.is_property_constructed("pressure_osm_phase")
            and self.config.has_full_reporting
        ):
            var_dict["Osmotic Pressure @Outlet,Bulk"] = feed_outlet.pressure_osm_phase[
                "Liq"
            ]
        if (
            interface_inlet.is_property_constructed("pressure_osm_phase")
            and self.config.has_full_reporting
        ):
            var_dict["Osmotic Pressure @Inlet,Membrane-Interface"] = (
                interface_inlet.pressure_osm_phase["Liq"]
            )
        if (
            feed_inlet.is_property_constructed("pressure_osm_phase")
            and self.config.has_full_reporting
        ):
            var_dict["Osmotic Pressure @Inlet,Bulk"] = feed_inlet.pressure_osm_phase[
                "Liq"
            ]
        if (
            feed_inlet.is_property_constructed("flow_vol_phase")
            and self.config.has_full_reporting
        ):
            var_dict["Volumetric Flowrate @Inlet"] = feed_inlet.flow_vol_phase["Liq"]
        if (
            feed_outlet.is_property_constructed("flow_vol_phase")
            and self.config.has_full_reporting
        ):
            var_dict["Volumetric Flowrate @Outlet"] = feed_outlet.flow_vol_phase["Liq"]
        if hasattr(self.feed_side, "dh") and self.config.has_full_reporting:
            var_dict["Hydraulic Diameter"] = self.feed_side.dh

        if self.config.has_full_reporting:
            expr_dict["Average Solvent Mass Flux"] = self.flux_mass_phase_comp_avg[
                time_point, "Liq", "H2O"
            ]
            if hasattr(self.feed_side, "N_Re_avg"):
                expr_dict["Average Reynolds Number"] = self.feed_side.N_Re_avg[
                    time_point
                ]
            for j in self.config.property_package.solute_set:
                expr_dict[f"{j} Average Solute Mass Flux"] = (
                    self.flux_mass_phase_comp_avg[time_point, "Liq", j]
                )
                if hasattr(self.feed_side, "K_avg"):
                    expr_dict[f"{j} Average Mass Transfer Coefficient"] = (
                        self.feed_side.K_avg[time_point, j]
                    )

        # TODO: add more vars
        return {"vars": var_dict, "exprs": expr_dict}

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

        if self.config.transport_model == TransportModel.SKK:
            if iscale.get_scaling_factor(self.alpha) is None:
                iscale.set_scaling_factor(self.alpha, 1e-8)

            if iscale.get_scaling_factor(self.reflect_coeff) is None:
                iscale.set_scaling_factor(self.reflect_coeff, 1)

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

        if not hasattr(self, "_permeate_scaled_properties"):
            self._permeate_scaled_properties = ComponentSet()

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

    @property
    def default_costing_method(self):
        return cost_reverse_osmosis
