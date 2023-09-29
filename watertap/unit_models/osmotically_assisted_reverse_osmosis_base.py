#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import math
from pyomo.common.config import Bool, ConfigValue
from pyomo.environ import (
    NonNegativeReals,
    Param,
    Suffix,
    Var,
    Constraint,
    check_optimal_termination,
    exp,
    units as pyunits,
    value,
)
from idaes.core import UnitModelBlockData, FlowDirection
from idaes.core.solvers import get_solver
from idaes.core.util import scaling as iscale
from idaes.core.util.exceptions import ConfigurationError, InitializationError
from idaes.core.util.misc import add_object_reference
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.tables import create_stream_table_dataframe
import idaes.logger as idaeslog

from watertap.core.membrane_channel_base import (
    validate_membrane_config_args,
    ConcentrationPolarizationType,
)

from watertap.core import InitializationMixin
from watertap.costing.unit_models.osmotically_assisted_reverse_osmosis import (
    cost_osmotically_assisted_reverse_osmosis,
)


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


class OsmoticallyAssistedReverseOsmosisBaseData(
    InitializationMixin, UnitModelBlockData
):
    """
    Osmotically Assisted Reverse Osmosis base class

    """

    def build(self):
        """
        Common variables and constraints for an OARO unit model

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

        # Raise exception if any of configuration arguments are provided incorrectly
        validate_membrane_config_args(self)

        # --------------------------------------------------------------
        # Add feed side MembraneChannel Control Volume and setup
        self._add_membrane_channel_and_geometry(
            side="feed_side", flow_direction=FlowDirection.forward
        )

        self.feed_side.add_state_blocks(has_phase_equilibrium=False)

        self.feed_side.add_material_balances(
            balance_type=self.config.material_balance_type, has_mass_transfer=True
        )

        self.feed_side.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            pressure_change_type=self.config.pressure_change_type,
            has_pressure_change=self.config.has_pressure_change,
            friction_factor=self.config.friction_factor,
        )

        # Add constraint for equal temperature between bulk and interface
        self.feed_side.add_interface_isothermal_conditions()

        # Add constraint for equal temperature between inlet and outlet
        self.feed_side.add_control_volume_isothermal_conditions()

        # Add constraint for volumetric flow equality between interface and bulk
        self.feed_side.add_extensive_flow_to_interface()

        # Concentration polarization constraint is not accounted for in the below method; it is
        # written later in the base model (eq_concentration_polarization)

        self.feed_side.add_concentration_polarization(
            concentration_polarization_type=self.config.concentration_polarization_type,
            mass_transfer_coefficient=self.config.mass_transfer_coefficient,
        )

        # Pass in 0D, applied in 1D
        self.feed_side.apply_transformation()

        self.feed_side.add_expressions()

        # --------------------------------------------------------------
        # Add permeate side MembraneChannel Control Volume and setup
        self._add_membrane_channel_and_geometry(
            side="permeate_side", flow_direction=FlowDirection.backward
        )

        self.permeate_side.add_state_blocks(has_phase_equilibrium=False)

        self.permeate_side.add_material_balances(
            balance_type=self.config.material_balance_type, has_mass_transfer=True
        )

        self.permeate_side.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            pressure_change_type=self.config.pressure_change_type,
            has_pressure_change=self.config.has_pressure_change,
            friction_factor=self.config.friction_factor,
        )

        # Add constraint for equal temperature between bulk and interface
        self.permeate_side.add_interface_isothermal_conditions()

        # NOTE: We do *not* add a constraint for equal temperature
        #       between inlet and outlet, that is handled by
        #       eq_permeate_isothermal below and checks in initialization

        # Add constraint for volumetric flow equality between interface and bulk
        self.permeate_side.add_extensive_flow_to_interface()

        # Concentration polarization constraint is not accounted for in the below method; it is
        # written later in the base model (eq_concentration_polarization)
        self.permeate_side.add_concentration_polarization(
            concentration_polarization_type=self.config.concentration_polarization_type,
            mass_transfer_coefficient=self.config.mass_transfer_coefficient,
        )

        # Pass in 0D, applied in 1D
        self.permeate_side.apply_transformation()

        self.permeate_side.add_expressions()

        add_object_reference(self, "length_domain", self.feed_side.length_domain)
        add_object_reference(
            self, "difference_elements", self.feed_side.difference_elements
        )
        add_object_reference(self, "first_element", self.feed_side.first_element)
        add_object_reference(self, "nfe", self.feed_side.nfe)

        # Add Ports
        self.add_inlet_port(name="feed_inlet", block=self.feed_side)
        self.add_outlet_port(name="feed_outlet", block=self.feed_side)
        self.add_inlet_port(name="permeate_inlet", block=self.permeate_side)
        self.add_outlet_port(name="permeate_outlet", block=self.permeate_side)

        # # ==========================================================================
        # Feed and permeate-side isothermal conditions
        @self.Constraint(
            self.flowsheet().config.time,
            self.length_domain,
            doc="Isothermal assumption for permeate",
        )
        def eq_permeate_isothermal(b, t, x):
            if x == self.length_domain.last():
                return Constraint.Skip
            return (
                b.permeate_side.properties[t, x].temperature
                == 0.5
                * b.feed_side.properties[
                    t, b.feed_side.length_domain.first()
                ].temperature
                + 0.5
                * b.permeate_side.properties[
                    t, b.permeate_side.length_domain.last()
                ].temperature
            )

        # ==========================================================================
        # Volumetric Recovery rate
        self.recovery_vol_phase = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            initialize=0.4,
            bounds=(1e-2, 1 - 1e-6),
            units=pyunits.dimensionless,
            doc="Volumetric recovery rate",
        )

        # ==========================================================================
        @self.Constraint(self.flowsheet().config.time)
        def eq_recovery_vol_phase(b, t):
            return (
                b.recovery_vol_phase[t, "Liq"]
                == (
                    b.permeate_side.properties[t, b.first_element].flow_vol_phase["Liq"]
                    - b.permeate_side.properties[
                        t, b.permeate_side.length_domain.last()
                    ].flow_vol_phase["Liq"]
                )
                / b.feed_side.properties[t, b.first_element].flow_vol_phase["Liq"]
            )

        solvent_set = self.config.property_package.solvent_set
        solute_set = self.config.property_package.solute_set

        self.recovery_mass_phase_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=lambda b, t, p, j: 0.4037 if j in solvent_set else 0.0033,
            bounds=lambda b, t, p, j: (0, 1 - 1e-6)
            if j in solvent_set
            else (1e-5, 1 - 1e-6),
            units=pyunits.dimensionless,
            doc="Mass-based component recovery",
        )

        self.rejection_phase_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            solute_set,
            initialize=0.9,
            bounds=(0, 1 - 1e-6),
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
                == b.permeate_side.properties[t, b.first_element].flow_mass_phase_comp[
                    "Liq", j
                ]
                - b.permeate_side.properties[
                    t, b.permeate_side.length_domain.last()
                ].flow_mass_phase_comp["Liq", j]
            )

        # rejection
        # TODO: consider importance of rejection in OARO; for now,
        # using outlet permeate concentration in the constraint; could calculate
        # as a function of length but need to decide on whether that's overkill
        @self.Constraint(self.flowsheet().config.time, solute_set)
        def eq_rejection_phase_comp(b, t, j):
            return b.rejection_phase_comp[t, "Liq", j] == 1 - (
                b.permeate_side.properties[t, b.first_element].conc_mass_phase_comp[
                    "Liq", j
                ]
                / b.feed_side.properties[t, b.first_element].conc_mass_phase_comp[
                    "Liq", j
                ]
            )

        self._add_flux_balance()

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
        if not hasattr(self, "area"):
            self.area = Var(
                initialize=10,
                bounds=(1e-1, 1e5),
                domain=NonNegativeReals,
                units=units_meta("length") ** 2,
                doc="Total Membrane area",
            )

        if include_constraint:
            if not hasattr(self, "eq_area"):
                # Membrane area equation
                @self.Constraint(doc="Total Membrane area")
                def eq_area(b):
                    return b.area == b.length * b.width

            else:
                raise ValueError(
                    "include_constraint was set to True inside of _add_area(), but area constraint already exists."
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

        self.flux_mass_phase_comp = Var(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=lambda b, t, x, p, j: 5e-4 if j in solvent_set else 1e-6,
            bounds=lambda b, t, x, p, j: (0, 3e-2) if j in solvent_set else (0, 1e-3),
            units=units_meta("mass")
            * units_meta("length") ** -2
            * units_meta("time") ** -1,
            doc="Mass flux across membrane at inlet and outlet",
        )

        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Solvent and solute mass flux",
        )
        def eq_flux_mass(b, t, x, p, j):
            prop_feed = b.feed_side.properties[t, x]
            prop_perm = b.permeate_side.properties[t, x]
            interface_feed = b.feed_side.properties_interface[t, x]
            interface_perm = b.permeate_side.properties_interface[t, x]
            comp = self.config.property_package.get_component(j)
            if comp.is_solvent():
                return b.flux_mass_phase_comp[t, x, p, j] == b.A_comp[
                    t, j
                ] * b.dens_solvent * (
                    (prop_feed.pressure - prop_perm.pressure)
                    - (
                        interface_feed.pressure_osm_phase[p]
                        - interface_perm.pressure_osm_phase[p]
                    )
                )
            elif comp.is_solute():
                return b.flux_mass_phase_comp[t, x, p, j] == b.B_comp[t, j] * (
                    interface_feed.conc_mass_phase_comp[p, j]
                    - interface_perm.conc_mass_phase_comp[p, j]
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
                doc="External Concentration polarization in feed side",
            )
            def eq_concentration_polarization_feed(b, t, x, j):
                jw = b.flux_mass_phase_comp[t, x, "Liq", "H2O"] / self.dens_solvent
                js = b.flux_mass_phase_comp[t, x, "Liq", j]
                exponent = jw / self.feed_side.K[t, x, j]
                return b.feed_side.properties_interface[t, x].conc_mass_phase_comp[
                    "Liq", j
                ] == b.feed_side.properties[t, x].conc_mass_phase_comp["Liq", j] * exp(
                    exponent
                ) - js / jw * (
                    exp(exponent) - 1
                )

            self.structural_parameter = Var(
                initialize=1200e-6,
                units=pyunits.m,
                domain=NonNegativeReals,
                doc="Membrane structural parameter (i.e., effective thickness)",
            )

            @self.Constraint(
                self.flowsheet().config.time,
                self.difference_elements,
                solute_set,
                doc="Internal Concentration polarization in permeate side",
            )
            def eq_concentration_polarization_permeate(b, t, x, j):
                jw = b.flux_mass_phase_comp[t, x, "Liq", "H2O"] / self.dens_solvent
                js = b.flux_mass_phase_comp[t, x, "Liq", j]
                exponent = -jw * (
                    b.structural_parameter
                    / b.permeate_side.properties[t, x].diffus_phase_comp["Liq", j]
                    + 1 / b.permeate_side.K[t, x, j]
                )
                return b.permeate_side.properties_interface[t, x].conc_mass_phase_comp[
                    "Liq", j
                ] == b.permeate_side.properties[t, x].conc_mass_phase_comp[
                    "Liq", j
                ] * exp(
                    exponent
                ) - js / jw * (
                    exp(exponent) - 1
                )

        return self.eq_flux_mass

    def initialize_build(
        self,
        initialize_guess=None,
        state_args_feed=None,
        state_args_permeate=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
        raise_on_isothermal_violation=True,
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
            state_args_feed : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for the inlet
                         feed side state block (see documentation of the specific
                         property package) (default = None).
            state_args_permeate : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for the inlet
                         permeate side state block (see documentation of the specific
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

        for t in self.flowsheet().config.time:
            if not math.isclose(
                value(self.feed_inlet.temperature[t]),
                value(self.permeate_inlet.temperature[t]),
            ):
                msg = f"Feed temperatures are different at time {t}, but OARO makes isothermal assumption"

                if raise_on_isothermal_violation:
                    raise InitializationError(msg)
                else:
                    init_log.warning(msg)

            # set them equal
            self.permeate_inlet.temperature[t].set_value(
                value(self.feed_inlet.temperature[t])
            )

        feed_flags = self.feed_side.initialize(
            state_args=state_args_feed,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            initialize_guess=initialize_guess,
        )

        init_log.info_high("Initialization Step 1a (feed side) Complete")

        permeate_flags = self.permeate_side.initialize(
            state_args=state_args_permeate,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            initialize_guess=initialize_guess,
        )

        init_log.info_high("Initialization Step 1b (permeate side) Complete")

        if degrees_of_freedom(self) != 0:
            # TODO: should we have a separate error for DoF?
            raise Exception(
                f"{self.name} degrees of freedom were not 0 at the beginning "
                f"of initialization. DoF = {degrees_of_freedom(self)}"
            )

        # Create solver
        opt = get_solver(solver, optarg)

        # Solve unit *without* flux equation
        self.eq_flux_mass.deactivate()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high(f"Initialization Step 2 {idaeslog.condition(res)}")

        # Solve unit *with* flux equation
        self.eq_flux_mass.activate()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high(f"Initialization Step 3 {idaeslog.condition(res)}")

        # release inlet state, in case this error is caught
        self.permeate_side.release_state(permeate_flags, outlvl)
        self.feed_side.release_state(feed_flags, outlvl)

        init_log.info(f"Initialization Complete: {idaeslog.condition(res)}")

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize")

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Feed Inlet": self.feed_inlet,
                "Feed Outlet": self.feed_outlet,
                "Permeate Inlet": self.permeate_inlet,
                "Permeate Outlet": self.permeate_outlet,
            },
            time_point=time_point,
        )

    # TODO: sort out first/last element for feed/permeate
    def _get_performance_contents(self, time_point=0):
        x_in = self.first_element
        x_interface_in = self.difference_elements.first()
        x_out = self.length_domain.last()
        feed_inlet = self.feed_side.properties[time_point, x_in]
        feed_outlet = self.feed_side.properties[time_point, x_out]
        feed_interface_inlet = self.feed_side.properties_interface[
            time_point, x_interface_in
        ]
        feed_interface_outlet = self.feed_side.properties_interface[time_point, x_out]
        permeate_inlet = self.permeate_side.properties[time_point, x_in]
        permeate_outlet = self.permeate_side.properties[time_point, x_out]
        permeate_interface_inlet = self.permeate_side.properties_interface[
            time_point, x_interface_in
        ]
        permeate_interface_outlet = self.permeate_side.properties_interface[
            time_point, x_out
        ]
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
                feed_interface_inlet.is_property_constructed("conc_mass_phase_comp")
                and self.config.has_full_reporting
            ):
                var_dict[
                    f"{j} Feed Concentration @Inlet,Membrane-Interface "
                ] = feed_interface_inlet.conc_mass_phase_comp["Liq", j]
            if (
                feed_interface_outlet.is_property_constructed("conc_mass_phase_comp")
                and self.config.has_full_reporting
            ):
                var_dict[
                    f"{j} Feed Concentration @Outlet,Membrane-Interface "
                ] = feed_interface_outlet.conc_mass_phase_comp["Liq", j]
            if (
                feed_inlet.is_property_constructed("conc_mass_phase_comp")
                and self.config.has_full_reporting
            ):
                var_dict[
                    f"{j} Feed Concentration @Inlet,Bulk"
                ] = feed_inlet.conc_mass_phase_comp["Liq", j]
            if (
                feed_outlet.is_property_constructed("conc_mass_phase_comp")
                and self.config.has_full_reporting
            ):
                var_dict[
                    f"{j} Feed Concentration @Outlet,Bulk"
                ] = feed_outlet.conc_mass_phase_comp["Liq", j]
            if (
                permeate_interface_inlet.is_property_constructed("conc_mass_phase_comp")
                and self.config.has_full_reporting
            ):
                var_dict[
                    f"{j} Permeate Concentration @Inlet,Membrane-Interface "
                ] = permeate_interface_inlet.conc_mass_phase_comp["Liq", j]
            if (
                permeate_interface_outlet.is_property_constructed(
                    "conc_mass_phase_comp"
                )
                and self.config.has_full_reporting
            ):
                var_dict[
                    f"{j} Permeate Concentration @Outlet,Membrane-Interface "
                ] = permeate_interface_outlet.conc_mass_phase_comp["Liq", j]
            if (
                permeate_inlet.is_property_constructed("conc_mass_phase_comp")
                and self.config.has_full_reporting
            ):
                var_dict[
                    f"{j} Permeate Concentration @Inlet,Bulk"
                ] = permeate_inlet.conc_mass_phase_comp["Liq", j]
            if (
                permeate_outlet.is_property_constructed("conc_mass_phase_comp")
                and self.config.has_full_reporting
            ):
                var_dict[
                    f"{j} Permeate Concentration @Outlet,Bulk"
                ] = permeate_outlet.conc_mass_phase_comp["Liq", j]
        if (
            feed_interface_outlet.is_property_constructed("pressure_osm_phase")
            and self.config.has_full_reporting
        ):
            var_dict[
                "Feed Osmotic Pressure @Outlet,Membrane-Interface "
            ] = feed_interface_outlet.pressure_osm_phase["Liq"]
        if (
            permeate_outlet.is_property_constructed("pressure_osm_phase")
            and self.config.has_full_reporting
        ):
            var_dict[
                "Feed Osmotic Pressure @Outlet,Bulk"
            ] = feed_outlet.pressure_osm_phase["Liq"]
        if (
            feed_interface_inlet.is_property_constructed("pressure_osm_phase")
            and self.config.has_full_reporting
        ):
            var_dict[
                "Feed Osmotic Pressure @Inlet,Membrane-Interface"
            ] = feed_interface_inlet.pressure_osm_phase["Liq"]
        if (
            feed_inlet.is_property_constructed("pressure_osm_phase")
            and self.config.has_full_reporting
        ):
            var_dict[
                "Feed Osmotic Pressure @Inlet,Bulk"
            ] = feed_inlet.pressure_osm_phase["Liq"]
        # TODO: add all corresponding values for permeate side for relevant
        #  vars/expressions from osmotic pressure and whatever is below
        if (
            feed_inlet.is_property_constructed("flow_vol_phase")
            and self.config.has_full_reporting
        ):
            var_dict["Feed Volumetric Flowrate @Inlet"] = feed_inlet.flow_vol_phase[
                "Liq"
            ]
        if (
            feed_outlet.is_property_constructed("flow_vol_phase")
            and self.config.has_full_reporting
        ):
            var_dict["Feed Volumetric Flowrate @Outlet"] = feed_outlet.flow_vol_phase[
                "Liq"
            ]
        if hasattr(self.feed_side, "dh") and self.config.has_full_reporting:
            var_dict["Feed-side Hydraulic Diameter"] = self.feed_side.dh

        if self.config.has_full_reporting:
            expr_dict["Average Solvent Flux (LMH)"] = (
                self.flux_mass_phase_comp_avg[time_point, "Liq", "H2O"] * 3.6e3
            )
            if hasattr(self.feed_side, "N_Re_avg"):
                expr_dict[
                    "Average Feed-side Reynolds Number"
                ] = self.feed_side.N_Re_avg[time_point]
            for j in self.config.property_package.solute_set:
                expr_dict[f"{j} Average Solute Flux (GMH)"] = (
                    self.flux_mass_phase_comp_avg[time_point, "Liq", j] * 3.6e6
                )
                if hasattr(self.feed_side, "K_avg"):
                    expr_dict[
                        f"{j} Average Feed-side Mass Transfer Coefficient (mm/h)"
                    ] = (self.feed_side.K_avg[time_point, j] * 3.6e6)

        # TODO: add more vars
        return {"vars": var_dict, "exprs": expr_dict}

    def calculate_scaling_factors(self):

        if iscale.get_scaling_factor(self.dens_solvent) is None:
            sf = iscale.get_scaling_factor(
                self.feed_side.properties[0, self.first_element].dens_mass_phase["Liq"]
            )
            iscale.set_scaling_factor(self.dens_solvent, sf)

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

        if hasattr(self, "structural_parameter"):
            if iscale.get_scaling_factor(self.structural_parameter) is None:
                # Structural parameter expected to be ~ 1200 microns
                iscale.set_scaling_factor(self.structural_parameter, 1e3)

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
                sf = iscale.get_scaling_factor(v)
                iscale.constraint_scaling_transform(self.eq_flux_mass[t, x, p, j], sf)

    @property
    def default_costing_method(self):
        return cost_osmotically_assisted_reverse_osmosis
