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


from pyomo.environ import (
    NonNegativeReals,
    Param,
    Var,
    check_optimal_termination,
    exp,
    units as pyunits,
)

from idaes.core import UnitModelBlockData
from idaes.core.solvers import get_solver
from idaes.core.util import scaling as iscale
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.misc import add_object_reference
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.tables import create_stream_table_dataframe
import idaes.logger as idaeslog
from watertap.core import InitializationMixin

from .MD_channel_base import (
    ConcentrationPolarizationType,
    TemperaturePolarizationType,
    MassTransferCoefficient,
)

__author__ = "Elmira Shamlou"


class MembraneDistillationBaseData(InitializationMixin, UnitModelBlockData):
    def build(self):
        """
        Common variables and constraints for a DCMD unit model

        """

        super().build()

        self._make_MD_channel_control_volume("hot_ch", self.config.hot_ch)

        self.hot_ch.add_state_blocks(
            has_phase_equilibrium=False,
            flow_direction=self.config.hot_ch.flow_direction,
            property_package_vapor=self.config.hot_ch.property_package_vapor,
            property_package_args_vapor=self.config.hot_ch.property_package_args_vapor,
        )

        # self.hot_ch.set_config(self.config.hot_ch)

        self.hot_ch.add_material_balances(
            balance_type=self.config.hot_ch.material_balance_type,
            has_mass_transfer=True,
        )

        self.hot_ch.add_momentum_balances(
            balance_type=self.config.hot_ch.momentum_balance_type,
            pressure_change_type=self.config.hot_ch.pressure_change_type,
            has_pressure_change=self.config.hot_ch.has_pressure_change,
            friction_factor=self.config.hot_ch.friction_factor,
        )

        self.hot_ch.add_energy_balances(
            balance_type=self.config.hot_ch.energy_balance_type,
            has_heat_transfer=True,
            has_enthalpy_transfer=True,
        )

        # Add constraint for volumetric flow equality between interface and bulk
        self.hot_ch.add_extensive_flow_to_interface()

        self.hot_ch.add_concentration_polarization(
            concentration_polarization_type=self.config.hot_ch.concentration_polarization_type,
            mass_transfer_coefficient=self.config.hot_ch.mass_transfer_coefficient,
        )

        self.hot_ch.add_temperature_polarization(
            temperature_polarization_type=self.config.hot_ch.temperature_polarization_type,
        )

        self._make_MD_channel_control_volume("cold_ch", self.config.cold_ch)

        self.cold_ch.add_state_blocks(
            has_phase_equilibrium=False,
            flow_direction=self.config.cold_ch.flow_direction,
            property_package_vapor=self.config.cold_ch.property_package_vapor,
            property_package_args_vapor=self.config.cold_ch.property_package_args_vapor,
        )

        # self.cold_ch.set_config(self.config.cold_ch)

        self.cold_ch.add_material_balances(
            balance_type=self.config.cold_ch.material_balance_type,
            has_mass_transfer=True,
        )

        self.cold_ch.add_momentum_balances(
            balance_type=self.config.cold_ch.momentum_balance_type,
            pressure_change_type=self.config.cold_ch.pressure_change_type,
            has_pressure_change=self.config.cold_ch.has_pressure_change,
            friction_factor=self.config.cold_ch.friction_factor,
        )

        self.cold_ch.add_energy_balances(
            balance_type=self.config.cold_ch.energy_balance_type,
            has_heat_transfer=True,
            has_enthalpy_transfer=True,
        )

        # Add constraint for volumetric flow equality between interface and bulk
        self.cold_ch.add_extensive_flow_to_interface()

        # Concentration polarization constraint is not accounted for in the below method; it is
        # written later in the base model (eq_concentration_polarization)
        self.cold_ch.add_concentration_polarization(
            concentration_polarization_type=ConcentrationPolarizationType.none,
            mass_transfer_coefficient=MassTransferCoefficient.none,
        )

        self.cold_ch.add_temperature_polarization(
            temperature_polarization_type=self.config.cold_ch.temperature_polarization_type,
        )
        add_object_reference(self, "length_domain", self.hot_ch.length_domain)
        add_object_reference(
            self, "difference_elements", self.hot_ch.difference_elements
        )
        add_object_reference(self, "first_element", self.hot_ch.first_element)
        add_object_reference(self, "nfe", self.hot_ch.nfe)

        # Add Ports
        self.add_inlet_port(name="hot_ch_inlet", block=self.hot_ch)
        self.add_outlet_port(name="hot_ch_outlet", block=self.hot_ch)
        self.add_inlet_port(name="cold_ch_inlet", block=self.cold_ch)
        self.add_outlet_port(name="cold_ch_outlet", block=self.cold_ch)

        self._add_heat_flux()
        self._add_mass_flux()

        self.recovery_mass = Var(
            self.flowsheet().config.time,
            # self.config.hot_ch.property_package.phase_list,
            # self.config.hot_ch.property_package.component_list,
            initialize=0.1,
            bounds=(1e-11, 0.99),
            units=pyunits.dimensionless,
            doc="Mass-based recovery",
        )

        @self.Constraint(
            self.flowsheet().config.time,
            # self.config.hot_ch.property_package.component_list,
        )
        def eq_recovery_mass(b, t):
            return (
                b.recovery_mass[t]
                * b.hot_ch.properties[t, b.first_element].flow_mass_phase_comp[
                    "Liq", "H2O"
                ]
                == b.cold_ch.properties[t, b.first_element].flow_mass_phase_comp[
                    "Liq", "H2O"
                ]
                - b.cold_ch.properties[
                    t, b.cold_ch.length_domain.last()
                ].flow_mass_phase_comp["Liq", "H2O"]
            )

        self._add_mass_transfer()
        self._add_heat_transfer()

        @self.Expression(
            self.flowsheet().config.time,
            doc="Average thermal efficiency",
        )
        def thermal_efficiency(b, t):
            total_enth_flux = sum(
                b.flux_enth_hot[t, x] for x in self.difference_elements
            )
            total_cond_heat_flux = sum(
                b.flux_conduction_heat[t, x] for x in self.difference_elements
            )

            return total_enth_flux / (total_enth_flux + total_cond_heat_flux)

        @self.Expression(
            self.flowsheet().config.time,
            doc="module heat recovery: ratio of the actual heat recovered by cold side to the maximum ideal heat recovery",
        )
        def effectiveness(b, t):

            return (
                b.cold_ch.properties[t, b.first_element].temperature
                - b.cold_ch.properties[t, b.cold_ch.length_domain.last()].temperature
            ) / (
                b.hot_ch.properties[t, b.first_element].temperature
                - b.cold_ch.properties[t, b.cold_ch.length_domain.last()].temperature
            )

    def _add_mass_flux(self):

        solute_set = self.config.hot_ch.property_package.solute_set

        units_meta = (
            self.config.hot_ch.property_package.get_metadata().get_derived_units
        )

        # todo: add calculation method for permeability coefficient
        self.permeability_coef = Var(
            self.flowsheet().config.time,
            initialize=1e-10,
            bounds=(1e-11, 1e-9),
            units=units_meta("mass")
            * units_meta("time") ** -1
            * units_meta("length") ** -1
            * units_meta("pressure") ** -1,
            doc="membrane permeabiity coeffcient",
        )

        self.flux_mass = Var(
            self.flowsheet().config.time,
            self.difference_elements,
            initialize=1e-10,
            bounds=(1e-12, 1e-1),
            units=units_meta("mass")
            * units_meta("time") ** -1
            * units_meta("length") ** -2,
            doc="solvent vapor mass flux across the membrane",
        )

        self.flux_enth_hot = Var(
            self.flowsheet().config.time,
            self.difference_elements,
            initialize=1e3,
            bounds=(1e-10, 1e5),
            units=pyunits.J * pyunits.s**-1 * pyunits.m**-2,
            doc="hot side evaporation enthalpy flux",
        )

        self.flux_enth_cold = Var(
            self.flowsheet().config.time,
            self.difference_elements,
            initialize=1e3,
            bounds=(1e-10, 1e5),
            units=pyunits.J * pyunits.s**-1 * pyunits.m**-2,
            doc="cold side condensation enthalpy flux",
        )

        self.dens_solvent = Param(
            initialize=1000,
            units=units_meta("mass") * units_meta("length") ** -3,
            doc="Pure water density",
        )

        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            doc="Solvent mass flux",
        )
        def eq_flux_mass(b, t, x):
            return b.flux_mass[t, x] == b.permeability_coef[
                t
            ] / b.membrane_thickness * (
                b.hot_ch.properties_interface[t, x].pressure_sat
                - b.cold_ch.properties_interface[t, x].pressure_sat
            )

        @self.Expression(
            self.flowsheet().config.time,
            doc="Average mass flux expression",
        )
        def flux_mass_avg(b, t):
            return sum(b.flux_mass[t, x] for x in self.difference_elements) / self.nfe

        if (
            self.config.hot_ch.concentration_polarization_type
            == ConcentrationPolarizationType.calculated
        ):

            @self.Constraint(
                self.flowsheet().config.time,
                self.difference_elements,
                solute_set,
                doc="External Concentration polarization in hot side",
            )
            def eq_concentration_polarization_hot_ch(b, t, x, j):
                jw = b.flux_mass[t, x] / self.dens_solvent
                exponent = jw / self.hot_ch.K[t, x, j]
                return b.hot_ch.properties_interface[t, x].conc_mass_phase_comp[
                    "Liq", j
                ] == b.hot_ch.properties[t, x].conc_mass_phase_comp["Liq", j] * exp(
                    exponent
                )

        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            doc="hot side evaporation enthalpy flux",
        )
        def eq_flux_hot_enth(b, t, x):
            return (
                b.flux_enth_hot[t, x]
                == b.flux_mass[t, x]
                * b.hot_ch.properties_vapor[t, x].enth_mass_phase["Vap"]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            doc="cold side evaporation enthalpy flux",
        )
        def eq_flux_cold_enth(b, t, x):
            return (
                b.flux_enth_cold[t, x]
                == b.flux_mass[t, x]
                * b.cold_ch.properties_vapor[t, x].enth_mass_phase["Vap"]
            )

        @self.Expression(
            self.flowsheet().config.time,
            doc="Average hot side enthalpy flux expression",
        )
        def flux_enth_hot_avg(b, t):
            return (
                sum(b.flux_enth_hot[t, x] for x in self.difference_elements) / self.nfe
            )

        @self.Expression(
            self.flowsheet().config.time,
            doc="Average cold side enthalpy flux expression",
        )
        def flux_enth_cold_avg(b, t):
            return (
                sum(b.flux_enth_cold[t, x] for x in self.difference_elements) / self.nfe
            )

        @self.Constraint(
            self.flowsheet().time,
            self.difference_elements,
            doc="hot side Vapor temperature",
        )
        def eq_vapor_temperature_hot(b, t, x):
            return (
                b.hot_ch.properties_vapor[t, x].temperature
                == b.hot_ch.properties_interface[t, x].temperature
            )

        @self.Constraint(
            self.flowsheet().time,
            self.difference_elements,
            doc="hot  side Vapor pressure",
        )
        def eq_vapor_pressure_hot(b, t, x):
            return (
                b.hot_ch.properties_vapor[t, x].pressure
                == b.hot_ch.properties_interface[t, x].pressure_sat
            )

        @self.Constraint(
            self.flowsheet().time,
            self.difference_elements,
            doc="cold side Vapor temperature",
        )
        def eq_vapor_temperature_cold(b, t, x):
            return (
                b.cold_ch.properties_vapor[t, x].temperature
                == b.cold_ch.properties_interface[t, x].temperature
            )

        @self.Constraint(
            self.flowsheet().time,
            self.difference_elements,
            doc="cold  side Vapor pressure",
        )
        def eq_vapor_pressure_cold(b, t, x):
            return (
                b.cold_ch.properties_vapor[t, x].pressure
                == b.cold_ch.properties_interface[t, x].pressure_sat
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            doc="Vapor flow rate and no liquid in vapor state block",
        )
        def eq_vapor_flow(b, t, x):
            lb = b.hot_ch.properties_vapor[t, x].flow_mass_phase_comp["Liq", "H2O"].lb
            b.hot_ch.properties_vapor[t, x].flow_mass_phase_comp["Liq", "H2O"].fix(lb)

            return (
                b.hot_ch.properties_vapor[t, x].flow_mass_phase_comp["Vap", "H2O"]
                == b.flux_mass[t, x] * b.area / self.nfe
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            doc="Equal vapor flow rate between hot_ch and cold_ch",
        )
        def eq_vapor_flow_equal(b, t, x):
            lb = b.hot_ch.properties_vapor[t, x].flow_mass_phase_comp["Liq", "H2O"].lb
            b.hot_ch.properties_vapor[t, x].flow_mass_phase_comp["Liq", "H2O"].fix(lb)
            return (
                b.hot_ch.properties_vapor[t, x].flow_mass_phase_comp["Vap", "H2O"]
                == b.cold_ch.properties_vapor[t, x].flow_mass_phase_comp["Vap", "H2O"]
            )

        # Check for hot channel temperature polarization type
        if (
            self.config.hot_ch.temperature_polarization_type
            != TemperaturePolarizationType.none
        ):

            @self.Constraint(
                self.flowsheet().config.time,
                self.difference_elements,
                doc="Temperature polarization in hot channel",
            )
            def eq_temperature_polarization_hot(b, t, x):
                return (
                    b.hot_ch.h_conv[t, x]
                    * (
                        b.hot_ch.properties[t, x].temperature
                        - b.hot_ch.properties_interface[t, x].temperature
                    )
                    == b.flux_conduction_heat[t, x]
                    + b.flux_enth_hot[t, x]
                    - b.flux_mass[t, x]
                    * b.hot_ch.properties[t, x].enth_mass_phase["Liq"]
                )

        # Check for cold channel temperature polarization type
        if (
            self.config.cold_ch.temperature_polarization_type
            != TemperaturePolarizationType.none
        ):

            @self.Constraint(
                self.flowsheet().config.time,
                self.difference_elements,
                doc="Temperature polarization in cold channel",
            )
            def eq_temperature_polarization_cold(b, t, x):
                return (
                    b.cold_ch.h_conv[t, x]
                    * (
                        -b.cold_ch.properties[t, x].temperature
                        + b.cold_ch.properties_interface[t, x].temperature
                    )
                    == b.flux_conduction_heat[t, x]
                    + b.flux_enth_cold[t, x]
                    - b.flux_mass[t, x]
                    * b.cold_ch.properties[t, x].enth_mass_phase["Liq"]
                )

        return self.eq_flux_mass

    def _add_heat_flux(self):

        # todo: add calculation method for permeability coefficient
        self.membrane_thickness = Var(
            initialize=1e-4,
            bounds=(1e-5, 1e-2),
            doc="membrane thickness",
            units=pyunits.m,
        )

        self.membrane_tc = Var(
            initialize=0.2,
            bounds=(0, 1),
            units=pyunits.J * pyunits.s**-1 * pyunits.K**-1 * pyunits.m**-1,
            doc="Thermal conductivity coefficient of the membrane",
        )

        self.flux_conduction_heat = Var(
            self.flowsheet().config.time,
            self.difference_elements,
            initialize=1e-4,
            bounds=(0, 1e10),
            units=pyunits.J * pyunits.s**-1 * pyunits.m**-2,
            doc="conduction heat flux",
        )

        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            doc="conduction heat flux",
        )
        def eq_flux_heat(b, t, x):
            return b.flux_conduction_heat[
                t, x
            ] == b.membrane_tc / b.membrane_thickness * (
                b.hot_ch.properties_interface[t, x].temperature
                - b.cold_ch.properties_interface[t, x].temperature
            )

        @self.Expression(
            self.flowsheet().config.time,
            doc="Average conduction heat flux expression",
        )
        def flux_conduction_heat_avg(b, t):
            return (
                sum(b.flux_conduction_heat[t, x] for x in self.difference_elements)
                / self.nfe
            )

        return self.eq_flux_heat

    def _add_length_and_width(self):
        units_meta = (
            self.config.hot_ch.property_package.get_metadata().get_derived_units
        )
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

        if not hasattr(self, "area"):
            self.area = Var(
                initialize=10,
                bounds=(1e-1, 1e3),
                domain=NonNegativeReals,
                units=pyunits.m**2,
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

    def initialize_build(
        self,
        initialize_guess=None,
        state_args_hot_ch=None,
        state_args_cold_ch=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):

        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        hot_ch_flags = self.hot_ch.initialize(
            state_args=state_args_hot_ch,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            initialize_guess=initialize_guess,
            type="hot_ch",
        )

        init_log.info_high("Initialization Step 1a (hot channel) Complete")

        cold_ch_flags = self.cold_ch.initialize(
            state_args=state_args_cold_ch,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            initialize_guess=initialize_guess,
            type="cold_ch",
        )

        init_log.info_high("Initialization Step 1b (cold channel) Complete")

        if degrees_of_freedom(self) != 0:
            raise Exception(
                f"{self.name} degrees of freedom were not 0 at the beginning "
                f"of initialization. DoF = {degrees_of_freedom(self)}"
            )

        # solver
        opt = get_solver(solver, optarg)

        # Solve unit *without* any flux equations
        self.eq_flux_mass.deactivate()
        self.eq_flux_heat.deactivate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high(f"Initialization Step 2 (No Flux) {idaeslog.condition(res)}")

        # Activate only the heat flux equations
        self.eq_flux_heat.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high(
            f"Initialization Step 3 (Heat Flux Only) {idaeslog.condition(res)}"
        )

        # Activate mass flux equations as well
        self.eq_flux_mass.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high(
            f"Initialization Step 4 (Heat and Mass Flux) {idaeslog.condition(res)}"
        )

        # Release inlet state
        self.cold_ch.release_state(cold_ch_flags, outlvl)
        self.hot_ch.release_state(hot_ch_flags, outlvl)

        init_log.info(f"Initialization Complete: {idaeslog.condition(res)}")

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize")

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Hot channel Inlet": self.hot_ch_inlet,
                "Hot channel Outlet": self.hot_ch_outlet,
                "Cold channel Inlet": self.cold_ch_inlet,
                "Cold channel Outlet": self.cold_ch_outlet,
            },
            time_point=time_point,
        )

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        expr_dict = {}

        var_dict["Mass Recovery Rate"] = self.recovery_mass[time_point]
        var_dict["Membrane Area"] = self.area
        if hasattr(self, "length"):
            var_dict["Membrane Length"] = self.length
        if hasattr(self, "width"):
            var_dict["Membrane Width"] = self.width

        expr_dict["Average Solute Flux (LMH)"] = self.flux_mass_avg[time_point] * 1000
        expr_dict["Thermal efficiency (%)"] = self.thermal_efficiency[time_point] * 100
        expr_dict["Effectiveness (%)"] = self.effectiveness[time_point] * 100

        return {"vars": var_dict, "exprs": expr_dict}

    def calculate_scaling_factors(self):

        if iscale.get_scaling_factor(self.dens_solvent) is None:
            sf = iscale.get_scaling_factor(
                self.hot_ch.properties[0, self.first_element].dens_mass_phase["Liq"]
            )
            iscale.set_scaling_factor(self.dens_solvent, sf)

        super().calculate_scaling_factors()
        for t, v in self.recovery_mass.items():
            sf = 1
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)

        if iscale.get_scaling_factor(self.area) is None:
            sf = iscale.get_scaling_factor(self.area, default=10, warning=True)
            iscale.set_scaling_factor(self.area, sf)

        if iscale.get_scaling_factor(self.permeability_coef) is None:
            iscale.set_scaling_factor(self.permeability_coef, 1e10)

        if iscale.get_scaling_factor(self.membrane_thickness) is None:
            iscale.set_scaling_factor(self.membrane_thickness, 1e4)

        if iscale.get_scaling_factor(self.membrane_tc) is None:
            iscale.set_scaling_factor(self.membrane_tc, 10)

        for (t, x), v in self.flux_mass.items():
            if iscale.get_scaling_factor(v) is None:

                sf_pressure_hot = iscale.get_scaling_factor(
                    self.hot_ch.properties_interface[t, x].pressure_sat
                )

                sf_permeability = iscale.get_scaling_factor(self.permeability_coef)
                sf_thickness = iscale.get_scaling_factor(self.membrane_thickness)

                sf_flux = sf_permeability * (sf_pressure_hot) / sf_thickness

                iscale.set_scaling_factor(v, sf_flux)

                iscale.constraint_scaling_transform(self.eq_flux_mass[t, x], sf_flux)

        for (t, x), v in self.flux_enth_hot.items():
            if iscale.get_scaling_factor(v) is None:
                sf_flux_enth = sf_flux * iscale.get_scaling_factor(
                    self.hot_ch.properties_vapor[t, x].enth_mass_phase["Vap"]
                )
                iscale.set_scaling_factor(v, sf_flux_enth)

        for (t, x), v in self.flux_enth_cold.items():
            if iscale.get_scaling_factor(v) is None:
                sf_flux_enth = sf_flux * iscale.get_scaling_factor(
                    self.cold_ch.properties_vapor[t, x].enth_mass_phase["Vap"]
                )
                iscale.set_scaling_factor(v, sf_flux_enth)

        for (t, x), v in self.flux_conduction_heat.items():
            if iscale.get_scaling_factor(v) is None:
                sf_flux_cond = (
                    iscale.get_scaling_factor(self.membrane_tc)
                    / iscale.get_scaling_factor(self.membrane_thickness)
                    * iscale.get_scaling_factor(
                        self.hot_ch.properties_interface[t, x].temperature
                    )
                )
                iscale.set_scaling_factor(v, sf_flux_cond)

        if hasattr(self, "length"):
            if iscale.get_scaling_factor(self.length) is None:
                iscale.set_scaling_factor(self.length, 1)

        if hasattr(self, "width"):
            if iscale.get_scaling_factor(self.width) is None:
                iscale.set_scaling_factor(self.width, 1)
