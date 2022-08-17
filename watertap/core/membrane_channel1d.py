
from pyomo.common.config import Bool, ConfigValue, In
from pyomo.environ import (
    NonNegativeReals,
    Param,
    Set,
    Var,
    value,
    units as pyunits,
)
from idaes.core import declare_process_block_class, DistributedVars, FlowDirection, useDefault
from idaes.core.base.control_volume1d import ControlVolume1DBlockData
from idaes.core.util import scaling as iscale
from idaes.core.util.misc import add_object_reference
from watertap.core.membrane_channel_base import MembraneChannelMixin, PressureChangeType, CONFIG_Template as Base_CONFIG_Template

CONFIG_Template = Base_CONFIG_Template()

CONFIG_Template.declare(
    "area_definition",
    ConfigValue(
        default=DistributedVars.uniform,
        domain=In(DistributedVars),
        description="Argument for defining form of area variable",
        doc="""Argument defining whether area variable should be spatially
variant or not. **default** - DistributedVars.uniform.
**Valid values:** {
DistributedVars.uniform - area does not vary across spatial domain,
DistributedVars.variant - area can vary over the domain and is indexed
by time and space.}""",
    ),
)

CONFIG_Template.declare(
    "transformation_method",
    ConfigValue(
        default=useDefault,
        description="Discretization method to use for DAE transformation",
        doc="""Discretization method to use for DAE transformation. See Pyomo
documentation for supported transformations.""",
    ),
)

CONFIG_Template.declare(
    "transformation_scheme",
    ConfigValue(
        default=useDefault,
        description="Discretization scheme to use for DAE transformation",
        doc="""Discretization scheme to use when transforming domain. See
Pyomo documentation for supported schemes.""",
    ),
)

CONFIG_Template.declare(
    "finite_elements",
    ConfigValue(
        default=20,
        domain=int,
        description="Number of finite elements in length domain",
        doc="""Number of finite elements to use when discretizing length 
        domain (default=20)""",
    ),
)

CONFIG_Template.declare(
    "collocation_points",
    ConfigValue(
        default=5,
        domain=int,
        description="Number of collocation points per finite element",
        doc="""Number of collocation points to use per finite element when
        discretizing length domain (default=5)""",
    ),
)

@declare_process_block_class("MembraneChannel1D")
class MembraneChannel0DBlockData(MembraneChannelMixin, ControlVolume1DBlockData):

    def add_geometry(self,**kwargs):
        super().add_geometry(**kwargs)

        units_meta = self.config.property_package.get_metadata().get_derived_units

        self.width = Var(
            initialize=1,
            bounds=(1e-1, 1e3),
            domain=NonNegativeReals,
            units=units_meta("length"),
            doc="Membrane width",
        )

        self._add_area_total()
        self._add_area_total_equation()


    def add_state_blocks(
        self, information_flow=FlowDirection.forward, has_phase_equilibrium=None
    ):
        """
        This method constructs the state blocks for the
        control volume.

        Args:
            information_flow: a FlowDirection Enum indicating whether
                               information flows from inlet-to-outlet or
                               outlet-to-inlet
            has_phase_equilibrium: indicates whether equilibrium calculations
                                    will be required in state blocks
            package_arguments: dict-like object of arguments to be passed to
                                state blocks as construction arguments
        Returns:
            None
        """
        super().add_state_blocks(information_flow, has_phase_equilibrium)
        self.first_element = self.length_domain.first()
        self.difference_elements = Set(
            ordered=True,
            initialize=(x for x in self.length_domain if x != self.first_element),
        )

        self.nfe = Param(
            initialize=(len(self.difference_elements)),
            units=pyunits.dimensionless,
            doc="Number of finite elements",
        )

        self._add_interface_blocks(has_phase_equilibrium)

    def add_total_enthalpy_balances(self,**kwrags):
        # make this a no-op for MC1D
        return None

    def add_mass_transfer(self):
        self._add_recovery_rejection()
        
        units_meta = self.config.property_package.get_metadata().get_derived_units

        def mass_transfer_phase_comp_initialize(b, t, x, p, j):
            return value(
                self.properties[t, x].get_material_flow_terms("Liq", j)
                * self.recovery_mass_phase_comp[t, "Liq", j]
            )

        self.mass_transfer_phase_comp = Var(
            self.flowsheet().config.time,
            self.length_domain,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=mass_transfer_phase_comp_initialize,
            bounds=(0.0, 1e6),
            domain=NonNegativeReals,
            units=units_meta("mass")
            * units_meta("time") ** -1
            * units_meta("length") ** -1,
            doc="Mass transfer to permeate",
        )

        # ==========================================================================
        # Mass transfer term equation
        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass transfer term",
        )
        def eq_mass_transfer_term(b, t, x, p, j):
            return (
                b.mass_transfer_phase_comp[t, x, p, j]
                == -b.mass_transfer_term[t, x, p, j]
            )

        # ==========================================================================
        # Feed and permeate-side mass transfer connection --> Mp,j = Mf,transfer = Jj * W * L/n
        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass transfer from feed to permeate",
        )
        def eq_connect_mass_transfer(b, t, x, p, j):
            return (
                b.permeate_side[t, x].get_material_flow_terms(p, j)
                == -b.mass_transfer_term[t, x, p, j] * b.length / b.nfe
            )

        # ==========================================================================
        # Mass flux = feed mass transfer equation
        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass transfer term",
        )
        def eq_mass_flux_equal_mass_transfer(b, t, x, p, j):
            return (
                b.flux_mass_phase_comp[t, x, p, j] * b.width
                == -b.mass_transfer_term[t, x, p, j]
            )

        # ==========================================================================
        # Final permeate mass flow rate (of solvent and solute) --> Mp,j, final = sum(Mp,j)

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Permeate mass flow rates exiting unit",
        )
        def eq_permeate_production(b, t, p, j):
            return b.mixed_permeate[t].get_material_flow_terms(p, j) == sum(
                b.permeate_side[t, x].get_material_flow_terms(p, j)
                for x in b.difference_elements
            )


    def add_isothermal_conditions(self, **kwargs):

        super().add_isothermal_conditions(**kwargs)

        ## ==========================================================================
        # Feed-side isothermal conditions
        @self.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            doc="Isothermal assumption for feed channel",
        )
        def eq_feed_isothermal(b, t, x):
            return (
                b.properties[t, b.first_element].temperature
                == b.properties[t, x].temperature
            )

    def _add_pressure_change(self, pressure_change_type=PressureChangeType.calculated):
        add_object_reference(self, "dP_dx", self.deltaP)
        units_meta = self.config.property_package.get_metadata().get_derived_units
        self.pressure_change_total = Var(
            self.flowsheet().config.time,
            initialize=-1e5,
            bounds=(-1e6, 0),
            domain=NegativeReals,
            units=units_meta("pressure"),
            doc="Pressure drop across unit",
        )
        if pressure_change_type in (PressureChangeType.fixed_per_unit_length, PressureChangeType.calculated):
            self._add_pressure_change_equation()
        elif pressure_change_type == PressureChangeType.fixed_per_stage:
            return
        else:
            raise ConfigurationError(f"Unrecognized pressure_change_type {pressure_change_type}")

    def calculate_scaling_factors(self):
        if iscale.get_scaling_factor(self.dens_solvent) is None:
            sf = iscale.get_scaling_factor(
                self.properties[0, 0].dens_mass_phase["Liq"]
            )
            iscale.set_scaling_factor(self.dens_solvent, sf)

        super().calculate_scaling_factors()

        # these variables should have user input, if not there will be a warning
        if iscale.get_scaling_factor(self.width) is None:
            sf = iscale.get_scaling_factor(self.width, default=1, warning=True)
            iscale.set_scaling_factor(self.width, sf)

        if iscale.get_scaling_factor(self.length) is None:
            sf = iscale.get_scaling_factor(self.length, default=10, warning=True)
            iscale.set_scaling_factor(self.length, sf)

        # setting scaling factors for variables

        # will not override if the user provides the scaling factor
        ## default of 1 set by ControlVolume1D
        if iscale.get_scaling_factor(self.area) == 1:
            iscale.set_scaling_factor(self.area, 100)

        for (t, x, p, j), v in self.mass_transfer_phase_comp.items():
            sf = (
                iscale.get_scaling_factor(
                    self.properties[t, x].get_material_flow_terms(p, j)
                )
                / iscale.get_scaling_factor(self.length)
            ) * value(self.nfe)
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)
            v = self.mass_transfer_term[t, x, p, j]
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)

        if hasattr(self, "pressure_change_total"):
            for v in self.pressure_change_total.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, 1e-4)

        if hasattr(self, "dP_dx"):
            for v in self.pressure_dx.values():
                iscale.set_scaling_factor(v, 1e-5)
        else:
            for v in self.pressure_dx.values():
                iscale.set_scaling_factor(v, 1e5)
