
from pyomo.common.config import Bool, ConfigValue, In
from pyomo.environ import (
    NonNegativeReals,
    NegativeReals,
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



@declare_process_block_class("MembraneChannel1DBlock")
class MembraneChannel1DBlockData(MembraneChannelMixin, ControlVolume1DBlockData):

    def add_geometry(self, length_var, width_var,
            flow_direction=FlowDirection.forward, **kwargs):
        """
        Method to create spatial domain and volume Var in ControlVolume.

        Args:
            length_var - external variable to use for the length of
                         the spatial domain. If a variable is provided, a
                         reference will be made to this in place of the length
                         Var.
            width_var - external variable to use for the width of
                         the spatial domain. If a variable is provided, a
                         reference will be made to this in place of the width 
                         Var.
            flow_direction - argument indicating direction of material flow
                            relative to length domain. Valid values:
                                - FlowDirection.forward (default), flow goes
                                  from 0 to 1.
                                - FlowDirection.backward, flow goes from 1 to 0
            length_domain - (optional) a ContinuousSet to use as the length
                            domain for the ControlVolume. If not provided, a
                            new ContinuousSet will be created (default=None).
                            ContinuousSet should be normalized to run between
                            0 and 1.
            length_domain_set - (optional) list of point to use to initialize
                            a new ContinuousSet if length_domain is not
                            provided (default = [0.0, 1.0]).

        Returns:
            None
        """
        super().add_geometry(length_var=length_var, flow_direction=flow_direction, **kwargs)
            
        # ControlVolume1D.add_geometry will add length reference
        self._add_reference("width", width_var)


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

        self._add_interface_stateblock(has_phase_equilibrium)

    def add_total_enthalpy_balances(self,**kwrags):
        # make this a no-op for MC1D
        return None

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
        super().calculate_scaling_factors()

        # setting scaling factors for variables

        # will not override if the user provides the scaling factor
        ## default of 1 set by ControlVolume1D
        if iscale.get_scaling_factor(self.area) == 1:
            iscale.set_scaling_factor(self.area, 100)

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
