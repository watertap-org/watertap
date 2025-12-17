# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    useDefault,
)

# Import Pyomo libraries
from pyomo.environ import (
    Var,
    NonNegativeReals,
    NegativeReals,
    PercentFraction,
    value,
    Block,
    Set,
    Expression,
    units as pyunits,
)
from watertap.unit_models.reverse_osmosis_base import (
    _add_has_full_reporting,
)
from watertap.unit_models.reverse_osmosis_1D import ReverseOsmosis1DData
from pyomo.common.config import Bool, ConfigBlock, ConfigValue, In

from watertap.unit_models.pseudo_steady_state.dead_volume_0D import DeadVolume0D

from watertap.core import (  # noqa # pylint: disable=unused-import
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    MembraneChannel1DBlock,
    PressureChangeType,
)

from idaes.core.util import scaling as iscale


@declare_process_block_class("ReverseOsmosis1DwithHoldUp")
class ReverseOsmosis1DwithHoldUpData(ReverseOsmosis1DData):
    """1D Reverse Osmosis Unit Model Class with HoldUp."""

    CONFIG = ReverseOsmosis1DData.CONFIG()
    CONFIG["delay_build"] = True

    def build(self):
        super().build()
        self._build_feed_side_state_blocks()

        units_meta = self.config.property_package.get_metadata().get_derived_units

        self.feed_side.volume = Var(
            initialize=1,
            units=units_meta("volume"),
            doc="Volume of dead space",
        )
        self.feed_side.accumulation_time = Var(
            self.flowsheet().config.time,
            initialize=1,
            units=units_meta("time"),
            doc="Time for accumulation",
        )
        self.dfe = Set(
            initialize=[
                round(i / self.config.finite_elements, 6)
                for i in range(1, self.config.finite_elements + 1)
            ],
            domain=PercentFraction,
        )
        self.dfe.display()
        self.feed_side.accumulation_mass_transfer_term = Var(
            self.flowsheet().config.time,
            self.dfe,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=0,
            units=units_meta("mass") / units_meta("length") / units_meta("time"),
            doc="Accumulation mass transfer term due to hold-up",
        )

        def get_custom_mass_transfer_term(t, x, j):
            # if self.feed_side.find_component("accumulation_mass_transfer_term") is None:
            #     self.feed_side.accumulation_mass_transfer_term = Var(
            #         self.flowsheet().config.time,
            #         self.difference_elements,
            #         self.config.property_package.phase_list,
            #         self.config.property_package.component_list,
            #         initialize=0,
            #         units=units_meta("mass")
            #         / units_meta("length")
            #         / units_meta("time"),
            #         doc="Accumulation mass transfer term due to hold-up",
            #     )
            print(t, x, j)
            return self.feed_side.accumulation_mass_transfer_term[t, x, "Liq", j]

        # @self.feed_side.Expression(0
        #     self.flowsheet().config.time,
        #     self.dfe,
        #     self.config.property_package.phase_list,
        #     self.config.property_package.component_list,
        # )
        # def eq_accumass_transfer_term(b, t, x, p, j):
        #     return b.accumulation_mass_transfer_term[t, x, p, j]

        ## Build standard ro model with custom mass transport term
        self._build_feed_side_transport_balances(
            custom_mass_transfer_term=get_custom_mass_transfer_term
        )
        self._build_permeate_side_state_blocks()
        self._build_mixed_permeate_state_blocks()
        self._build_connected_system()

        self.feed_side.node_mass_phase_comp = Var(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=1,
            units=units_meta("mass"),
            doc="Mass in each node due to hold-up",
        )

        self.feed_side.delta_state = Block()

        self.feed_side.delta_state.node_mass_phase_comp = Var(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=1,
            units=units_meta("mass"),
            doc="Prior mass in dead volume",
        )

        self.feed_side.delta_state.node_mass_frac_phase_comp = Var(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=1,
            units=pyunits.dimensionless,
            doc="Prior mass fraction in dead volume",
        )
        self.feed_side.delta_state.node_dens_mass_phase = Var(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
            initialize=1,
            units=units_meta("mass") / units_meta("volume"),
            doc="Prior density in dead volume",
        )

        ###### build hold-up related constraints ######
        self.feed_side.node_volume = Expression(
            expr=self.feed_side.volume / value(self.nfe)
        )

        # # the accumulation at each node is simply
        # # amount of mass transferred through the interface times the accumulation time
        # @self.feed_side.Constraint(
        #     self.flowsheet().config.time,
        #     self.difference_elements,
        #     self.config.property_package.phase_list,
        #     self.config.property_package.component_list,
        # )
        # def eq_node_mass_phase_comp(b, t, x, p, j):
        #     return (
        #         b.node_mass_phase_comp[t, x, p, j]
        #         == (
        #             b.material_flow_dx[0.0, 0.5, "Liq", "H2O"]
        #             + b.mass_transfer_term[t, x, p, j]
        #         )
        #         * b.length
        #         / b.nfe
        #         * b.accumulation_time[t]
        #         + b.delta_state.node_mass_phase_comp[t, x, p, j]
        #     )

        # the transfer in the volume due to accumulation
        @self.feed_side.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
        )
        def eq_accumulation_mass_transfer_term(b, t, x, p, j):
            return (
                -b.accumulation_mass_transfer_term[t, x, p, j]
                == (
                    b.node_mass_phase_comp[t, x, p, j]
                    - b.delta_state.node_mass_phase_comp[t, x, p, j]
                )
                / b.accumulation_time[t]
                / b.length
                * b.nfe
            )

        # volume at each node is mass / density
        @self.feed_side.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
        )
        def eq_node_mass_frac_equality(b, t, x, p, j):
            return b.properties[t, x].mass_frac_phase_comp[
                p, j
            ] == b.node_mass_phase_comp[t, x, p, j] / sum(
                b.node_mass_phase_comp[t, x, p, j]
                for j in self.config.property_package.component_list
            )

        # volume at each node is mass / density
        @self.feed_side.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
        )
        def eq_node_volume(b, t, x, p):
            return (
                self.feed_side.node_volume
                == sum(
                    b.node_mass_phase_comp[t, x, p, j]
                    for j in self.config.property_package.component_list
                )
                / b.properties[t, x].dens_mass_phase[p]
            )

        # volume at each node is mass / density
        # @self.feed_side.delta_state.Constraint(
        #     self.flowsheet().config.time,
        #     self.difference_elements,
        #     self.config.property_package.phase_list,
        # )
        # def eq_node_volume(b, t, x, p):
        #     return (
        #         self.feed_side.node_volume
        #         == sum(
        #             b.node_mass_phase_comp[t, x, p, j]
        #             for j in self.config.property_package.component_list
        #         )
        #         / b.node_dens_mass_phase[t, x, p]
        #     )

        @self.feed_side.delta_state.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass fractions after accumulation",
        )
        def eq_mass_frac_phase_comp(b, t, x, p, j):
            mass_sum = []
            for comp in self.config.property_package.component_list:
                mass_sum.append(b.node_mass_phase_comp[t, x, p, comp])
            return b.node_mass_frac_phase_comp[t, x, p, j] == (
                b.node_mass_phase_comp[t, x, p, j] / sum(mass_sum)
            )

        # self.feed_side.material_flow_linking_constraints.pprint()
        # self.feed_side.del_component(self.feed_side.material_balances)

        # @self.feed_side.Constraint(
        #     self.flowsheet().config.time,
        #     self.difference_elements,
        #     self.config.property_package.component_list,
        #     doc="Mass fractions after accumulation",
        # )
        # def material_balances(b, t, x, j):
        #     p = "Liq"
        #     return (
        #         b.material_flow_dx[t, x, p, j]
        #         == b.mass_transfer_term[t, x, p, j]
        #         + b.accumulation_mass_transfer_term[t, x, p, j]
        #     )

    def calculate_scaling_factors(self):

        super().calculate_scaling_factors()
        for t, x, p, j in self.feed_side.delta_state.node_mass_phase_comp:
            sf = iscale.get_scaling_factor(
                self.feed_side.mass_transfer_term[t, x, p, j]
            )
            iscale.set_scaling_factor(
                self.feed_side.accumulation_mass_transfer_term[t, x, p, j], sf
            )

        # @self.feed_side.delta_state.Expression(
        #     self.flowsheet().config.time,
        #     self.feed_side.length_domain,
        #     self.config.property_package.phase_list,
        #     self.config.property_package.component_list,
        # )
        # def conc_mass_phase_comp(b, t, x, p, j):
        #     return (
        #         b.node_mass_frac_phase_comp[t, x, p, j]
        #         * b.node_dens_mass_phase[t, x, p]
        #     )

    # def _add_feed_side_membrane_channel_and_geometry(self):
    #     # Check configuration errors
    #     self.accumulation_mass_transfer_term = Var(
    #         range(self.config.finite_elements),
    #         self.config.property_package.phase_list,
    #         self.config.property_package.component_list,
    #         initialize=0,
    #         units=self.config.property_package.get_metadata().get_derived_units("mass")
    #         / self.config.property_package.get_metadata().get_derived_units("time"),
    #         doc="Accumulation mass transfer term due to hold-up",
    #     )
    #     self._process_config()

    #     # Build 1D Membrane Channel
    #     self.feed_side = MembraneChannel1DBlock(
    #         dynamic=self.config.dynamic,
    #         has_holdup=self.config.has_holdup,
    #         area_definition=self.config.area_definition,
    #         property_package=self.config.property_package,
    #         property_package_args=self.config.property_package_args,
    #         transformation_method=self.config.transformation_method,
    #         transformation_scheme=self.config.transformation_scheme,
    #         finite_elements=self.config.finite_elements,
    #         collocation_points=self.config.collocation_points,
    #         # =self.accumulation_mass_transfer_term,
    #     )

    #     self.feed_side.add_geometry(length_var=self.length, width_var=self.width)
    #     self._add_area(include_constraint=True)

    # units_meta = self.config.property_package.get_metadata().get_derived_units

    # self.volume = Var(
    #     initialize=1,
    #     units=units_meta("volume"),
    #     doc="Volume of dead space",
    # )
    # self.accumulation_time = Var(
    #     self.flowsheet().config.time,
    #     initialize=1,
    #     units=units_meta("time"),
    #     doc="Time for accumulation",
    # )

    # self.feed_side.material_flow_linking_constraints.pprint()
    # self.feed_side.del_component(self.feed_side.material_flow_linking_constraints)

    # self.dead_volume = DeadVolume0D(
    #     self.feed_side.length_domain,
    #     property_package=self.config.property_package,
    # )
    # self.node_volume = Expression(expr=self.volume / value(self.nfe))

    # @self.Constraint(
    #     self.flowsheet().config.time,
    #     self.feed_side.length_domain,
    #     self.config.property_package.phase_list,
    # )
    # def equal_accumulation_time_constraint(b, t, x, p):
    #     return self.dead_volume[x].accumulation_time[t] == self.accumulation_time[t]

    # @self.Constraint(
    #     self.flowsheet().config.time,
    #     self.feed_side.length_domain,
    #     self.config.property_package.phase_list,
    # )
    # def equal_volume_constraint(b, t, x, p):
    #     return self.dead_volume[x].volume[t, p] == self.node_volume

    # @self.Constraint(
    #     self.flowsheet().config.time,
    #     self.feed_side.length_domain,
    #     self.config.property_package.phase_list,
    # )
    # def equal_volume_delta_state_constraint(b, t, x, p):
    #     return self.dead_volume[x].delta_state.volume[t, p] == self.node_volume

    # @self.feed_side.Constraint(
    #     self.flowsheet().config.time,
    #     self.feed_side.length_domain,
    #     self.config.property_package.phase_list,
    #     self.config.property_package.component_list,
    # )
    # def dead_volume_material_flow_linking_constraints(b, t, x, p, j):
    #     return (
    #         self.feed_side.properties[t, x].flow_mass_phase_comp[p, j]
    #         == self.dead_volume[x].inlet.flow_mass_phase_comp[t, p, j]
    #     )

    # @self.feed_side.Constraint(
    #     self.flowsheet().config.time,
    #     self.feed_side.length_domain,
    #     self.config.property_package.phase_list,
    #     self.config.property_package.component_list,
    # )
    # def material_flow_linking_constraints(b, t, x, p, j):
    #     return (
    #         self.feed_side._flow_terms[t, x, p, j]
    #         == self.dead_volume[x]
    #         .dead_volume.properties_out[t]
    #         .flow_mass_phase_comp[p, j]
    #     )

    # self.feed_side.material_flow_linking_constraints.pprint()

    # @self.feed_side.Constraint(
    #     self.flowsheet().config.time,
    #     self.feed_side.length_domain,
    # )
    # def dead_volume_isothermal_link(b, t, x):
    #     return (
    #         self.feed_side.properties[t, x].temperature
    #         == self.dead_volume[x].inlet.temperature[t]
    #     )

    # @self.feed_side.Constraint(
    #     self.flowsheet().config.time,
    #     self.feed_side.length_domain,
    # )
    # def dead_volume_pressure_link(b, t, x):
    #     return (
    #         self.feed_side.properties[t, x].pressure
    #         == self.dead_volume[x].inlet.pressure[t]
    #     )
