# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    useDefault,
)

# Import Pyomo libraries
from kiwisolver import Constraint
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
    Constraint,
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
            # This is the DT term, should be on left side side of mass balance
            # dM/dt= delta_dx+mass_transfer
            # but in mass_balance this is added on right
            return -1 * self.feed_side.accumulation_mass_transfer_term[t, x, "Liq", j]

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

        # the transfer in the volume due to accumulation
        @self.feed_side.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
        )
        def eq_accumulation_mass_transfer_term(b, t, x, p, j):
            return (
                b.accumulation_mass_transfer_term[t, x, p, j]
                == (
                    b.node_mass_phase_comp[t, x, p, j]
                    - b.delta_state.node_mass_phase_comp[t, x, p, j]
                )
                / b.accumulation_time[t]
                / b.length
                * b.nfe
            )

        # mass fraction equality at each node is mass / density
        @self.feed_side.Constraint(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
        )
        def eq_node_mass_frac_equality(b, t, x, p, j):
            # H2O mass fraction is taken care off by density in
            # property package
            if j == "H2O":
                return Constraint.Skip
            else:
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
        @self.feed_side.delta_state.Constraint(
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
                / b.node_dens_mass_phase[t, x, p]
            )

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

        @self.feed_side.delta_state.Expression(
            self.flowsheet().config.time,
            self.difference_elements,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
        )
        def conc_mass_phase_comp(b, t, x, p, j):
            return (
                b.node_mass_frac_phase_comp[t, x, p, j]
                * b.node_dens_mass_phase[t, x, p]
            )

    def calculate_scaling_factors(self):

        super().calculate_scaling_factors()
        sf = iscale.get_scaling_factor(self.feed_side.accumulation_time[0])
        if sf is None:
            sf = 1 / self.feed_side.accumulation_time[0].value
            iscale.set_scaling_factor(
                self.feed_side.accumulation_time[0],
                sf,
            )
        sf = iscale.get_scaling_factor(self.feed_side.volume)
        if sf is None:
            sf = 1 / self.feed_side.volume.value
            iscale.set_scaling_factor(
                self.feed_side.volume,
                sf,
            )

        for t, x, p in self.feed_side.delta_state.node_dens_mass_phase:
            sf = iscale.get_scaling_factor(
                self.feed_side.properties[t, x].dens_mass_phase[p]
            )
            iscale.set_scaling_factor(
                self.feed_side.delta_state.node_dens_mass_phase[t, x, p],
                sf,
            )
        sf = 1
        for t, x, p, j in self.feed_side.delta_state.node_mass_frac_phase_comp:
            iscale.set_scaling_factor(
                self.feed_side.delta_state.node_mass_frac_phase_comp[t, x, p, j],
                sf,
            )
            iscale.constraint_scaling_transform(
                self.feed_side.delta_state.eq_mass_frac_phase_comp[t, x, p, j], sf
            )
            if "H2O" not in j:
                iscale.constraint_scaling_transform(
                    self.feed_side.eq_node_mass_frac_equality[t, x, p, j], sf
                )
        for t, x, p, j in self.feed_side.delta_state.node_mass_phase_comp:
            sf = iscale.get_scaling_factor(
                self.feed_side.delta_state.node_mass_phase_comp[t, x, p, j]
            )
            if sf is None:
                sf_vol = iscale.get_scaling_factor(self.feed_side.volume)
                sf_dense = iscale.get_scaling_factor(
                    self.feed_side.properties[t, x].dens_mass_phase[p]
                )
                sf_mass_frac = sf_vol * sf_dense

                iscale.set_scaling_factor(
                    self.feed_side.node_mass_phase_comp[t, x, p, j],
                    sf_mass_frac,
                )
                iscale.set_scaling_factor(
                    self.feed_side.delta_state.node_mass_phase_comp[t, x, p, j],
                    sf_mass_frac,
                )

        for t, x, p, j in self.feed_side.delta_state.node_mass_phase_comp:
            sf = iscale.get_scaling_factor(
                self.feed_side.mass_transfer_term[t, x, p, j]
            )
            sf_time = iscale.get_scaling_factor(self.feed_side.accumulation_time[t])
            sf = sf * sf_time
            iscale.set_scaling_factor(
                self.feed_side.accumulation_mass_transfer_term[t, x, p, j], sf
            )
