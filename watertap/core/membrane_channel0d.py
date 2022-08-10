
from idaes.core import (declare_process_block_class, ControlVolume0DBlock,
    EnergyBalanceType,)
from watertap.core.membrane_channel_base import (
    MembraneChannelMixin, MassTransferCoefficient, PressureChangeType, ConcentrationPolarizationType) 


@declare_process_block_class("MembraneChannel0D")
class MembraneChannel0DBlockData(ControlVolume0DBlockData, MembraneChannelMixin):

    def build(self):
        super().build()

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
        # quack like a 1D model
        self.length_domain = Set(ordered=True, initialize=(0.0, 1.0))
        add_object_reference(self, "difference_elements", self.length_domain)
        self.first_element = self.length_domain.first()
        add_object_reference(
            self,
            "properties",
            {
                **{
                    (t, 0.0): self.properties_in[t]
                    for t in self.flowsheet().config.time
                },
                **{
                    (t, 1.0): self.properties_out[t]
                    for t in self.flowsheet().config.time
                },
            },
        )

        self._add_interface_blocks(information_flow, has_phase_equilibrium)

    def add_total_enthalpy_balances(self,**kwargs):
        eb = super().add_total_enthalpy_balances(**kwargs)

        if self._constructed_material_balance_type not in (EnergyBalanceType.enthalpyTotal, EnergyBalanceType.none):
            raise BalanceTypeNotSupportedError(
                "{self.name} OD membrane channels do not support {self._constructed_material_balance_type}")

        if kwargs.get("has_enthalpy_transfer", False) and eb is not None:
            # Non-existent in MembraneChannel1D 
            @self.Constraint(
                self.flowsheet().config.time, doc="Enthalpy transfer from feed to permeate"
            )
            def eq_connect_enthalpy_transfer(b, t):
                return (
                    b.mixed_permeate[t].get_enthalpy_flow_terms("Liq")
                    == -b.enthalpy_transfer[t]
                )

        return eb

    def _add_calculated_pressure_change_mass_transfer_components(self):
        if hasattr(self, "channel_height"):
            return

        # not optional in 1DRO
        self.width = Var(
            initialize=1,
            bounds=(0.1, 5e2),
            domain=NonNegativeReals,
            units=units_meta("length"),
            doc="Effective feed-channel width",
        )
        # comes from ControlVolume1D in 1DRO
        self.area = Var(
            initialize=1e-3 * 1 * 0.95,
            bounds=(0, 1e3),
            domain=NonNegativeReals,
            units=units_meta("length") ** 2,
            doc="Cross sectional area",
        )

        super()._add_calculated_pressure_change_mass_transfer_components()

    def _add_calculated_pressure_change(self):
        # already created by ControlVolume1D
        self.deltaP = Var(
            self.flowsheet().config.time,
            self.length_domain,
            initialize=-5e4,
            bounds=(-2e5, -1e3),
            domain=NegativeReals,
            units=units_meta("pressure") * units_meta("length") ** -1,
            doc="Pressure drop per unit length of feed channel at inlet and outlet",
        )
        super()._add_calculated_pressure_change()
