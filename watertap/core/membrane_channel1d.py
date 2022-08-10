
from idaes.core import declare_process_block_class, ControlVolume1DBlock
from watertap.core.membrane_channel_base import MembraneChannelMixin


@declare_process_block_class("MembraneChannel1D")
class MembraneChannel0DBlockData(ControlVolume1DBlockData, MembraneChannelMixin):

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
        self.first_element = self.length_domain.first()
        self.difference_elements = Set(
            ordered=True,
            initialize=(x for x in self.length_domain if x != self.first_element),
        )

        self._add_interface_blocks(information_flow, has_phase_equilibrium)
