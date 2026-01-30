from pyomo.environ import (
    Var,
    units as pyunits,
)

from pyomo.common.config import ConfigBlock, ConfigValue, In

from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
)
from watertap.costing.unit_models.conduit import cost_conduit
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog


__author__ = "Alexander V. Dudchenko"


@declare_process_block_class("Conduit")
class ConduitData(UnitModelBlockData):
    """
    Trivial conduit model for flushing operation in CCRO process.
    Simply used to track volume and cost.
    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Solar energy surrogate models are steady-state only""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Solar energy models do not include holdup""",
        ),
    )

    def build(self):

        super().build()

        self.volume = Var(
            initialize=1,
            units=pyunits.m**3,
            doc="Volume of the conduit representing the accumulation volume",
        )
        self.volume.fix()

    def initialize_build(
        self,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        pass

    def calculate_scaling_factors(self):
        iscale.set_scaling_factor(self.volume, 1)

    @property
    def default_costing_method(self):
        return cost_conduit
