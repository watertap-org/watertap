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
"""
Thickener unit model for BSM2 and plant-wide wastewater treatment modeling. 
This unit inherits from the IDAES separator unit.

Model based on

J. Alex, L. Benedetti, J.B. Copp, K.V. Gernaey, U. Jeppsson,
I. Nopens, M.N. Pons, C. Rosen, J.P. Steyer and
P. A. Vanrolleghem
Benchmark Simulation Model no. 2 (BSM2)

Modifications made to TSS formulation based on ASM type.
"""
from enum import Enum, auto

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
)
from idaes.models.unit_models.separator import SeparatorData, SplittingType
from idaes.core.util.constants import Constants
from idaes.core.util.tables import create_stream_table_dataframe
import idaes.logger as idaeslog

from pyomo.environ import (
    Constraint,
    Param,
    Var,
    NonNegativeReals,
    units as pyunits,
)
from pyomo.common.config import ConfigValue, In

from idaes.core.util.exceptions import (
    ConfigurationError,
)
from idaes.core.scaling import CustomScalerBase, ConstraintScalingScheme
from watertap.costing.unit_models.thickener import cost_thickener

__author__ = "Alejandro Garciadiego, Adam Atia"


# Set up logger
_log = idaeslog.getLogger(__name__)


class ActivatedSludgeModelType(Enum):
    """
    ASM1: ASM1 model
    ASM2D: ASM2D model
    modified_ASM2D: modified ASM2D model for ADM1 compatibility
    """

    ASM1 = auto()
    ASM2D = auto()
    modified_ASM2D = auto()


class ThickenerScaler(CustomScalerBase):
    """
    Default modular scaler for the thickener unit model.
    This Scaler relies on the associated property and reaction packages,
    either through user provided options (submodel_scalers argument) or by default
    Scalers assigned to the packages.
    """

    DEFAULT_SCALING_FACTORS = {
        "volume": 1e-3,
        "electricity_consumption": 1e1,
        "height": 1,
        "diameter": 1,
        "hydraulic_retention_time": 1e-4,
    }

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        """
        Routine to apply scaling factors to variables in model.
        Args:
            model: model to be scaled
            overwrite: whether to overwrite existing scaling factors
            submodel_scalers: dict of Scalers to use for sub-models, keyed by submodel local name
        Returns:
            None
        """
        # Call scaling methods for sub-models
        self.call_submodel_scaler_method(
            submodel=model.mixed_state,
            method="variable_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )
        self.propagate_state_scaling(
            target_state=model.underflow_state,
            source_state=model.mixed_state,
            overwrite=overwrite,
        )
        self.propagate_state_scaling(
            target_state=model.overflow_state,
            source_state=model.mixed_state,
            overwrite=overwrite,
        )

        self.call_submodel_scaler_method(
            submodel=model.underflow_state,
            method="variable_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )
        self.call_submodel_scaler_method(
            submodel=model.overflow_state,
            method="variable_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )

        # Scale unit level variables
        self.scale_variable_by_default(model.volume[0], overwrite=overwrite)
        self.scale_variable_by_default(
            model.electricity_consumption[0], overwrite=overwrite
        )
        self.scale_variable_by_default(model.diameter, overwrite=overwrite)
        self.scale_variable_by_default(model.height, overwrite=overwrite)
        self.scale_variable_by_default(
            model.hydraulic_retention_time[0], overwrite=overwrite
        )

    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        """
        Routine to apply scaling factors to constraints in model.
        Submodel Scalers are called for the property and reaction blocks. All other constraints
        are scaled using the inverse maximum scheme.
        Args:
            model: model to be scaled
            overwrite: whether to overwrite existing scaling factors
            submodel_scalers: dict of Scalers to use for sub-models, keyed by submodel local name
        Returns:
            None
        """
        # Call scaling methods for sub-models
        self.call_submodel_scaler_method(
            submodel=model.mixed_state,
            method="constraint_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )
        self.call_submodel_scaler_method(
            submodel=model.underflow_state,
            method="constraint_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )
        self.call_submodel_scaler_method(
            submodel=model.overflow_state,
            method="constraint_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )

        # Scale unit level constraints
        for c in model.component_data_objects(Constraint, descend_into=True):
            self.scale_constraint_by_nominal_value(
                c,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )


@declare_process_block_class("Thickener")
class ThickenerData(SeparatorData):
    """
    Thickener unit model for BSM2
    """

    default_scaler = ThickenerScaler

    CONFIG = SeparatorData.CONFIG()
    CONFIG.outlet_list = ["underflow", "overflow"]
    CONFIG.split_basis = SplittingType.componentFlow

    CONFIG.declare(
        "activated_sludge_model",
        ConfigValue(
            default=ActivatedSludgeModelType.ASM1,
            domain=In(ActivatedSludgeModelType),
            description="Activated Sludge Model used with unit",
            doc="""
        Options to account for version of activated sludge model property package.

        **default** - ``ActivatedSludgeModelType.ASM1``

    .. csv-table::
        :header: "Configuration Options", "Description"

        "``ActivatedSludgeModelType.ASM1``", "ASM1 model"
        "``ActivatedSludgeModelType.ASM2D``", "ASM2D model"
        "``ActivatedSludgeModelType.modified_ASM2D``", "modified ASM2D model for ADM1 compatibility"
    """,
        ),
    )

    def build(self):
        """
        Begin building model.
        Args:
            None
        Returns:
            None
        """

        # Call UnitModel.build to set up dynamics
        super(ThickenerData, self).build()

        if "underflow" and "overflow" not in self.config.outlet_list:
            raise ConfigurationError(
                "{} encountered unrecognised "
                "outlet_list. This should not "
                "occur - please use overflow "
                "and underflow as outlets.".format(self.name)
            )

        self.p_thick = Param(
            initialize=0.07,
            units=pyunits.dimensionless,
            mutable=True,
            doc="Fraction of suspended solids in the underflow",
        )

        self.TSS_rem = Param(
            initialize=0.98,
            units=pyunits.dimensionless,
            mutable=True,
            doc="Fraction of suspended solids removed",
        )

        self.electricity_consumption = Var(
            self.flowsheet().time,
            units=pyunits.kW,
            bounds=(0, None),
            doc="Electricity consumption of unit",
        )

        # 0.01255 kWh/m3 average for centrifuge thickening in relation to flow capacity
        self.energy_electric_flow_vol_inlet = Param(
            self.flowsheet().time,
            units=pyunits.kWh / (pyunits.m**3),
            initialize=0.01255,
            mutable=True,
            doc="Specific electricity intensity of unit",
        )

        @self.Constraint(self.flowsheet().time, doc="Electricity consumption equation")
        def eq_electricity_consumption(blk, t):
            return blk.electricity_consumption[t] == pyunits.convert(
                blk.energy_electric_flow_vol_inlet[t] * blk.inlet.flow_vol[t],
                to_units=pyunits.kW,
            )

        self.hydraulic_retention_time = Var(
            self.flowsheet().time,
            initialize=86400,
            domain=NonNegativeReals,
            units=pyunits.s,
            doc="Hydraulic retention time",
        )
        self.volume = Var(
            self.flowsheet().time,
            initialize=1800,
            domain=NonNegativeReals,
            units=pyunits.m**3,
            doc="Volume",
        )

        self.diameter = Var(
            initialize=10,
            domain=NonNegativeReals,
            units=pyunits.m,
            doc="Thickener diameter",
        )
        self.height = Var(
            initialize=5,
            domain=NonNegativeReals,
            units=pyunits.m,
            doc="Thickener height",
        )

        @self.Constraint(self.flowsheet().time, doc="Hydraulic retention time equation")
        def eq_hydraulic_retention(blk, t):
            return (
                self.hydraulic_retention_time[t] * self.inlet.flow_vol[t]
                == self.volume[t]
            )

        @self.Expression(doc="Surface area of circular thickener")
        def surface_area(blk):
            return self.diameter**2 * Constants.pi / 4

        @self.Constraint(self.flowsheet().time, doc="Total volume equation")
        def eq_volume(blk, t):
            return self.surface_area * self.height == self.volume[t]

        @self.Expression(self.flowsheet().time, doc="Suspended solids concentration")
        def TSS_in(blk, t):
            if blk.config.activated_sludge_model == ActivatedSludgeModelType.ASM1:
                return 0.75 * (
                    sum(
                        blk.inlet.conc_mass_comp[t, i]
                        for i in blk.config.property_package.tss_component_set
                    )
                )
            elif blk.config.activated_sludge_model == ActivatedSludgeModelType.ASM2D:
                return blk.inlet.conc_mass_comp[
                    t, blk.config.property_package.tss_component_set.first()
                ]
            elif (
                blk.config.activated_sludge_model
                == ActivatedSludgeModelType.modified_ASM2D
            ):
                return blk.mixed_state[t].TSS
            else:
                raise ConfigurationError(
                    "The activated_sludge_model was not specified properly in configuration options."
                )

        @self.Expression(self.flowsheet().time, doc="Thickening factor")
        def f_thick(blk, t):
            return blk.p_thick * (10 / (blk.TSS_in[t]))

        @self.Expression(self.flowsheet().time, doc="Remove factor")
        def f_q_du(blk, t):
            return blk.TSS_rem / (pyunits.kg / pyunits.m**3) / 100 / blk.f_thick[t]

        @self.Constraint(
            self.flowsheet().time,
            self.config.property_package.particulate_component_set,
            doc="particulate fraction",
        )
        def overflow_particulate_fraction(blk, t, i):
            return blk.split_fraction[t, "overflow", i] == 1 - blk.TSS_rem

        @self.Constraint(
            self.flowsheet().time,
            self.config.property_package.non_particulate_component_set,
            doc="soluble fraction",
        )
        def non_particulate_components(blk, t, i):
            return blk.split_fraction[t, "overflow", i] == 1 - blk.f_q_du[t]

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        expr_dict = {}
        param_dict = {}
        if hasattr(self, "split_fraction"):
            for k in self.split_fraction.keys():
                if k[0] == time_point:
                    var_dict[f"Split Fraction [{str(k[1:])}]"] = self.split_fraction[k]
        var_dict["Electricity consumption"] = self.electricity_consumption[time_point]
        param_dict["Specific electricity consumption"] = (
            self.energy_electric_flow_vol_inlet[time_point]
        )
        var_dict["Unit Volume"] = self.volume[time_point]
        var_dict["Hydraulic Retention Time"] = self.hydraulic_retention_time[time_point]
        var_dict["Unit Height"] = self.height
        var_dict["Unit Diameter"] = self.diameter
        expr_dict["Surface Area"] = self.surface_area
        return {"vars": var_dict, "params": param_dict, "exprs": expr_dict}

    def _get_stream_table_contents(self, time_point=0):
        outlet_list = self.create_outlet_list()

        io_dict = {}
        if self.config.mixed_state_block is None:
            io_dict["Inlet"] = self.mixed_state
        else:
            io_dict["Inlet"] = self.config.mixed_state_block

        for o in outlet_list:
            io_dict[o] = getattr(self, o + "_state")

        return create_stream_table_dataframe(io_dict, time_point=time_point)

    @property
    def default_costing_method(self):
        return cost_thickener
