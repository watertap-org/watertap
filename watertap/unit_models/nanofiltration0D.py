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
Nanofiltration unit model with assumed 100% rejection of divalent ions.
"""

from pyomo.environ import Block, Constraint, Set, Var
from pyomo.common.config import ConfigDict, ConfigValue, In, Bool
from pyomo.network import Port

from idaes.core import (
    declare_process_block_class,
    MaterialBalanceType,
    MomentumBalanceType,
    useDefault,
    UnitModelBlockData,
)
from idaes.core.initialization import ModularInitializerBase
from idaes.core.scaling import CustomScalerBase, ConstraintScalingScheme
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
from idaes.models.unit_models.separator import EnergySplittingType
import idaes.logger as idaeslog

from watertap.flowsheets.generic_desalination_train.unit_operations.separator import initialize

_log = idaeslog.getLogger(__name__)


EPS = 1e-10  # Default value for multivalent ion rejection


class Nanofiltration0DInitializer(ModularInitializerBase):
    """
    Initializer for 0D Nanofiltration models.

    """
    CONFIG = ModularInitializerBase.CONFIG()

    def initialize_main_model(
        self,
        model: Block,
    ):
        """
        Initialization routine for main Nanofiltration0D models.

        Args:
            model: Pyomo Block to be initialized.
            copy_inlet_state: bool (default=False). Whether to copy inlet state to other states or not.
                Copying will generally be faster, but inlet states may not contain
                all properties required elsewhere.

        Returns:
            Pyomo solver results object.

        """
        # Get loggers
        init_log = idaeslog.getInitLogger(
            model.name, self.get_output_level(), tag="unit"
        )
        solve_log = idaeslog.getSolveLogger(
            model.name, self.get_output_level(), tag="unit"
        )

        # Create solver
        solver = self._get_solver()

        # ---------------------------------------------------------------------
        # Initialize inlet state
        prop_init = self.get_submodel_initializer(model.properties_in)

        if prop_init is not None:
            prop_init.initialize(
                model=model.properties_in,
                output_level=self.get_output_level(),
            )
        else:
            raise ValueError(
                "No Initializer found for property package. Please set provide "
                "sub-model initializers or assign a default Initializer."
            )

        init_log.info_high("Initialization Step 1 (Inlet State) Complete.")

        # ---------------------------------------------------------------------
        # Estimate outlet states from inlet and initialize
        for t, in_state in model.properties_in.items():
            state_vars = in_state.define_state_vars()
            state_vars_ret = model.properties_retentate[t].define_state_vars()
            state_vars_per = model.properties_permeate[t].define_state_vars()

            for n, sv in state_vars.items():
                sv_ret = state_vars_ret[n]
                sv_per = state_vars_per[n]

                if sv.local_name == "pressure" and not sv_ret.fixed:
                    # For pressure, only retentate is linked to inlet
                    if hasattr(model, "deltaP"):
                        sv_ret.set_value(sv + model.deltaP)
                    else:
                        sv_ret.set_value(sv)
                elif "flow" in sv.local_name:
                    self._init_flow(in_state, sv, sv_ret, sv_per)
                elif any(sv.startswith(i) for i in ["mass_frac", "mole_frac"]):
                    self._init_frac(in_state, sv, sv_ret, sv_per)
                elif sv.startswith("conc"):
                    self._init_conc(in_state, sv, sv_ret, sv_per)
                else:
                    # For everything else, assume outlet similar to inlet
                    for k, svd in sv.items():
                        if not sv_ret[k].fixed:
                            sv_ret[k].set_value(svd)
                        if not sv_per[k].fixed:
                            sv_per[k].set_value(svd)

        prop_init.initialize(
            model=model.properties_retentate,
            output_level=self.get_output_level(),
        )
        prop_init.initialize(
            model=model.properties_permeate,
            output_level=self.get_output_level(),
        )

        init_log.info_high("Initialization Step 2 (Outlet States) Complete.")

        # ---------------------------------------------------------------------
        # Solve full model
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver.solve(model, tee=slc.tee)

        init_log.info("Initialization Completed, {}".format(idaeslog.condition(res)))

        return res

    def _init_conc(self, in_state, sv, sv_ret, sv_per):
        # Component indexed flow - need to apply split fractions
        for k, svd in sv.items():
            # Component is always the last index
            j = k[-1]

            # Determine split based on type
            comp = in_state.params.get_component(j)

            if comp.is_solvent():
                # Assume solvent concentrations remain roughly constants
                split = 1
            else:
                try:
                    charge = comp.config.charge
                except AttributeError:
                    charge = 0

                if abs(charge) > 1:
                    split = model.multivalent_recovery
                else:
                    split = model.solute_recovery[j]

            if not sv_ret[k].fixed:
                sv_ret[k].set_value(svd * (1 - split))
            if not sv_per[k].fixed:
                sv_per[k].set_value(svd * split)

    def _init_flow(self, in_state, sv, sv_ret, sv_per):
        # Determine if it is a component flow or not
        if sv.local_name.endswith("_comp"):
            # Component indexed flow - need to apply split fractions
            for k, svd in sv.items():
                # Component is always the last index
                j = k[-1]

                # Determine split based on type
                comp = in_state.params.get_component(j)

                if comp.is_solvent():
                    split = model.solvent_recovery
                else:
                    try:
                        charge = comp.config.charge
                    except AttributeError:
                        charge = 0

                    if abs(charge) > 1:
                        split = model.multivalent_recovery
                    else:
                        split = model.solute_recovery[j]

                if not sv_ret[k].fixed:
                    sv_ret[k].set_value(svd * (1 - split))
                if not sv_per[k].fixed:
                    sv_per[k].set_value(svd * split)
        else:
            # Assume a total flow basis, and use solvent recovery
            split = model.solvent_recovery
            for k, svd in sv.items():
                if not sv_ret[k].fixed:
                    sv_ret[k].set_value(svd * (1 - split))
                if not sv_per[k].fixed:
                    sv_per[k].set_value(svd * split)

    def _init_frac(self, in_state, sv, sv_ret, sv_per):
        # First need to iterate over all indices to collect normalized flows
        # Assume a basis of 1 mass or mole unit
        # Also, only single phase property packages supported, so we only need to
        # worry about the component index
        nom_ret_comp_flow = {}
        nom_per_comp_flow = {}

        for k, svd in sv.items():
            # Component is always the last index
            j = k[-1]

            # Determine split based on type
            comp = in_state.params.get_component(j)

            if comp.is_solvent():
                split = model.solvent_recovery
            else:
                try:
                    charge = comp.config.charge
                except AttributeError:
                    charge = 0

                if abs(charge) > 1:
                    split = model.multivalent_recovery
                else:
                    split = model.solute_recovery[j]

            nom_ret_comp_flow[j] = value(svd * (1 - split))
            nom_per_comp_flow[j] = value(svd * (1 - split))

        nom_ret_tot_flow = sum(f for f in nom_ret_comp_flow.values())
        nom_per_tot_flow = sum(f for f in nom_per_comp_flow.values())

        # # Next set value based on fraction of nominal total flow
        for k, svd in sv.items():
            # Component is always the last index
            j = k[-1]

            if not sv_ret[k].fixed:
                sv_ret[k].set_value(nom_ret_comp_flow[j]/nom_ret_tot_flow)
            if not sv_per[k].fixed:
                sv_per[k].set_value(nom_per_comp_flow[j]/nom_per_tot_flow)


class Nanofiltration0DScaler(CustomScalerBase):
    DEFAULT_SCALING_FACTORS = {
        "multivalent_recovery": 1e2,
        "solvent_recovery": 10,
        "solute_recovery": 10,
    }

    def variable_scaling_routine(
            self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # Scale inlet state
        self.call_submodel_scaler_method(
            submodel=model.properties_in,
            method="variable_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )

        # Propagate scaling from inlet state
        self.propagate_state_scaling(
            target_state=model.properties_retentate,
            source_state=model.properties_in,
            overwrite=overwrite,
        )
        self.propagate_state_scaling(
            target_state=model.properties_permeate,
            source_state=model.properties_in,
            overwrite=overwrite,
        )

        # Scale remaining states
        self.call_submodel_scaler_method(
            submodel=model.properties_retentate,
            method="variable_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )
        self.call_submodel_scaler_method(
            submodel=model.properties_permeate,
            method="variable_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )

        # Scale remaining variables
        # Variables with default scaling factors
        for n in self.DEFAULT_SCALING_FACTORS.keys():
            c = model.find_component(n)
            for v in c.values():
                self.scale_variable_by_default(v, overwrite=overwrite)

    def constraint_scaling_routine(
            self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # Call scaling methods for sub-models
        self.call_submodel_scaler_method(
            submodel=model.properties_in,
            method="constraint_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )
        self.call_submodel_scaler_method(
            submodel=model.properties_retentate,
            method="constraint_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )
        self.call_submodel_scaler_method(
            submodel=model.properties_permeate,
            method="constraint_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )

        # Scale remaining constraints
        for c in model.component_data_objects(Constraint, descend_into=False):
            self.scale_constraint_by_nominal_value(
                c,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )


@declare_process_block_class("Nanofiltration0D")
class Nanofiltration0DData(UnitModelBlockData):
    """
    Zero order nanofiltration model.
    """

    default_initializer = Nanofiltration0DInitializer
    default_scaler = Nanofiltration0DScaler

    CONFIG = ConfigDict()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False. Product blocks are always steady-state.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Product blocks do not contain holdup, thus this must be
    False.""",
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for mixer",
            doc="""Property parameter object used to define property
    calculations,
    **default** - useDefault.
    **Valid values:** {
    **useDefault** - use default package from parent model or flowsheet,
    **PropertyParameterObject** - a PropertyParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigDict(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigDict with arguments to be passed to a property
    block(s) and used when constructing these,
    **default** - None.
    **Valid values:** {
    see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "include_pressure_balance",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Whether to include pressure balance for retentate side.",
        )
    )
    CONFIG.declare(
        "has_retentate_pressure_drop",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Whether to include pressure drop in the retentate side.",
        )
    )
    CONFIG.declare(
        "include_temperature_equality",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Whether to include temperature equality constraints",
            doc="""Argument indicating whether temperature equality should be
            inlcluded. If this is False, no energy balance constraint will be
            written.""",
        ),
    )
    CONFIG.declare(
        "electroneutrality_ion",
        ConfigValue(
            default="Cl_-",
            domain=str,
            description="Balancing ion to use to maintain electroneutrality",
            doc="""Name of ion to use to maintain electroneutrality in
            permeate stream. If an ion is provided, the split fraction
            for this ion will be solved implicitly to ensure electroneutrality
            is maintained in the permeate stream. If None, user must provide
            spl;it fractions for all monovalent ions and electroneutrality is
            not guaranteed..""",
        ),
    )

    def build(self):
        super().build()

        prop_params = self.config.property_package

        # Check that property package only supports single phase
        # TODO: Extend the model to handle phase equilibrium
        if len(prop_params.phase_list) > 1:
            raise ConfigurationError(
                "Nanofiltration0D model only supports single phase "
                "property packages."
            )

        # Check that the electroneutrality ion is a valid monovalent ion
        if self.config.electroneutrality_ion is not None:
            try:
                bal_ion = prop_params.get_component(self.config.electroneutrality_ion)
                # Check that balancing ion is monovalent
                if not hasattr(bal_ion.config, "charge") or abs(bal_ion.config.charge) > 1:
                    raise ConfigurationError(
                        f"electroneutrality_ion ({self.config.electroneutrality_ion}) "
                        "must be a monovalent ion."
                    )
            except AttributeError:
                raise ConfigurationError(
                    f"electroneutrality_ion ({self.config.electroneutrality_ion}) "
                    "is not a valid component in property package."
                )
        # Check that pressure balance arguments are consistent
        if self.config.has_retentate_pressure_drop and not self.config.include_pressure_balance:
            raise ConfigurationError(
                "Inconsistent configuration arguments. has_retentate_pressure_drop=True "
                "requires that include_pressure_balance=True."
            )

        tmp_dict_in = dict(**self.config.property_package_args)
        tmp_dict_in["has_phase_equilibrium"] = False
        tmp_dict_in["defined_state"] = True

        self.properties_in = self.config.property_package.build_state_block(
            self.flowsheet().time, doc="Material properties at inlet", **tmp_dict_in
        )

        tmp_dict_out = dict(**self.config.property_package_args)
        tmp_dict_out["has_phase_equilibrium"] = False
        tmp_dict_out["defined_state"] = False

        self.properties_retentate = self.config.property_package.build_state_block(
            self.flowsheet().time, doc="Material properties at inlet", **tmp_dict_out
        )

        self.properties_permeate = self.config.property_package.build_state_block(
            self.flowsheet().time, doc="Material properties at inlet", **tmp_dict_out
        )

        self.add_port("inlet", self.properties_in, doc="Inlet Port")
        self.add_port("retentate", self.properties_retentate, doc="Retentate Port")
        self.add_port("permeate", self.properties_permeate, doc="Permeate Port")

        # NF separation variables
        self.solvent_recovery = Var(initialize=0.8)
        self.multivalent_recovery = Var(initialize=EPS)
        self.multivalent_recovery.fix()

        other_solutes = []
        for j in prop_params.component_list:
            comp = prop_params.get_component(j)
            if j in prop_params.ion_set:
                # Check charge
                if j == self.config.electroneutrality_ion:
                    # This is the ion ot use to balance electronegativity - leave out of set
                    continue
                if abs(comp.config.charge) == 1:
                    # Monovalent - include in set
                    other_solutes.append(j)

                    # TODO: check for balance species
            elif not comp.is_solvent():
                # No ionic solute - include in set
                other_solutes.append(j)

        self.split_species = Set(initialize=other_solutes)
        self.solute_recovery = Var(self.split_species, initialize=0.9)

        units = self.config.property_package.get_metadata().get_derived_units
        if self.config.has_retentate_pressure_drop:
            self.deltaP = Var(self.flowsheet().time, initialize=0, units=units.PRESSURE, doc="Retentate side pressure drop")

        # Material balance
        @self.Constraint(self.flowsheet().time, prop_params.component_list)
        def material_balances(b, t, j):
            pset = b.properties_in[t].phase_list
            pcset = b.properties_in[t].phase_component_set
            return (
                sum(b.properties_in[t].get_material_flow_terms(p, j) for p in pset if (p, j) in pcset)
                == sum(b.properties_retentate[t].get_material_flow_terms(p, j) for p in pset if (p, j) in pcset)
                + sum(b.properties_permeate[t].get_material_flow_terms(p, j) for p in pset if (p, j) in pcset)
            )

        @self.Constraint(self.flowsheet().time, prop_params.component_list)
        def recovery_constraint(b, t, j):
            pset = b.properties_in[t].phase_list
            pcset = b.properties_in[t].phase_component_set
            comp = prop_params.get_component(j)

            if j in b.split_species:
                rec = b.solute_recovery[j]
            elif comp.is_solvent():
                rec =  b.solvent_recovery
            elif j == self.config.electroneutrality_ion:
                return Constraint.Skip
            else:
                rec = b.multivalent_recovery

            # Permeate flows equal to recovery * inlet flow
            return (
                rec * sum(b.properties_in[t].get_material_flow_terms(p, j) for p in pset if (p, j) in pcset)
                == sum(b.properties_permeate[t].get_material_flow_terms(p, j) for p in pset if (p, j) in pcset)
            )

        if self.config.electroneutrality_ion is not None:
            # Add electroneutrality constraint for permeate stream
            @self.Constraint(self.flowsheet().time)
            def permeate_electronegativity(b, t):
                perm_state = b.properties_permeate[t]
                return 0 == sum(
                    perm_state.params.get_component(j).config.charge
                    * perm_state.get_material_flow_terms(p, j)
                    for p, j in perm_state.phase_component_set
                    if j in perm_state.params.ion_set
                )

        # Retentate pressure balance
        if self.config.include_pressure_balance:
            @self.Constraint(self.flowsheet().time, doc="Retentate pressure balance")
            def retentate_pressure_balance(b, t):
                expr = b.properties_retentate[t].pressure
                if b.config.has_retentate_pressure_drop:
                    expr += -b.deltaP[t]
                return b.properties_in[t].pressure == expr

        # Temperature equalities
        if self.config.include_temperature_equality:
            @self.Constraint(self.flowsheet().time, doc="Retentate temperature equality")
            def retentate_temperature_equality(b, t):
                return b.properties_in[t].temperature == b.properties_retentate[t].temperature

            @self.Constraint(self.flowsheet().time, doc="Permeate temperature equality")
            def permeate_temperature_equality(b, t):
                return b.properties_in[t].temperature == b.properties_permeate[t].temperature
