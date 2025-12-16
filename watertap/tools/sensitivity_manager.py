from pyomo.environ import (
    Var,
    value,
    Constraint,
    units as pyunits,
)
from idaes.core import (
    ProcessBlockData,
    declare_process_block_class,
)

import idaes.core.util.scaling as iscale
import yaml


class SensitivityType:
    equality = "equality"
    upper_bound = "upper_bound"
    lower_bound = "lower_bound"


class SensitivityData:
    """class for storing sensitivity data"""

    def __init__(
        self,
        global_variable,
        model_variables,
        multiplier_variable,
        default_value,
        sensitivity_constraints,
        existing_constraints=None,
        sensitivity_type=SensitivityType.equality,
        remove_bounds=False,
    ):
        """Constructor for SensitivityData that stores all model variables,
        global variable, multiplier variable, default values, and constraints.

        Useing provided information this will be used to build a constraint where

        model_variable = multiplier_variable * global_variable for each model variable

        global_variable value will be used as the baseline for the model variables.

        If sense is equality, model variables will be fixed to the global variable value
        if sense is upper_bound, model variables will have an upper bound of the global variable value
        if sense is lower_bound, model variables will have a lower bound of the global variable value

        Args:
            global_variable: pyomo Var, the global variable to be manipulated
            model_variables: list of pyomo Vars, the model variables to be manipulated
            multiplier_variable: pyomo Var, the multiplier variable
            default_value: float, the default value for the variable
            constraints: pyomo Constraint, the constraints to be activated/deactivated
            existing_constraints: pyomo Constraint or list of Constraints, existing constraints to be deactivated/reactivated
            sensitivity_type: SensitivityType, the type of sensitivity
            remove_bounds: bool, whether to remove bounds on model variables when activating sensitivity
        """
        self.global_variable = global_variable
        self.model_variables = model_variables
        self.multiplier_variable = multiplier_variable
        self.default_value = default_value
        self.sensitivity_constraints = sensitivity_constraints
        self.existing_constraints = existing_constraints
        self.sensitivity_type = sensitivity_type
        self.remove_bounds = remove_bounds

    def set_default_value(self):
        """Set the default value for the global variable and fix model variables"""

        self.global_variable.fix(self.default_value.value)
        for var in self.model_variables:
            if self.sensitivity_type == SensitivityType.equality:
                var.fix(self.default_value)
                print(
                    "Setting default value of " + str(self.default_value),
                    self.default_value.value,
                )
            if self.sensitivity_type == SensitivityType.upper_bound:
                var.setub(value(self.default_value.ub))
                print(
                    "Setting default upper bound to " + str(self.default_value),
                    self.default_value.ub,
                )
            if self.sensitivity_type == SensitivityType.lower_bound:
                var.setlb(value(self.default_value.lb))
                print(
                    "Setting default lower bound to " + str(self.default_value),
                    self.default_value.lb,
                )

    def deactivate_constraints(self):
        """Deactivate sensitivity constraint and fix all model variables to default value
        and reactivate existing constraint if provided"""
        self.sensitivity_constraints.deactivate()
        for var in self.model_variables:
            if self.sensitivity_type == SensitivityType.equality:
                var.fix()
        if self.existing_constraints is not None:
            if isinstance(self.existing_constraints, list):
                for c in self.existing_constraints:
                    print("Reactivating existing constraint", c.name)
                    c.activate()
            else:
                print(
                    "Reactivating existing constraint", self.existing_constraints.name
                )
                self.existing_constraints.activate()

    def activate_constraints(self):
        """Activate sensitivity constraint and unfix all model variables
        and deactivate existing constraints"""
        self.sensitivity_constraints.activate()
        for var in self.model_variables:
            var.unfix()
            print(
                "Unfixed variable " + var.name,
                var.value,
                "global var",
                self.global_variable.value,
            )
            if self.remove_bounds:
                var.setlb(None)
                var.setub(None)
                print("Removed bounds on variable " + var.name)
        if self.existing_constraints is not None:
            if isinstance(self.existing_constraints, list):
                for c in self.existing_constraints:
                    print("Deactivating existing constraint", c.name)
                    c.deactivate()
            else:
                print(
                    "Deactivating existing constraint", self.existing_constraints.name
                )
                self.existing_constraints.deactivate()

    def scale_vars_and_constraints(self):
        """scale constraints and model variables"""
        iscale.set_scaling_factor(self.multiplier_variable, 1)
        sf = iscale.get_scaling_factor(self.default_value)
        if sf is not None:
            iscale.set_scaling_factor(self.global_variable, sf)
            for i in self.sensitivity_constraints:
                print(
                    "Scaling constraint "
                    + self.sensitivity_constraints[i].name
                    + " by "
                    + str(sf)
                )
                iscale.constraint_scaling_transform(self.sensitivity_constraints[i], sf)


@declare_process_block_class("SensitivityManager")
class SensitivityManagerData(ProcessBlockData):
    """Sensitivity manager block for tracking sensitivity vars and creating global
    sensitivity constraints to manipulate variables on a flowsheets,

    This is designed to aid in performing sensetivity analysis and SVOI analysis by
    creating global sensitivity variables and constraints.

    Usage guide:
    Create Sensitivity block on model. m.fs.sensitivity = SensitivityManager()

    Register sensitivities using register_sensitivity method.

    m.fs.sensitivity.register_sensitivity(
        sensitivity_name="membrane_cost",
        model_variables=[m.fs.RO_unit[0].membrane_cost, m.fs.RO_unit[1].membrane_cost],
        default_value=m.fs.RO_unit[0].membrane_cost, # this will use the value from this pyomo var as default
        sense=SensitivityType.equality, # sets an equality constraint
        )

    The above will be used to create constraint of where each membrane cost is set to global membrance cost
    multiplied by a multiplier variable.
    This will build a membrane_cost variable on m.fs block and a membrane_cost_multiplier variable

    It will then create a constraint for each model variable in provided list (if a list is provided), other wise a single constraint
    of form
        m.fs.RO_unit[0].membrane_cost == m.fs.sensitivity.membrane_cost_multiplier * m.fs.sensitivity.membrane_cost


    m.fs.sensitivity.register_sensitivity(
        sensitivity_name="pump_pressure",
        model_variables=[m.fs.pump[0].outlet.pressure[0], m.fs.pump[1].outlet.pressure[0]],
        default_value=m.fs.pump[0].outlet.pressure[0], # this will use the value from this pyomo var as default
        existing_constraints=m.fs.global_pressure_constraint, # this will deactivate this constraint when sensitivity is activated
        sense=SensitivityType.upper_bound, # sets an upper bound constraint
        )

    The above will be used to create constraint where each pressure is set to be less then a global variable
    multiplied by a multiplier variable.
    This will build a pump_pressure variable on m.fs block and a pump_pressure_multiplier variable

    It will then create a constraint for each model variable in provided list (if a list is provided), other wise a single constraint
    of form
        m.fs.pump[0].outlet.pressure[0] <= m.fs.sensitivity.pump_pressure_multiplier * m.fs.sensitivity.pump_pressure

    Once the model is built and all sensitivities are registered, call fix_and_scale method to
    fix all global variables to default values and scale variables and constraints.

    NOTE: The model sensitivity constraints are deactivated by default to avoid
    interfering with model initialization.
    Use the activate_sensitivities method to activate all sensitivity constraints.

    The fix and scale function will fix all default global variables and multiplier variables to 1.
    WARNING: This will not change variable bounds, ensure to update those before changing multiplier or global variables.
    or enable bound removal by passing remove_bounds=True to register_sensitivity method.

    Before running sensitivity analysis use the activate_sensitivities method to activate all sensitivity constraints
    and deactivate any registered existing constraints.
    After sensitivity analysis use the deactivate_sensitivities method to deactivate all sensitivity constraints and reactivate existing constraints
    """

    CONFIG = ProcessBlockData.CONFIG()

    def build(self):
        super().build()
        self.sensitivities = {}

    def register_sensitivity(
        self,
        sensitivity_name,
        model_variables,
        default_value=None,
        sensitivity_type=SensitivityType.equality,
        existing_constraints=None,
        remove_bounds=False,
    ):
        """Register a sensitivity variable and constraint on the SensitivityManager block.
        Args:
            sensitivity: str, the name of the sensitivity
            model_variables: list of pyomo Vars, the model variables to be manipulated
            default_value: pyomo Var the default value for the variable, if None, first Var will be used for default
            sense: SensitivityType, the type of sensitivity
            existing_constraints: pyomo Constraint or list of Constraints, existing constraints to be deactivated/reactivated
            remove_bounds: bool, whether to remove bounds on model variables when activating sensitivity
        """
        if default_value is None:
            if isinstance(model_variables, list):
                default_value = model_variables[0]
            else:
                default_value = model_variables
        print(
            f"Registering sensitivity for {sensitivity_name} with default value {default_value} {model_variables}"
        )
        self.add_component(
            sensitivity_name,
            Var(initialize=1, units=default_value.get_units()),
        )
        self.find_component(sensitivity_name)
        self.add_component(
            f"{sensitivity_name}_multiplier",
            Var(initialize=1, units=pyunits.dimensionless),
        )
        self.find_component(f"{sensitivity_name}_multiplier").fix()

        if isinstance(model_variables, list) is False:
            model_variables = [model_variables]

        def indexed_constraint(m, i):
            if sensitivity_type == SensitivityType.equality:
                return model_variables[i] == self.find_component(
                    f"{sensitivity_name}_multiplier"
                ) * self.find_component(sensitivity_name)

            if sensitivity_type == SensitivityType.upper_bound:
                return model_variables[i] <= self.find_component(
                    sensitivity_name
                ) * self.find_component(f"{sensitivity_name}_multiplier")
            if sensitivity_type == SensitivityType.lower_bound:
                return model_variables[i] >= self.find_component(
                    f"{sensitivity_name}_multiplier"
                ) * self.find_component(sensitivity_name)

        self.add_component(
            f"eq_{sensitivity_name}",
            Constraint(list(range(len(model_variables))), rule=indexed_constraint),
        )
        self.sensitivities[sensitivity_name] = SensitivityData(
            self.find_component(sensitivity_name),
            model_variables,
            self.find_component(f"{sensitivity_name}_multiplier"),
            default_value,
            self.find_component(f"eq_{sensitivity_name}"),
            existing_constraints=existing_constraints,
            sensitivity_type=sensitivity_type,
            remove_bounds=remove_bounds,
        )
        self.sensitivities[sensitivity_name].deactivate_constraints()
        print("Registered sensitivity for " + sensitivity_name)

    def fix_and_scale(self):
        for sense in self.sensitivities.values():
            sense.set_default_value()
            sense.scale_vars_and_constraints()

    def activate_sensitivities(self):
        for sense in self.sensitivities.values():
            sense.activate_constraints()

    def deactivate_sensitivities(self):
        for sense in self.sensitivities.values():
            sense.deactivate_constraints()

    def display_sensitivities(self):
        print("----------------Showing Sensitivities----------------")

        for sense in self.sensitivities:
            print(
                f"{sense}: default value {self.sensitivities[sense].global_variable.value}, multiplier {self.sensitivities[sense].multiplier_variable.name} {self.sensitivities[sense].multiplier_variable.value} "
            )
        print("-----------------------------------------------------")

    def generate_multiplier_yaml_template(
        self,
        filename="sweep_template_multiplier.yaml",
        multiplier_lb=0.8,
        multiplier_ub=1.2,
    ):
        sweep_dict = {}
        for sense in self.sensitivities:
            sweep_dict[f"{sense}_multiplier"] = {
                "type": "LinearSample",
                "param": self.sensitivities[sense].multiplier_variable.name,
                "lower_limit": multiplier_lb,
                "upper_limit": multiplier_ub,
                "num_samples": 3,
            }
        with open(filename, "w") as f:
            yaml.dump(sweep_dict, f, sort_keys=False)
        print(f"Generated sweep template file: {filename}")

    def generate_absolute_yaml_template(
        self,
        filename="sweep_template_absolute.yaml",
        multiplier_lb=0.8,
        multiplier_ub=1.2,
    ):
        sweep_dict = {}
        for sense in self.sensitivities:
            sweep_dict[f"{sense}"] = {
                "type": "LinearSample",
                "param": self.sensitivities[sense].global_variable.name,
                "lower_limit": self.sensitivities[sense].global_variable.value
                * multiplier_lb,
                "upper_limit": self.sensitivities[sense].global_variable.value
                * multiplier_ub,
                "num_samples": 3,
            }
        with open(filename, "w") as f:
            yaml.dump(sweep_dict, f, sort_keys=False)
        print(f"Generated sweep template file: {filename}")

    def generate_multiplier_svoi_template(
        self, filename="sweep_template_multiplier.yaml"
    ):
        sweep_dict = {}
        sweep_dict["diff_param_loop"] = {}
        for sense in self.sensitivities:
            sweep_dict["diff_param_loop"][f"{sense}_multiplier"] = {
                "diff_mode": "percentile",
                "diff_sample_type": "UniformSample",
                "param": self.sensitivities[sense].multiplier_variable.name,
                "relative_lb": -0.01,
                "relative_ub": -0.01,
                "nominal_lb": 0.1,
                "nominal_ub": 0.2,
                "num_samples": 10,
            }
        sweep_dict["diff_param_loop"]["sweep_reference_params"] = {}
        for sense in self.sensitivities:
            sweep_dict["diff_param_loop"]["sweep_reference_params"][
                f"{sense}_multiplier"
            ] = {
                "type": "UniformSample",
                "param": self.sensitivities[sense].multiplier_variable.name,
                "lower_limit": 1,
                "upper_limit": 1,
            }
        with open(filename, "w") as f:
            yaml.dump(sweep_dict, f, sort_keys=False)
        print(f"Generated sweep template file: {filename}")

    def generate_absolute_svoi_template(
        self, filename="svoi_sweep_template_absolute.yaml"
    ):
        sweep_dict = {}
        sweep_dict["diff_param_loop"] = {}
        for sense in self.sensitivities:
            sweep_dict["diff_param_loop"][f"{sense}"] = {
                "diff_mode": "percentile",
                "diff_sample_type": "UniformSample",
                "param": self.sensitivities[sense].global_variable.name,
                "relative_lb": -0.01,
                "relative_ub": -0.01,
                "nominal_lb": self.sensitivities[sense].global_variable.value,
                "nominal_ub": self.sensitivities[sense].global_variable.value,
                "num_samples": 1,
            }
        sweep_dict["diff_param_loop"]["sweep_reference_params"] = {}
        for sense in self.sensitivities:
            sweep_dict["diff_param_loop"]["sweep_reference_params"][f"{sense}"] = {
                "type": "UniformSample",
                "param": self.sensitivities[sense].global_variable.name,
                "lower_limit": self.sensitivities[sense].global_variable.value,
                "upper_limit": self.sensitivities[sense].global_variable.value,
            }
        with open(filename, "w") as f:
            yaml.dump(sweep_dict, f, sort_keys=False)
        print(f"Generated sweep template file: {filename}")
