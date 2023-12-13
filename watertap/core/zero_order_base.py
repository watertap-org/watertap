#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
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
This module contains the base class for all zero order unit models.
"""

from idaes.core import UnitModelBlockData, useDefault, declare_process_block_class
from idaes.core.util.config import is_physical_parameter_block
import idaes.logger as idaeslog
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.tables import create_stream_table_dataframe

from pyomo.common.config import ConfigBlock, ConfigValue, In
import pyomo.environ as pyo


# Some more inforation about this module
__author__ = "Andrew Lee"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("ZeroOrderBase")
class ZeroOrderBaseData(UnitModelBlockData):
    """
    Standard base class for zero order unit models.

    This class contains the basic consistency checks and common methods for
    zero order type models.
    """

    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""All zero-order models are steady-state only""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Zero order models do not include holdup""",
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property  calculations,
        **default** - useDefault.
        **Valid values:** {
        **useDefault** - use default package from parent model or flowsheet,
        **PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
        and used when constructing these, **default** - None.
        **Valid values:** {see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "database",
        ConfigValue(
            description="An instance of a WaterTAP Database to use for parameters."
        ),
    )
    CONFIG.declare(
        "process_subtype",
        ConfigValue(
            description="Process subtype to use when looking up parameters from database."
        ),
    )

    def build(self):
        super().build()

        # Set a placeholder attributes
        # Placeholder for technology type string
        self._tech_type = None

        # Attribute indicating what parameters need to be fixed in the model
        self._has_recovery_removal = False
        self._fixed_perf_vars = []

        # Attributed for storing contents of reporting output
        self._stream_table_dict = {}
        self._perf_var_dict = {}

        # Check that property package meets requirements
        if self.config.property_package.phase_list != ["Liq"]:
            raise ConfigurationError(
                f"{self.name} configured with invalid property package. "
                f"Zero-order models only support property packages with a "
                f"single phase named 'Liq'."
            )
        if not hasattr(
            self.config.property_package, "solvent_set"
        ) or self.config.property_package.solvent_set != ["H2O"]:
            raise ConfigurationError(
                f"{self.name} configured with invalid property package. "
                f"Zero-order models only support property packages which "
                f"include 'H2O' as the only Solvent."
            )
        if not hasattr(self.config.property_package, "solute_set"):
            raise ConfigurationError(
                f"{self.name} configured with invalid property package. "
                f"Zero-order models require property packages to declare all "
                f"dissolved species as Solutes."
            )
        if (
            len(self.config.property_package.solute_set)
            != len(self.config.property_package.component_list) - 1
        ):
            raise ConfigurationError(
                f"{self.name} configured with invalid property package. "
                f"Zero-order models only support `H2O` as a solvent and all "
                f"other species as Solutes."
            )

    # Place holders for assigning methods
    # used to link to initization routine
    def _initialize(self, state_args, outlvl, solver, optarg):
        raise NotImplementedError()

    # used to link to scaling routine
    def _scaling(self):
        raise NotImplementedError()

    # used to provide inlet volumetric flow
    def _get_Q(self, t):
        raise NotImplementedError()

    def initialize_build(
        self, state_args=None, outlvl=idaeslog.NOTSET, solver=None, optarg=None
    ):
        """
        Passthrough initialization routine, raises NotImplementedError if
        the unit model does not have an `_initialize` function.
        """
        self._initialize(
            state_args=state_args, outlvl=outlvl, solver=solver, optarg=optarg
        )

    def calculate_scaling_factors(self):
        """
        Placeholder scaling routine, should be overloaded by derived classes
        """
        super().calculate_scaling_factors()
        self._scaling()

    def load_parameters_from_database(self, use_default_removal=False):
        """
        Method to load parameters for from database.

        Args:
            use_default_removal - (optional) indicate whether to use defined
                                  default removal fraction if no specific value
                                  defined in database

        Returns:
            None
        """
        # Get parameter dict from database
        if self._tech_type is None:
            raise NotImplementedError(
                f"{self.name} derived zero order unit model has not "
                f"implemented the _tech_type attribute. This is required "
                f"to identify the database file to load parameters from."
            )

        # Get parameter dict from database
        pdict = self.config.database.get_unit_operation_parameters(
            self._tech_type, subtype=self.config.process_subtype
        )

        if self._has_recovery_removal:
            self.set_recovery_and_removal(pdict, use_default_removal)

        for v in self._fixed_perf_vars:
            self.set_param_from_data(v, pdict)

    def set_recovery_and_removal(self, data, use_default_removal=False):
        """
        Common utility method for setting values of recovery and removal
        fractions.

        Args:
            data - dict of parameter values to use when fixing variables
            use_default_removal - (optional) indicate whether to use defined
                                  default removal fraction if no specific value
                                  defined in database

        Returns:
            None
        """
        try:
            self.set_param_from_data(self.recovery_frac_mass_H2O, data)
        except KeyError:
            if self.recovery_frac_mass_H2O[:].fixed:
                pass
            else:
                raise

        for t, j in self.removal_frac_mass_comp:
            self.set_param_from_data(
                self.removal_frac_mass_comp[t, j],
                data,
                index=j,
                use_default_removal=use_default_removal,
            )

    def set_param_from_data(
        self, parameter, data, index=None, use_default_removal=False
    ):
        """
        General method for setting parameter values from a dict of data
        returned from a database.

        Args:
            parameter - a Pyomo Var to be fixed to value from database
            data - dict of parameter values from database
            index - (optional) index to fix if parameter is an IndexedVar
            use_default_removal - (optional) indicate whether to use defined
                                  default removal fraction if no specific value
                                  defined in database

        Returns:
            None

        Raises:
            KeyError if values cannot be found for parameter in data dict

        """

        pname = parameter.parent_component().local_name

        try:
            pdata = data[pname]
        except KeyError:
            raise KeyError(
                f"{self.name} - database provided does not contain an entry "
                f"for {pname} for technology."
            )

        if index is not None:
            try:
                pdata = pdata[index]
            except KeyError:
                if pname == "removal_frac_mass_comp" and use_default_removal:
                    try:
                        pdata = data["default_removal_frac_mass_comp"]
                        index = "default"
                    except KeyError:
                        raise KeyError(
                            f"{self.name} - database provided does not "
                            f"contain an entry for {pname} with index {index} "
                            f"for technology and no default removal was "
                            f"specified."
                        )
                else:
                    raise KeyError(
                        f"{self.name} - database provided does not contain "
                        f"an entry for {pname} with index {index} for "
                        f"technology."
                    )

        try:
            val = pdata["value"]
        except KeyError:
            raise KeyError(
                f"{self.name} - no value provided for {pname} (index: "
                f"{index}) in database."
            )
        try:
            units = getattr(pyo.units, pdata["units"])
        except KeyError:
            raise KeyError(
                f"{self.name} - no units provided for {pname} (index: "
                f"{index}) in database."
            )

        parameter.fix(val * units)
        _log.info_high(f"{parameter.name} fixed to value {val} {str(units)}")

    def get_inlet_flow(self, t):
        return self._get_Q(t)

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            self._stream_table_dict, time_point=time_point
        )

    def _get_performance_contents(self, time_point=0):
        var_dict = {}

        for k, v in self._perf_var_dict.items():
            if k in ["Solute Removal", "Reaction Extent", "Rejection"]:
                for j, vd in v[time_point, :].wildcard_items():
                    var_dict[f"{k} [{j}]"] = vd
            elif v.is_indexed():
                var_dict[k] = v[time_point]
            else:
                var_dict[k] = v

        return {"vars": var_dict}

    # -------------------------------------------------------------------------
    # Unit operation costing methods
    @staticmethod
    def _get_unit_cost_method(blk):
        """
        Get a specified cost_method if one is defined in the YAML file.
        This is meant for units with different cost methods between subtypes.
        """
        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        if "cost_method" not in parameter_dict["capital_cost"]:
            raise KeyError(
                f"Costing for {blk.unit_model._tech_type} requires a cost_method argument, however "
                f"this was not defined for process sub-type {blk.unit_model.config.process_subtype}."
            )

        return parameter_dict["capital_cost"]["cost_method"]

    @staticmethod
    def _get_tech_parameters(blk, parameter_dict, subtype, param_list):
        """
        First, need to check to see if a Block with parameters for this technology
        exists.
        Second, to handle technology subtypes all parameters need to be indexed by
        subtype. We will dynamically add subtypes to the indexing set and Vars as
        required.
        """
        # Check to see in parameter Block already exists
        try:
            # Try to get parameter Block from costing package
            pblock = getattr(
                blk.config.flowsheet_costing_block, blk.unit_model._tech_type
            )
        except AttributeError:
            # Parameter Block for this technology hasn't been added yet. Create it.
            pblock = pyo.Block()

            # Add block to FlowsheetCostingBlock
            blk.config.flowsheet_costing_block.add_component(
                blk.unit_model._tech_type, pblock
            )

            # Add subtype Set to Block
            pblock.subtype_set = pyo.Set()

            # Add required Vars
            for p in param_list:
                try:
                    vobj = pyo.Var(
                        pblock.subtype_set,
                        units=getattr(
                            pyo.units, parameter_dict["capital_cost"][p]["units"]
                        ),
                    )
                    pblock.add_component(p, vobj)
                except KeyError:
                    raise KeyError(
                        "Error when trying to retrieve costing parameter "
                        "for {p}. Please check the YAML "
                        "file for this technology for errors.".format(p=p)
                    )

        # Check to see if required subtype is in subtype_set
        vlist = []
        if subtype not in pblock.subtype_set:
            # Need to add subtype and set Vars
            pblock.subtype_set.add(subtype)

            # Set vars
            for p in param_list:
                vobj = getattr(pblock, p)
                vobj[subtype].fix(
                    float(parameter_dict["capital_cost"][p]["value"])
                    * getattr(pyo.units, parameter_dict["capital_cost"][p]["units"])
                )
                vlist.append(vobj[subtype])
        else:
            for p in param_list:
                vobj = getattr(pblock, p)
                vlist.append(vobj[subtype])

        # add conditional for cases where there is only one parameter returned
        if len(vlist) == 1:
            return vlist[0]
        else:
            return tuple(x for x in vlist)

    @staticmethod
    def _general_power_law_form(
        blk, A, B, sizing_term, factor=None, number_of_parallel_units=1
    ):
        """
        General method for building power law costing expressions.
        Args:
            number_of_parallel_units (int, optional) - cost this unit as
                        number_of_parallel_units parallel units (default: 1)
        """
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = pyo.units.convert(
            A
            * (
                pyo.units.convert(sizing_term, to_units=pyo.units.dimensionless)
                / number_of_parallel_units
            )
            ** B,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        expr *= number_of_parallel_units

        blk.costing_package.add_cost_factor(blk, factor)

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

    @property
    def default_costing_method(self):
        return self.cost_power_law_flow

    @staticmethod
    def cost_power_law_flow(blk, number_of_parallel_units=1):
        """
        General method for costing equipment based on power law form. This is
        the most common costing form for zero-order models.
        CapCost = A*(F/Fref)**B
        This method also registers electricity demand as a costed flow (if
        present in the unit operation model).
        Args:
            number_of_parallel_units (int, optional) - cost this unit as
                        number_of_parallel_units parallel units (default: 1)
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B, state_ref = blk.unit_model._get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["capital_a_parameter", "capital_b_parameter", "reference_state"],
        )

        # Get state block for flow bases
        basis = parameter_dict["capital_cost"]["basis"]
        try:
            sblock = blk.unit_model.properties_in[t0]
        except AttributeError:
            # Pass-through case
            sblock = blk.unit_model.properties[t0]

        if basis == "flow_vol":
            state = sblock.flow_vol
            sizing_term = state / state_ref
        elif basis == "flow_mass":
            state = sum(sblock.flow_mass_comp[j] for j in sblock.component_list)
            sizing_term = state / state_ref
        else:
            raise ValueError(
                f"{blk.name} - unrecognized basis in parameter "
                f"declaration: {basis}."
            )

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        # Call general power law costing method
        blk.unit_model._general_power_law_form(
            blk, A, B, sizing_term, factor, number_of_parallel_units
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )
