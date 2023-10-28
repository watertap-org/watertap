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

import pyomo.environ as pyo
from pyomo.common.config import ConfigBlock, ConfigValue, ListOf
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from idaes.core import declare_process_block_class, ProcessBlockData, UnitModelBlockData
from idaes.core.util.misc import add_object_reference
from idaes.core.base.costing_base import (
    UnitModelCostingBlockData,
    UnitModelCostingBlock,
    assert_flowsheet_costing_block,
    DefaultCostingComponents,
)


@declare_process_block_class("MultipleChoiceCostingBlock")
class MultipleChoiceCostingBlockData(UnitModelCostingBlockData):
    """
    Class for constructing several costing blocks on the same
    unit model and then allowing for choice between them

    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "flowsheet_costing_block",
        ConfigValue(
            domain=assert_flowsheet_costing_block,
            doc="Reference to associated FlowsheetCostingBlock to use.",
        ),
    )
    CONFIG.declare(
        "costing_methods",
        ConfigValue(
            doc="Costing methods to use for unit",
        ),  # TODO:figure out domain (ListOf(callable) doesn't work)
    )
    CONFIG.declare(
        "costing_methods_arguments",
        ConfigValue(
            default={},
            domain=dict,
            doc="Arguments to be passed to each of the costing methods. "
            "Should be a dict whose keys are the costing methods and whose "
            "values are the keyword arguments for each costing method. "
            "Costing methods without an entry will be treated as if they "
            "do not have any keyword arguments.",
        ),
    )
    CONFIG.declare(
        "initial_active_method",
        ConfigValue(
            doc="Costing method to be initially active",
            default=None,
        ),
    )

    def build(self):
        ProcessBlockData.build(self)

        # Alias flowsheet costing block reference
        fcb = self.config.flowsheet_costing_block

        # Get reference to unit model
        unit_model = self.parent_block()

        # Check that parent is an instance of a UnitModelBlockData
        if UnitModelBlockData not in unit_model.__class__.__mro__:
            raise TypeError(
                f"{self.name} - parent object ({unit_model.name}) is not an "
                f"instance of a UnitModelBlockData object. "
                "UnitModelCostingBlocks can only be added to UnitModelBlocks."
            )

        # Check to see if unit model already has costing
        for b in unit_model.component_objects(pyo.Block, descend_into=False):
            if b is not self and isinstance(
                b, (UnitModelCostingBlock, MultipleChoiceCostingBlock)
            ):
                # Block already has costing, clean up and raise exception
                raise RuntimeError(
                    f"Unit model {unit_model.name} already has a costing block"
                    f" registered: {b.name}. Each unit may only have a single "
                    "UnitModelCostingBlock associated with it."
                )
        # Add block to unit model initialization order
        unit_model._initialization_order.append(self)

        # Register unit model with this costing package
        fcb._registered_unit_costing.append(self)

        # Assign object references for costing package and unit model
        add_object_reference(self, "costing_package", fcb)
        add_object_reference(self, "unit_model", unit_model)

        self.costing_blocks = pyo.Block(self.config.costing_methods)
        self.costing_selector = pyo.Param(
            self.config.costing_methods, domain=pyo.Boolean, default=0, mutable=True
        )

        if self.config.initial_active_method is None:
            self.costing_selector[self.config.costing_methods[0]].set_value(1)
        else:
            self.costing_selector[self.config.initial_active_method].set_value(1)

        # Get costing method if not provided
        for idx, method in enumerate(self.config.costing_methods):
            blk = self.costing_blocks[method]

            # Assign object references for costing package and unit model
            add_object_reference(blk, "costing_package", fcb)
            add_object_reference(blk, "unit_model", unit_model)

            try:
                kwds = self.config.costing_methods_arguments[method]
            except KeyError:
                try:
                    kwds = self.config.costing_methods_arguments[idx]
                except KeyError:
                    kwds = {}

            # Call unit costing method
            method(blk, **kwds)

            # Check that costs are Vars and have lower bound of 0
            cost_vars = DefaultCostingComponents
            for v in cost_vars:
                try:
                    cvar = getattr(self, v)
                    if not isinstance(cvar, pyo.Var):
                        raise TypeError(
                            f"{unit_model.name} {v} component must be a Var. "
                            "Please check the costing package you are using to "
                            "ensure that all costing components are declared as "
                            "variables."
                        )
                    elif cvar.lb is None or cvar.lb < 0:
                        _log.warning(
                            f"{unit_model.name} {v} component has a lower bound "
                            "less than zero. Be aware that this may result in "
                            "negative costs during optimization."
                        )
                except AttributeError:
                    pass

        # Now we need to tie them all together
        cost_vars = list(DefaultCostingComponents) + ["direct_capital_cost"]
        for vname in cost_vars:
            for blk in self.costing_blocks.values():
                if hasattr(blk, vname):
                    break
            else:  # no break
                continue

            expr = 0.0
            for name, blk in self.costing_blocks.items():
                cvar = blk.component(vname)
                if cvar is None:
                    continue
                expr += self.costing_selector[name] * cvar

            self.add_component(vname, pyo.Expression(expr=expr))

    def select_costing_method(self, costing_method):
        """
        Set the active costing method
        """
        # zero out everything else
        self.costing_selector[:].set_value(0)
        self.costing_selector[costing_method].set_value(1)

    def initialize(self, *args, **kwargs):
        """
        See who is active, initialize that block
        ??? maybe we should just initialize them all ???
        """
        # TODO: Implement an initialization method
        # TODO: Need to have a general purpose method (block triangularisation?)
        # TODO: Should also allow registering custom methods

        # Vars and Constraints
        for blk in self.costing_blocks.values():
            for c in DefaultCostingComponents:
                if hasattr(blk, c):
                    var = getattr(blk, c)
                    cons = getattr(blk, f"{c}_constraint")
                    calculate_variable_from_constraint(var, cons)
