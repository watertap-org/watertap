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

__author__ = "Alexander V. Dudchenko"

from idaes.core.util.initialization import propagate_state
from pyomo.network import Arc
from pyomo.environ import (
    Constraint,
)
import idaes.core.util.scaling as iscale
from pyomo.network import Port
import idaes.logger as idaeslog

# Set up logger
_log = idaeslog.getLogger(__name__)


class PortContainer:
    """Port container will aggregate ports and variables
    user wants to connect between unit models
    Although, the variable can be registered to the port directly, it can lead
    to issues when using UI visualization tools"""

    def __init__(self, name, port, var_dict, unit_block_reference):
        self.name = name
        self.port = port
        self.unit_block_reference = unit_block_reference
        if var_dict != None and isinstance(var_dict, dict) == False:
            raise TypeError(
                "Var dict must be a dictionary with structure {'variable name': pyomo.var}"
            )
        self.var_dict = var_dict

    def connect_to(self, inlet):
        """connect current port to provided port, registering
        the generated connection container using provided function"""
        connection = ConnectionContainer(self, inlet)
        if hasattr(self.unit_block_reference, "register_outlet_connection"):
            self.unit_block_reference.register_outlet_connection(connection)
        else:
            self.connection = connection

    def fix(self):
        """this will fix the port and all variables in the var_dict"""
        self.port.fix()
        if self.var_dict != None:
            for var, obj in self.var_dict.items():
                obj.fix()

    def unfix(self):
        """this will unfix the port and all variables in the var_dict"""
        self.port.unfix()
        if self.var_dict != None:
            for var, obj in self.var_dict.items():
                obj.unfix()


class ConnectionContainer:
    """This serves as container for any ports and equality constraints for use between unit
    models"""

    def __init__(self, outlet, inlet):
        self.registered_equality_constraints = []
        if not isinstance(outlet, (PortContainer, Port)) or not isinstance(
            inlet, (PortContainer, Port)
        ):
            raise TypeError("Provided outlet and inlet must be PortContainer objects")
        self.build_arc(outlet, inlet)
        self.build_constraints(outlet, inlet)

    def get_port(self, possible_port_object):
        """Return port for port container or port object"""
        if isinstance(possible_port_object, PortContainer):
            return possible_port_object.port
        elif isinstance(possible_port_object, Port):
            return possible_port_object
        else:
            raise TypeError("Provided object is not a PortContainer or Port")

    def get_port_unit(self, possible_port_object):
        """Return unit for port container or port object"""
        if isinstance(possible_port_object, PortContainer):
            return possible_port_object.unit_block_reference
        elif isinstance(possible_port_object, Port):
            return possible_port_object.parent_block().name
        else:
            raise TypeError("Provided object is not a PortContainer or Port")

    def build_arc(self, outlet, inlet):
        """
        builds a standard arc while naming it with outlet and inlet name, this should
        be a unique pair always (e.g. for a unit model you should only havbe single outlet-> inlet conenction)
        """
        arc = Arc(source=self.get_port(outlet), destination=self.get_port(inlet))

        def get_safe_name(name):
            return name.replace(".", "_").replace("-", "_")

        arc_name = f"{get_safe_name(outlet.name)}_to_{get_safe_name(inlet.name)}"
        outlet.unit_block_reference.add_component(
            arc_name,
            arc,
        )

        # find it, if it was not created correctly, we will get an error
        self.registered_arc = outlet.unit_block_reference.find_component(arc_name)
        if self.registered_arc is None:
            raise ValueError(f"Arc was not created correctly for {arc_name}")
        self.unit_connection = f"{self.get_port_unit(outlet)}.{get_safe_name(outlet.name)}_to_{self.get_port_unit(inlet)}.{get_safe_name(inlet.name)}"
        _log.info("Created arc connection: %s", arc_name)

    def build_constraints(self, outlet, inlet):
        """
        builds equality constraints for provided variables and scales them.
        We ensure that outlet and inlet vars are the same.
        """
        if isinstance(outlet, PortContainer) and isinstance(inlet, PortContainer):
            none_test = [d != None for d in [outlet.var_dict, inlet.var_dict]]

            if all(none_test):
                if set(outlet.var_dict) != set(inlet.var_dict):
                    raise KeyError(
                        f"Provided inlet keys: {outlet.var_dict} do not match outlet keys: {inlet.var_dict}"
                    )

                for outlet_key in outlet.var_dict:
                    # Do not create constraint if the variable is the same)
                    if outlet.var_dict[outlet_key] is inlet.var_dict[outlet_key]:
                        pass
                    else:
                        # create equality constraint between outlet and inlet var dicts
                        outlet.unit_block_reference.add_component(
                            f"eq_{outlet_key}_{outlet.name}_to_{inlet.name}",
                            Constraint(
                                expr=outlet.var_dict[outlet_key]
                                == inlet.var_dict[outlet_key]
                            ),
                        )
                        constraint = outlet.unit_block_reference.find_component(
                            f"eq_{outlet_key}_{outlet.name}_to_{inlet.name}"
                        )
                        # ensure we register it for propagation later on
                        self.registered_equality_constraints.append(
                            (
                                constraint,
                                outlet.var_dict[outlet_key],
                                inlet.var_dict[outlet_key],
                            )
                        )
                        _log.info(
                            "Created arc constraint: %s",
                            f"eq_{outlet_key}_{outlet.name}_to_{inlet.name}",
                        )
                        # scale the constraint
                        sf = iscale.get_scaling_factor(outlet.var_dict[outlet_key])
                        if sf != None:
                            iscale.constraint_scaling_transform(
                                constraint,
                                sf,
                            )

    def propagate(self):
        """this should prop any arcs and also ensure all equality constraints are satisfied"""
        _log.info("Propagating connection: %s", self.unit_connection)
        propagate_state(self.registered_arc)
        self.propagate_equality_constraints()

    def propagate_equality_constraints(self):
        """this will ensure that all equality constraints are satisfied by setting the inlet var to the outlet var"""
        for (
            constraint,
            outlet_var,
            inlet_var,
        ) in self.registered_equality_constraints:
            # set the inlet var to the outlet var value
            inlet_var.value = outlet_var.value
