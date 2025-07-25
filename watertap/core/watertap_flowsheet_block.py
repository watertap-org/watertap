from idaes.core.base.flowsheet_model import FlowsheetBlockData

from idaes.core import (
    declare_process_block_class,
)
from pyomo.common.config import ConfigValue

from watertap.core.utils.connections import (
    PortContainer,
    ConnectionContainer,
)
from watertap.core.utils.report import (
    build_report_table,
)


@declare_process_block_class("WaterTapFlowsheetBlock")
class WaterTapFlowsheetBlockData(FlowsheetBlockData):
    CONFIG = FlowsheetBlockData.CONFIG()
    CONFIG.declare(
        "default_costing_package",
        ConfigValue(
            default=None,
            description="defines default costing package",
            doc="""
                defines default costing package
            """,
        ),
    )
    CONFIG.declare(
        "default_costing_package_kwargs",
        ConfigValue(
            default={},
            description="kwargs to pass into the Costing unit block",
            doc="""
                kwargs to pass into the Costing unit block,
            """,
        ),
    )

    def build(self):
        self.outlet_connections = []
        super().build()

    def fix_and_scale(self):
        self.set_fixed_operation()
        self.scale_before_initialization()

    def initialize_unit(self, **kwargs):
        """Developer should implement an initialize routine for their flowsheet model"""
        pass

    def initialize(self, **kwargs):
        """routine to initialize a unit and propagate its connections"""
        self.initialize_unit()
        self.propagate_outlets()

    def propagate_outlets(self):
        """propagates registered outlet connections"""
        for outlet in self.outlet_connections:
            outlet.propagate()

    def set_fixed_operation(self, **kwargs):
        """Developer should implement a routine to fix unit operation for initialization and 0DOF solving"""
        pass

    def set_optimization_operation(self, **kwargs):
        """Developer should implement a routine to unfix variables for unit to perform optimization"""
        pass

    def setup_optimization(self, **kwargs):
        """Developer should implement a routine to unfix unit operation for specific operation"""
        pass

    def scale_before_initialization(self, **kwargs):
        """Developer should implement scaling function to scale unit using
        default values before initialization routine is ran"""
        pass

    def scale_post_initialization(self, **kwargs):
        """Developer should implement scaling function to scale unit after initialization routine is ran"""
        pass

    def register_port(self, name, port=None, var_list=None):
        """Registers a port for the flowsheet unit, including variables that should
        be connected through equality constraints

        Args:
        name - Name of the port
        port - port from a unit model, can be none
        var_list - list of variables that should be connected through equality constraints
        """
        # create our var on the block so we can reference port and variables
        setattr(
            self,
            name,
            PortContainer(name, port, var_list, self),
        )

    def register_outlet_connection(self, connection):
        """registers outlet connections to enable automatic propagation"""
        if isinstance(connection, ConnectionContainer):
            self.outlet_connections.append(connection)
        else:
            raise TypeError("Outlet connection must be a ConnectionContainer")

    def get_model_state_dict(self):
        """Developer should define model_vars dict that includes
                a list of varaibles that should be displayed
                example:

                def get_model_state_dict(self):
                    model_state = {
                        "Composition": {},
                        "Physical state": {},
                    }
                    for phase, ion in self.feed.properties[0].conc_mass_phase_comp:
                        model_state["Composition"][ion] = self.feed.properties[
                            0
                        ].conc_mass_phase_comp[phase, ion]
                    model_state["Physical state"]["pH"] = self.feed.pH
                    model_state["Physical state"]["Temperature"] = self.feed.properties[
                        0
                    ].temperature
                    model_state["Physical state"]["Pressure"] = self.feed.properties[0].pressure
                    return self.name, model_state

                Output should look like:

        ------------------------------------------------------------------------------------
            fs.feed state

            Composition:
            Key    : Value     : Units                 : Fixed : Bounds
             Ca_2+ :   0.25800 : kilogram / meter ** 3 :  True : (0, 2000.0)
              Cl_- :   0.87000 : kilogram / meter ** 3 :  True : (0, 2000.0)
               H2O :    996.64 : kilogram / meter ** 3 : False : (0, 2000.0)
            HCO3_- :   0.38500 : kilogram / meter ** 3 :  True : (0, 2000.0)
               K_+ : 0.0090000 : kilogram / meter ** 3 :  True : (0, 2000.0)
             Mg_2+ :  0.090000 : kilogram / meter ** 3 :  True : (0, 2000.0)
              Na_+ :   0.73900 : kilogram / meter ** 3 :  True : (0, 2000.0)
            SO4_2- :    1.0110 : kilogram / meter ** 3 :  True : (0, 2000.0)

            Physical state:
            Key         : Value      : Units         : Fixed : Bounds
               Pressure : 1.0000e+05 :        pascal :  True : (100000.0, None)
            Temperature :     293.15 :        kelvin :  True : (273.15, 373.15)
                     pH :     7.0700 : dimensionless :  True : (None, None)

        ------------------------------------------------------------------------------------
        """
        return None, None

    def _get_stream_table_contents(self, time_point=0):
        """override default as developer should manually define model state"""
        if self.get_model_state_dict()[0] is None:
            return super()._get_stream_table_contents(time_point)
        else:
            return None

    def report(
        self, time_point=0, dof=False, ostream=None, prefix="", use_default_units=False
    ):
        unit_name, defined_variables = self.get_model_state_dict()

        if unit_name is not None:
            build_report_table(
                unit_name, defined_variables, ostream, prefix, use_default_units
            )
        if unit_name is None:
            super().report(time_point, dof, ostream, prefix)
