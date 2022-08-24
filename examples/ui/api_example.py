from watertap.ui.fsapi import FlowsheetInterface
from idaes.core.solvers import get_solver


def export_to_ui():
    return FlowsheetInterface(
        name="foo",
        do_export=export_variables,
        do_build=build_flowsheet,
        do_solve=solve_flowsheet,
    )


def export_variables(flowsheet=None, exports=None):
    for n in 1, 2:
        t = getattr(flowsheet, f"Tank{n}")
        exports.add(obj=t.inlet.flow_vol[0], name=f"Tank {n} inlet flowrate")
    # ..etc..


def build_flowsheet():
    model = example_flowsheet()
    return model.fs


def solve_flowsheet(flowsheet=None):
    solver = get_solver()
    # results = solver.solve(flowsheet)
    return {}


def example_flowsheet():
    """Example flowsheet from the IDAES tutorial."""
    # Import Pyomo libraries
    from pyomo.environ import ConcreteModel, TransformationFactory
    from pyomo.network import Arc

    # Import IDAES core
    from idaes.core import FlowsheetBlock

    # Import Unit Model Modules
    import idaes.models.properties.examples.saponification_thermo as thermo_props
    import idaes.models.properties.examples.saponification_reactions as reaction_props

    # Import Unit Model Modules
    from idaes.models.unit_models.cstr import CSTR

    m = ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # Add property packages to flowsheet library
    m.fs.thermo_params = thermo_props.SaponificationParameterBlock()
    m.fs.reaction_params = reaction_props.SaponificationReactionParameterBlock(
        default={"property_package": m.fs.thermo_params}
    )

    # Create unit models
    m.fs.Tank1 = CSTR(
        default={
            "property_package": m.fs.thermo_params,
            "reaction_package": m.fs.reaction_params,
            "has_equilibrium_reactions": False,
            "has_heat_of_reaction": True,
            "has_heat_transfer": True,
            "has_pressure_change": False,
        }
    )
    m.fs.Tank2 = CSTR(
        default={
            "property_package": m.fs.thermo_params,
            "reaction_package": m.fs.reaction_params,
            "has_equilibrium_reactions": False,
            "has_heat_of_reaction": True,
            "has_heat_transfer": True,
            "has_pressure_change": False,
        }
    )

    # Make Streams to connect units
    m.fs.stream = Arc(source=m.fs.Tank1.outlet, destination=m.fs.Tank2.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    # Set inlet and operating conditions, and some initial conditions.
    m.fs.Tank1.inlet.flow_vol[0].fix(1.0)
    m.fs.Tank1.inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
    m.fs.Tank1.inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
    m.fs.Tank1.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
    m.fs.Tank1.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
    m.fs.Tank1.inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

    m.fs.Tank1.inlet.temperature.fix(303.15)
    m.fs.Tank1.inlet.pressure.fix(101325.0)

    m.fs.Tank1.volume.fix(1.0)
    m.fs.Tank1.heat_duty.fix(0.0)

    m.fs.Tank2.volume.fix(1.0)
    m.fs.Tank2.heat_duty.fix(0.0)

    return m

    # Initialize Units
    # m.fs.Tank1.initialize()
    # m.fs.Tank2.initialize(state_args={
    #         "flow_vol": 1.0,
    #         "conc_mol_comp": {"H2O": 55388.0,
    #                           "NaOH": 100.0,
    #                           "EthylAcetate": 100.0,
    #                           "SodiumAcetate": 0.0,
    #                           "Ethanol": 0.0},
    #         "temperature": 303.15,
    #         "pressure": 101325.0})
    #
