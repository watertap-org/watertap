"""
Interface for flowsheet in :module:`metab`.
"""
from watertap.ui.api import FlowsheetInterface, WorkflowActions
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.metab import (
    metab,
)


def flowsheet_for_ui():
    model = metab.build()
    fsi = FlowsheetInterface(
        model.fs, {"display_name": "METAB treatment train", "variables": []}
    )
    fsi.set_action(WorkflowActions.build, build_flowsheet, model=model)
    fsi.set_action(WorkflowActions.solve, solve_flowsheet, model=model)
    return fsi


def build_flowsheet(fs, model=None, **kwargs):
    metab.set_operating_conditions(model)
    metab.assert_degrees_of_freedom(model, 0)
    metab.assert_units_consistent(model)


def solve_flowsheet(fs, model=None, **kwargs):
    "Solve the flowsheet"
    metab.initialize_system(model)

    results = metab.solve(model)
    metab.assert_optimal_termination(results)
    # metab.display_results(model)

    metab.add_costing(model)
    model.fs.costing.initialize()

    metab.adjust_default_parameters(model)

    metab.assert_degrees_of_freedom(model, 0)
    results = metab.solve(model)
    metab.assert_optimal_termination(results)
    #display_costing(m)


def cli_driver(fsi):
    """Dumb little interactive driver for CLI testing.
    """
    commands = {"print": print_json, "quit": exit_program,
                "build": run_build, "solve": run_solve}
    while True:
        try:
            cmd = input("Command> ")
        except EOFError:
            break
        cmd = cmd.strip().lower()
        if cmd in commands:
            print(f"Running command: {cmd}")
            commands[cmd](fsi)
        else:
            print_help(cmd, commands)
    exit_program()


def run_build(fsi: FlowsheetInterface):
    """Build the flowsheet"""
    fsi.run_action(WorkflowActions.build)


def run_solve(fsi: FlowsheetInterface):
    """Solve the flowsheet"""
    fsi.run_action(WorkflowActions.solve)


def print_help(user_cmd, commands):
    if user_cmd != "":
        print(f"'{user_cmd}' is not a command.")
    print(f"Commands:")
    for cmd_key, cmd_val in commands.items():
        desc = cmd_val.__doc__.strip()
        print(f"  {cmd_key} - {desc}")


def exit_program(*args):
    """Exit the program"""
    import sys
    print("Exit program.")
    sys.exit(0)


def print_json(fsi):
    """Print the flowsheet as JSON"""
    import json

    print("Flowsheet data")
    print("--------------")
    print(json.dumps(fsi.as_dict(), indent=2))
    print("--------------")


if __name__ == "__main__":
    fsi = flowsheet_for_ui()
    cli_driver(fsi)
