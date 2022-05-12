"""
Interface for flowsheet in :module:`metab`.
"""
from watertap.ui.api import FlowsheetInterface, WorkflowActions
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.metab import (
    metab,
)


def flowsheet_for_ui():
    fsi = FlowsheetInterface({"display_name": "METAB treatment train", "variables": []})
    fsi.set_action(WorkflowActions.build, build_flowsheet)
    fsi.set_action(WorkflowActions.solve, solve_flowsheet)
    # note: don't have any flowsheet block yet, will get that in build_flowsheet()
    return fsi


def build_flowsheet(fs, ui: FlowsheetInterface = None, **kwargs):
    model = metab.build()
    metab.set_operating_conditions(model)
    metab.assert_degrees_of_freedom(model, 0)
    metab.assert_units_consistent(model)
    # now we have the flowsheet block, so set it in the FlowsheetInterface object (ui)
    ui.set_block(model)


def solve_flowsheet(model, **kwargs):
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


# Terminal interactive testing
# -----------------------------

def cli_driver(fsi):
    """Dumb little interactive driver for CLI testing.
    """
    commands = {
        "json": print_json, "quit": exit_program,
        "build": run_build, "solve": run_solve,
        "results": print_results, "update": update,
    }
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


def update(fsi: FlowsheetInterface):
    """Pretend to update the values in the flowsheet"""
    data = fsi.as_dict()
    fsi.update(data)


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
        desc = cmd_val.__doc__
        desc = "" if desc is None else desc.strip()
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


def print_results(fsi: FlowsheetInterface):
    """Print the results of solving the flowsheet"""
    fs = fsi.block.fs

    print("Performance results")
    print("-------------------")
    metab.display_results(fs)

    print("Costing results")
    print("---------------")
    metab.display_costing(fs)


if __name__ == "__main__":
    fsi = flowsheet_for_ui()
    cli_driver(fsi)
