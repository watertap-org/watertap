"""
Interface for flowsheet in :module:`metab`.
"""
from watertap.ui.api import FlowsheetInterface, BlockInterface
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.metab import (
    metab,
)


def flowsheet_for_ui():
    model = metab.build()
    fsi = FlowsheetInterface(
        model.fs, {"display_name": "METAB treatment train", "variables": []}
    )


def cli_driver(fsi):
    """Dumb little interactive driver for CLI testing.
    """
    commands = {"print": print_json, "quit": exit_program}
    while True:
        try:
            cmd = input("Command> ")
        except EOFError:
            break
        cmd = cmd.strip().lower()
        if cmd in commands:
            commands[cmd](fsi)
        else:
            print_help(cmd, commands)
    exit_program()


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
