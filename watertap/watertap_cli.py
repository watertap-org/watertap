import click
import subprocess
import os

@click.group()
def cli():
    """WaterTAP CLI for creating project templates and more."""
    pass

@cli.command()
@click.option('--no-input', is_flag=True, help="Skip prompts and use default values from cookiecutter.json")
@click.option('--extra-context', type=str, help="Comma-separated extra context (e.g., project_name=Test,author_name=Zach)")
def project(no_input, extra_context):
    """Create a new WaterTAP project from the WaterTAP template."""
    template_path = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "project_template")
    )
    click.echo(f"Using template: {template_path}")

    cmd = ["cookiecutter", template_path]
    if no_input:
        cmd.append("--no-input")
        if extra_context:
            # Parse comma-separated values like key=value
            parts = extra_context.split(",")
            for part in parts:
                cmd.append(part.strip())

    subprocess.run(cmd)

if __name__ == '__main__':
    cli()
