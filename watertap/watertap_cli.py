import click
import subprocess
import os


@click.group()
def cli():
    """WaterTAP CLI for creating project templates and more."""


@cli.command()
@click.option("--no-input", is_flag=True, help="Skip prompts and use defaults")
@click.option(
    "--extra-context",
    type=str,
    help="Comma-separated context: project_name=...,author_name=...",
)
def project(no_input, extra_context):
    """Create a new WaterTAP project from GitHub (fallback to local if offline)"""
    GITHUB_TEMPLATE = "https://github.com/watertap-org/watertap_project_template"
    GITHUB_TEMPLATE_NAME = "watertap_project_template"
    CACHE_DIR = os.path.expanduser(f"~/.cookiecutters/{GITHUB_TEMPLATE_NAME}")

    local_template = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "project_template")
    )

    # Always remove the cached template
    if os.path.exists(CACHE_DIR):
        import shutil

        shutil.rmtree(CACHE_DIR)

    def build_cmd(template_path):
        cmd = ["cookiecutter", "--overwrite-if-exists", template_path]
        if no_input:
            cmd.append("--no-input")
            if extra_context:
                for part in extra_context.split(","):
                    cmd.append(part.strip())
        return cmd

    try:
        click.echo(f"Trying GitHub template: {GITHUB_TEMPLATE}")
        subprocess.run(build_cmd(GITHUB_TEMPLATE), check=True)
    except subprocess.CalledProcessError:
        click.echo("⚠️  GitHub template failed. Falling back to local template.")
        subprocess.run(build_cmd(local_template))


if __name__ == "__main__":
    cli()
