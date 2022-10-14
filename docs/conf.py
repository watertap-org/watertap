# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import sphinx_rtd_theme

sys.path.insert(0, os.path.abspath(".."))
# sys.path.insert(0,os.path.dirname(sys.path[0]))

# -- Project information -----------------------------------------------------

project = "WaterTAP"
copyright = "2021, NAWI"
author = "NAWI"

# The full version, including alpha/beta/rc tags
release = "0.6.0"
# The short X.Y version
version = "0.6.0"
# -- General configuration ---------------------------------------------------


# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx_rtd_theme",
    "sphinx.ext.napoleon",  # Google and NumPy-style docstrings
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.ifconfig",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.doctest",
    "myst_parser",
    "nbsphinx",  # Jupyter notebooks as docs
]

mathjax3_config = {"chtml": {"displayAlign": "left", "displayIndent": "2em"}}

autosectionlabel_prefix_document = True
autodoc_warningiserror = False  # suppress warnings during autodoc

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "apidoc/*tests*"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
html_css_files = ["custom.css"]

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
#
html_logo = "_static/NAWI_logo.png"

# The name of an image file (relative to this directory) to use as a favicon of
# the docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
#
html_favicon = "_static/favicon.ico"

## for MyST (Markdown)

myst_enable_extensions = [
    "dollarmath",
    "amsmath",
    "deflist",
    "html_admonition",
    "html_image",
    "colon_fence",
    "smartquotes",
    "replacements",
    "linkify",
    "substitution",
    "tasklist",
]
myst_heading_anchors = 2
myst_footnote_transition = True
myst_dmath_double_inline = True
panels_add_bootstrap_css = False


def run_apidoc(*args):
    # NOTE the env var must be set before importing apidoc, or the options
    # will have no effect
    os.environ["SPHINX_APIDOC_OPTIONS"] = "members,show-inheritance"
    from sphinx.ext import apidoc

    args = ["../watertap", "../watertap/*tests", "-o", "apidoc", "--force"]
    apidoc.main(args)


def skip(app, what, name, obj, would_skip, options):
    """Do not skip constructors!"""
    if name == "__init__":
        return False
    return would_skip


def setup(app):
    if os.environ.get("SKIP_APIDOC", False):
        print("Skipping apidoc")
    else:
        app.connect("builder-inited", run_apidoc)
    app.connect("autodoc-skip-member", skip)
