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
copyright = "2020-2024, NAWI"
author = "NAWI"

# The full version, including alpha/beta/rc tags
release = "1.3.dev0"
# The short X.Y version
version = "1.3.dev0"
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

# intersphinx mapping to idaes
intersphinx_mapping = {
    "idaes": ("https://idaes-pse.readthedocs.io/en/stable/", None),
    "parameter_sweep": ("https://parameter-sweep.readthedocs.io/en/latest/", None),
}

rst_epilog = """
.. |Binder launch button| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/watertap-org/watertap/main?labpath=tutorials%2F00-index.ipynb
"""


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
