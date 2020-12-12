# ProteusLib
The ProteusLib development repository.

ProteusLib is developed as part of the [National Alliance for Water Innovation](https://nawihub.org/) project.

## Developer notes

### Installation

To install from source (i.e., a `git clone` of the repository), 
first create and activate your python environment, e.g. `conda create -n proteus python=3`
and `conda activate proteus`. Next, change to the top-level directory with this README file and run pip
in "editable" mode, with the "dev" option:

    pip install -e .[dev]

### Tests

Tests are created for the "pytest" framework. This is a dependency and should be installed when you install
the rest of the code. To run, simply type at the top-level of the source code tree:

    pytest


