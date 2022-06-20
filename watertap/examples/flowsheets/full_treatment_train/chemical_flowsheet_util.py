###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

"""
    These are a set of utils for aiding in solving flowsheets with chemical
    reactions. Based on initiale testing, it was found that these tools may
    be necessary to get a flowsheet that includes electrolyte chemistry to
    converge properly.
"""

from pyomo.network import SequentialDecomposition
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver
from pyomo.environ import value

__author__ = "Austin Ladshaw"

# Get default solver for testing
solver = get_solver()


def set_H2O_molefraction(port):
    # Perform a summation of all non-H2O molefractions to find the H2O molefraction
    sum = 0
    for i in port.mole_frac_comp:
        # NOTE: i will be a tuple with format (time, component)
        if i[1] != "H2O":
            sum += value(port.mole_frac_comp[i[0], i[1]])

    port.mole_frac_comp[0, "H2O"].set_value(1 - sum)


def zero_out_non_H2O_molefractions(port):
    for i in port.mole_frac_comp:
        if i[1] != "H2O":
            port.mole_frac_comp[i[0], i[1]].set_value(0)


def fix_all_molefractions(port):
    for i in port.mole_frac_comp:
        port.mole_frac_comp[i[0], i[1]].fix()


def unfix_all_molefractions(port):
    for i in port.mole_frac_comp:
        port.mole_frac_comp[i[0], i[1]].unfix()


#   This tool for block decomp may be less useful moving forward.
#       probably better to initialize each manually.
def block_initializer(blk, tee=False):
    if tee:
        print(blk)
    results = solver.solve(blk, tee=tee)


def seq_decomp_initializer(model):
    seq = SequentialDecomposition(tol=1.0e-3)
    seq.options.select_tear_method = "heuristic"
    seq.run(model, block_initializer)
