###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
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
from idaes.core.util import get_solver

# Get default solver for testing
solver = get_solver()

def block_initializer(blk, tee=False):
    solver.options['bound_push'] = 1e-10
    solver.options['mu_init'] = 1e-6
    solver.options["nlp_scaling_method"] = "user-scaling"
    results = solver.solve(blk, tee=tee)
    iscale.constraint_autoscale_large_jac(blk)

def seq_decomp_initializer(model, tee=False):
    seq = SequentialDecomposition(tol=1.0E-3)
    seq.options.select_tear_method = "heuristic"
    seq.run(model, block_initializer)
    iscale.constraint_autoscale_large_jac(model)
