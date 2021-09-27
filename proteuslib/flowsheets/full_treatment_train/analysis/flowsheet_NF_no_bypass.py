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
import proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_NF as flowsheet_NF


# def build_components(m):
#     flowsheet_NF.build_components(m, has_bypass=False)
#
#
# def build(m):
#     flowsheet_NF.build(m, has_bypass=False)
#
#
# def scale(m):
#     flowsheet_NF.scale(m, has_bypass=False)
#
#
# def initialize(m):
#     flowsheet_NF.initialize(m, has_bypass=False)
#
#
# def report(m):
#     flowsheet_NF.report(m, has_bypass=False)


def solve_flowsheet():
    m = flowsheet_NF.solve_flowsheet(has_bypass=False)
    return m


def simulate(m, **kwargs):
    if kwargs is not None:
        if kwargs['unfix_nf_area']:
            m.fs.NF.area.unfix()

    flowsheet_NF.simulate(m)
    return m


if __name__ == "__main__":
    solve_flowsheet()
