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
    Simple example of a flowsheet containing an RO separator unit model and
    a simple NaOCl chlorination post-treatment unit model.

    inlet ---> [RO Separator] ---> permeate ---> (Translator) ---> [Chlorination] ---> outlet
                    |
                    |
                    v
                retentate (i.e., waste)

    ---------- WORK IN PROGRESS ---------
"""

from proteuslib.flowsheets.full_treatment_train.example_models.unit_separator import (
    build_RO_separator_example, run_RO_example)
from proteuslib.flowsheets.full_treatment_train.chemistry_flowsheets.PostTreatment_SimpleNaOCl_Chlorination import (
    build_simple_naocl_chlorination_unit,
    initialize_chlorination_example,
    display_results_of_chlorination, run_chlorination_example)


if __name__ == "__main__":
    model = run_RO_example()
