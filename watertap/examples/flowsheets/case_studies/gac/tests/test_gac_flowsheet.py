#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pytest

from idaes.core.solvers import get_solver
from watertap.examples.flowsheets.case_studies.gac import (
    gac_flowsheet as gac_fs,
)

__author__ = "Hunter Barber"

solver = get_solver()


class TestGACFlowsheet:
    @pytest.fixture(scope="class")
    def gac_frame(self):

        m = gac_fs.build(
            film_transfer_coefficient_type="calculated",
            surface_diffusion_coefficient_type="calculated",
            diffusivity_calculation="HaydukLaudie",
        )

        return m

    # placeholder
