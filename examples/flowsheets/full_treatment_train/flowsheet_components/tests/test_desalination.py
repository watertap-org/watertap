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
import pytest
from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock
from watertap.examples.flowsheets.full_treatment_train.flowsheet_components import (
    desalination,
)
from watertap.examples.flowsheets.full_treatment_train.model_components import (
    property_models,
)
from watertap.examples.flowsheets.full_treatment_train.util import (
    check_build,
    check_scaling,
)


@pytest.mark.unit
def test_build_and_scale_desalination():
    for has_desal_feed in [True]:
        for is_twostage in [True, False]:
            for has_ERD in [True, False]:
                for RO_type in ["Sep", "0D", "1D"]:
                    for RO_base in ["TDS"]:
                        for RO_level in ["simple", "detailed"]:
                            kwargs = {
                                "has_desal_feed": has_desal_feed,
                                "is_twostage": is_twostage,
                                "has_ERD": has_ERD,
                                "RO_type": RO_type,
                                "RO_base": RO_base,
                                "RO_level": RO_level,
                            }

                            if RO_type == "Sep" and is_twostage:
                                continue  # not a supported combination
                            elif RO_type == "Sep" and has_ERD:
                                continue  # not a supported combination
                            elif RO_type == "1D" and RO_level == "simple":
                                continue  # not a supported combination

                            m = ConcreteModel()
                            m.fs = FlowsheetBlock(dynamic=False)
                            property_models.build_prop(m, base=kwargs["RO_base"])

                            check_build(
                                m, build_func=desalination.build_desalination, **kwargs
                            )
                            assert hasattr(m.fs, "RO")
                            check_scaling(
                                m, scale_func=desalination.scale_desalination, **kwargs
                            )

                            desalination.display_desalination(m, **kwargs)
