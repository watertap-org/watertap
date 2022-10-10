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
    pretreatment_NF,
)
from watertap.examples.flowsheets.full_treatment_train.model_components import (
    property_models,
)
from watertap.examples.flowsheets.full_treatment_train.util import (
    check_build,
    check_scaling,
)


@pytest.mark.unit
def test_build_and_scale_pretreatment_NF():
    for has_bypass in [True, False]:
        for NF_type in ["Sep", "ZO"]:
            for NF_base in ["ion", "salt"]:
                kwargs = {
                    "has_bypass": has_bypass,
                    "NF_type": NF_type,
                    "NF_base": NF_base,
                }
                print("\n\n***kwargs****\n", kwargs)

                if NF_type == "ZO" and NF_base == "salt":
                    continue  # not a valid combination

                m = ConcreteModel()
                m.fs = FlowsheetBlock(dynamic=False)
                property_models.build_prop(m, base=kwargs["NF_base"])

                check_build(
                    m, build_func=pretreatment_NF.build_pretreatment_NF, **kwargs
                )
                assert hasattr(m.fs, "NF")
                check_scaling(
                    m, scale_func=pretreatment_NF.scale_pretreatment_NF, **kwargs
                )

                pretreatment_NF.display_pretreatment_NF(m, **kwargs)
