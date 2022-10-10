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
    translator_block,
)
from watertap.examples.flowsheets.full_treatment_train.model_components import (
    property_models,
)
from watertap.examples.flowsheets.full_treatment_train.util import check_scaling


@pytest.mark.unit
def test_build_and_scale_translator_block():
    tb_kwargs_list = [
        {"base_inlet": "ion", "base_outlet": "TDS", "name_str": "tb_pretrt_to_desal"},
        {"base_inlet": "salt", "base_outlet": "TDS", "name_str": "tb_pretrt_to_desal"},
    ]
    for kwargs in tb_kwargs_list:
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        property_models.build_prop(m, base=kwargs["base_inlet"])
        property_models.build_prop(m, base=kwargs["base_outlet"])

        translator_block.build_tb(m, **kwargs)
        assert hasattr(m.fs, kwargs["name_str"])

        check_scaling(m, **kwargs)
