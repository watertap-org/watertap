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
    Simple unit tests for example flowsheet of SepRO with Chlorination.

    NOTE: That flowsheet is not meant to be viewed as a final product, but
    a sample of how to incorporate more complex chemistry into a simple
    flowsheet.
"""
import pytest

from watertap.examples.flowsheets.full_treatment_train.flowsheet_components.chemistry.posttreatment_ideal_naocl_chlorination_block import (
    run_ideal_naocl_mixer_example,
    run_ideal_naocl_chlorination_example,
    run_chlorination_block_example,
    build_translator_from_RO_to_chlorination_block,
    unfix_ideal_naocl_mixer_inlet_stream,
)
from watertap.examples.flowsheets.full_treatment_train.model_components import (
    property_models,
)
from watertap.examples.flowsheets.full_treatment_train.flowsheet_components import (
    desalination,
)
from watertap.examples.flowsheets.full_treatment_train.util import check_dof
from pyomo.environ import TransformationFactory
from pyomo.network import Arc

__author__ = "Austin Ladshaw"


@pytest.mark.component
def test_ideal_naocl_mixer():
    model = run_ideal_naocl_mixer_example()
    assert model.fs.ideal_naocl_mixer_unit.dosing_rate.value == pytest.approx(
        0.4e-6, rel=1e-3
    )
    assert model.fs.ideal_naocl_mixer_unit.outlet.flow_mol[0].value == pytest.approx(
        10.00001074690966, rel=1e-3
    )
    assert model.fs.ideal_naocl_mixer_unit.outlet.mole_frac_comp[
        0, "OCl_-"
    ].value == pytest.approx(5.373448801531194e-07, rel=1e-3)


@pytest.mark.component
def test_ideal_naocl_chlorination():
    model = run_ideal_naocl_chlorination_example()
    assert model.fs.ideal_naocl_chlorination_unit.free_chlorine.value == pytest.approx(
        1.893849296278168, rel=1e-3
    )
    assert model.fs.ideal_naocl_chlorination_unit.outlet.mole_frac_comp[
        0, "OCl_-"
    ].value == pytest.approx(4.254858076511e-07, rel=1e-3)
    assert model.fs.ideal_naocl_chlorination_unit.outlet.mole_frac_comp[
        0, "H_+"
    ].value == pytest.approx(5.6407676871845223e-11, rel=1e-3)


@pytest.mark.component
@pytest.mark.requires_idaes_solver
def test_ideal_naocl_chlorination_full_block():
    model = run_chlorination_block_example(fix_free_chlorine=True)
    assert model.fs.ideal_naocl_mixer_unit.dosing_rate.value == pytest.approx(
        0.9504457542440085e-6, rel=1e-3
    )
    assert model.fs.ideal_naocl_mixer_unit.outlet.flow_mol[0].value == pytest.approx(
        25.000025535888078, rel=1e-3
    )
    assert model.fs.ideal_naocl_mixer_unit.outlet.mole_frac_comp[
        0, "OCl_-"
    ].value == pytest.approx(5.107172398859054e-07, rel=1e-3)

    assert model.fs.ideal_naocl_chlorination_unit.free_chlorine.value == pytest.approx(
        2, rel=1e-3
    )
    assert model.fs.ideal_naocl_chlorination_unit.outlet.mole_frac_comp[
        0, "OCl_-"
    ].value == pytest.approx(4.508812116261189e-07, rel=1e-3)
    assert model.fs.ideal_naocl_chlorination_unit.outlet.mole_frac_comp[
        0, "H_+"
    ].value == pytest.approx(5.47976004923939e-11, rel=1e-3)


@pytest.mark.component
def test_addition_of_translator():
    model = run_ideal_naocl_mixer_example()
    property_models.build_prop(model, base="TDS")
    build_translator_from_RO_to_chlorination_block(model)
    assert hasattr(model.fs, "RO_to_Chlor")


@pytest.mark.component
def test_build_flowsheet():
    model = run_chlorination_block_example(fix_free_chlorine=True)
    property_models.build_prop(model, base="TDS")
    build_translator_from_RO_to_chlorination_block(model)
    # Here, we set 'has_feed' to True because RO is our first block in the flowsheet
    kwargs_desal = {
        "has_desal_feed": True,
        "is_twostage": False,
        "has_ERD": False,
        "RO_type": "0D",
        "RO_base": "TDS",
        "RO_level": "detailed",
    }
    desal_port = desalination.build_desalination(model, **kwargs_desal)
    desalination.scale_desalination(model, **kwargs_desal)
    desalination.initialize_desalination(model, **kwargs_desal)

    # Add the connecting arcs
    model.fs.S1 = Arc(source=desal_port["out"], destination=model.fs.RO_to_Chlor.inlet)
    model.fs.S2 = Arc(
        source=model.fs.RO_to_Chlor.outlet,
        destination=model.fs.ideal_naocl_mixer_unit.inlet_stream,
    )
    TransformationFactory("network.expand_arcs").apply_to(model)

    # Unfix inlet stream for mixer
    unfix_ideal_naocl_mixer_inlet_stream(model)
    check_dof(model)
