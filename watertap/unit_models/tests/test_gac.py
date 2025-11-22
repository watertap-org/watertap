#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
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
import pyomo.environ as pyo
import idaes.core.util.scaling as iscale

from idaes.core import (
    FlowsheetBlock,
    UnitModelCostingBlock,
)
from watertap.core.solvers import get_solver
from idaes.core.util.exceptions import ConfigurationError
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.unit_models.gac import GAC
from watertap.costing import WaterTAPCosting
from watertap.unit_models.tests.unit_test_harness import UnitTestHarness

__author__ = "Hunter Barber"

solver = get_solver()
zero = 1e-8
relative_tolerance = 1e-3


def build_hand():
    # trial problem from Hand, 1984 for removal of trace DCE
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = MCASParameterBlock(
        solute_list=["DCE"],
        mw_data={"H2O": 0.018, "DCE": 0.09896},
        ignore_neutral_charge=True,
    )
    m.fs.unit = GAC(
        property_package=m.fs.properties,
        film_transfer_coefficient_type="fixed",
        surface_diffusion_coefficient_type="fixed",
    )

    # feed specifications
    unit_feed = m.fs.unit.process_flow.properties_in[0]
    unit_feed.pressure.fix(101325)
    unit_feed.temperature.fix(273.15 + 25)
    unit_feed.flow_mol_phase_comp["Liq", "H2O"].fix(55555.55426666667)
    unit_feed.flow_mol_phase_comp["Liq", "DCE"].fix(0.0002344381568310428)

    # adsorption isotherm
    m.fs.unit.freund_k.fix(37.9e-6 * (1e6**0.8316))
    m.fs.unit.freund_ninv.fix(0.8316)
    # gac particle specifications
    m.fs.unit.particle_dens_app.fix(722)
    m.fs.unit.particle_dia.fix(0.00106)
    # adsorber bed specifications
    m.fs.unit.ebct.fix(300)  # seconds
    m.fs.unit.bed_voidage.fix(0.449)
    m.fs.unit.bed_length.fix(6)  # assumed
    # design spec
    m.fs.unit.conc_ratio_replace.fix(0.50)
    # parameters
    m.fs.unit.kf.fix(3.29e-5)
    m.fs.unit.ds.fix(1.77e-13)
    m.fs.unit.a0.fix(3.68421)
    m.fs.unit.a1.fix(13.1579)
    m.fs.unit.b0.fix(0.784576)
    m.fs.unit.b1.fix(0.239663)
    m.fs.unit.b2.fix(0.484422)
    m.fs.unit.b3.fix(0.003206)
    m.fs.unit.b4.fix(0.134987)

    # scaling
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e-4, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e4, index=("Liq", "DCE")
    )
    iscale.calculate_scaling_factors(m)

    return m


class TestGACHand(UnitTestHarness):
    def configure(self):
        m = build_hand()

        # arguments for UnitTestHarness
        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance

        # Approx data pulled from graph in Hand, 1984 at ~30 days
        # 30 days adjusted to actual solution to account for web plot data extraction error within reason
        # values calculated by hand and match those reported in Hand, 1984
        self.unit_solutions[m.fs.unit.equil_conc] = 0.0005178
        self.unit_solutions[m.fs.unit.dg] = 19780
        self.unit_solutions[m.fs.unit.N_Bi] = 6.113
        self.unit_solutions[m.fs.unit.min_N_St] = 35.68
        self.unit_solutions[m.fs.unit.throughput] = 0.9882
        self.unit_solutions[m.fs.unit.min_residence_time] = 468.4
        self.unit_solutions[m.fs.unit.residence_time] = 134.7
        self.unit_solutions[m.fs.unit.min_operational_time] = 9153000
        self.unit_solutions[m.fs.unit.operational_time] = 2554000
        self.unit_solutions[m.fs.unit.bed_volumes_treated] = 8514

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
                    "Liq", "DCE"
                ],
                "out": m.fs.unit.process_flow.properties_out[0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_out[0].flow_mol_phase_comp[
                    "Liq", "DCE"
                ],
            },
        }

        return m


def build_hand_surrogate():
    # trial problem from Hand, 1984 for removal of trace DCE
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = MCASParameterBlock(
        solute_list=["DCE"],
        mw_data={"H2O": 0.018, "DCE": 0.09896},
        ignore_neutral_charge=True,
    )
    m.fs.unit = GAC(
        property_package=m.fs.properties,
        film_transfer_coefficient_type="fixed",
        surface_diffusion_coefficient_type="fixed",
        cphsdm_calaculation_method="surrogate",
    )

    # feed specifications
    unit_feed = m.fs.unit.process_flow.properties_in[0]
    unit_feed.pressure.fix(101325)
    unit_feed.temperature.fix(273.15 + 25)
    unit_feed.flow_mol_phase_comp["Liq", "H2O"].fix(55555.55426666667)
    unit_feed.flow_mol_phase_comp["Liq", "DCE"].fix(0.0002344381568310428)

    # adsorption isotherm
    m.fs.unit.freund_k.fix(37.9e-6 * (1e6**0.8316))
    m.fs.unit.freund_ninv.fix(0.8316)
    # gac particle specifications
    m.fs.unit.particle_dens_app.fix(722)
    m.fs.unit.particle_dia.fix(0.00106)
    # adsorber bed specifications
    m.fs.unit.ebct.fix(300)  # seconds
    m.fs.unit.bed_voidage.fix(0.449)
    m.fs.unit.bed_length.fix(6)  # assumed
    # design spec
    m.fs.unit.conc_ratio_replace.fix(0.50)
    # parameters
    m.fs.unit.kf.fix(3.29e-5)
    m.fs.unit.ds.fix(1.77e-13)

    # scaling
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e-4, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e4, index=("Liq", "DCE")
    )
    iscale.calculate_scaling_factors(m)

    return m


class TestGACHandSurrogate(UnitTestHarness):
    def configure(self):
        m = build_hand_surrogate()

        # arguments for UnitTestHarness
        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance

        # Approx data pulled from graph in Hand, 1984 at ~30 days
        # 30 days adjusted to actual solution to account for web plot data extraction error within reason
        # values calculated by hand and match those reported in Hand, 1984
        self.unit_solutions[m.fs.unit.equil_conc] = 0.0005178
        self.unit_solutions[m.fs.unit.dg] = 19775.77
        self.unit_solutions[m.fs.unit.N_Bi] = 6.113
        self.unit_solutions[m.fs.unit.min_N_St] = 47.760
        self.unit_solutions[m.fs.unit.throughput] = 0.9892
        self.unit_solutions[m.fs.unit.min_residence_time] = 626.96
        self.unit_solutions[m.fs.unit.residence_time] = 134.7
        self.unit_solutions[m.fs.unit.min_operational_time] = 12266003
        self.unit_solutions[m.fs.unit.operational_time] = 2530674
        self.unit_solutions[m.fs.unit.bed_volumes_treated] = 8435

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
                    "Liq", "DCE"
                ],
                "out": m.fs.unit.process_flow.properties_out[0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_out[0].flow_mol_phase_comp[
                    "Liq", "DCE"
                ],
            },
        }

        return m


def build_hand_surrogate_no_ss():
    # trial problem from Hand, 1984 for removal of trace DCE
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = MCASParameterBlock(
        solute_list=["DCE"],
        mw_data={"H2O": 0.018, "DCE": 0.09896},
        ignore_neutral_charge=True,
    )
    m.fs.unit = GAC(
        property_package=m.fs.properties,
        film_transfer_coefficient_type="fixed",
        surface_diffusion_coefficient_type="fixed",
        cphsdm_calaculation_method="surrogate",
        add_trapezoidal_effluent_approximation=False,
    )

    # feed specifications
    unit_feed = m.fs.unit.process_flow.properties_in[0]
    unit_feed.pressure.fix(101325)
    unit_feed.temperature.fix(273.15 + 25)
    unit_feed.flow_mol_phase_comp["Liq", "H2O"].fix(55555.55426666667)
    unit_feed.flow_mol_phase_comp["Liq", "DCE"].fix(0.0002344381568310428)

    # adsorption isotherm
    m.fs.unit.freund_k.fix(37.9e-6 * (1e6**0.8316))
    m.fs.unit.freund_ninv.fix(0.8316)
    # gac particle specifications
    m.fs.unit.particle_dens_app.fix(722)
    m.fs.unit.particle_dia.fix(0.00106)
    # adsorber bed specifications
    m.fs.unit.ebct.fix(300)  # seconds
    m.fs.unit.bed_voidage.fix(0.449)
    m.fs.unit.bed_length.fix(6)  # assumed
    # design spec
    m.fs.unit.conc_ratio_replace.fix(0.50)
    # parameters
    m.fs.unit.kf.fix(3.29e-5)
    m.fs.unit.ds.fix(1.77e-13)

    # scaling
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e-4, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e4, index=("Liq", "DCE")
    )
    iscale.calculate_scaling_factors(m)

    return m


class TestGACHandSurrogateNoSS(UnitTestHarness):
    def configure(self):
        m = build_hand_surrogate_no_ss()

        # arguments for UnitTestHarness
        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance

        # Approx data pulled from graph in Hand, 1984 at ~30 days
        # 30 days adjusted to actual solution to account for web plot data extraction error within reason
        # values calculated by hand and match those reported in Hand, 1984
        self.unit_solutions[m.fs.unit.equil_conc] = 0.0005178
        self.unit_solutions[m.fs.unit.dg] = 19775.77
        self.unit_solutions[m.fs.unit.N_Bi] = 6.113
        self.unit_solutions[m.fs.unit.min_N_St] = 47.760
        self.unit_solutions[m.fs.unit.throughput] = 0.9892
        self.unit_solutions[m.fs.unit.min_residence_time] = 626.96
        self.unit_solutions[m.fs.unit.residence_time] = 134.7
        self.unit_solutions[m.fs.unit.min_operational_time] = 12266003
        self.unit_solutions[m.fs.unit.operational_time] = 2530674
        self.unit_solutions[m.fs.unit.bed_volumes_treated] = 8435

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
                    "Liq", "DCE"
                ],
                "out": m.fs.unit.process_flow.properties_out[0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_out[0].flow_mol_phase_comp[
                    "Liq", "DCE"
                ],
            },
        }

        return m


# -----------------------------------------------------------------------------
def build_crittenden():
    # trial problem from Crittenden, 2012 for removal of TCE
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = MCASParameterBlock(
        solute_list=["TCE"],
        mw_data={"H2O": 0.018, "TCE": 0.1314},
        diffus_calculation="HaydukLaudie",
        molar_volume_data={("Liq", "TCE"): 9.81e-5},
        ignore_neutral_charge=True,
    )
    m.fs.properties.visc_d_phase["Liq"] = 1.3097e-3
    m.fs.properties.dens_mass_const = 999.7
    m.fs.unit = GAC(
        property_package=m.fs.properties,
        film_transfer_coefficient_type="fixed",
        surface_diffusion_coefficient_type="fixed",
        finite_elements_ss_approximation=9,
    )

    # feed specifications
    unit_feed = m.fs.unit.process_flow.properties_in[0]
    unit_feed.pressure.fix(101325)
    unit_feed.temperature.fix(273.15 + 25)
    unit_feed.flow_mol_phase_comp["Liq", "H2O"].fix(823.8)
    unit_feed.flow_mol_phase_comp["Liq", "TCE"].fix(5.6444e-05)

    # adsorption isotherm
    m.fs.unit.freund_k.fix(1062e-6 * (1e6**0.48))
    m.fs.unit.freund_ninv.fix(0.48)
    # gac particle specifications
    m.fs.unit.particle_dens_app.fix(803.4)
    m.fs.unit.particle_dia.fix(0.001026)
    # adsorber bed specifications
    m.fs.unit.ebct.fix(10 * 60)
    m.fs.unit.bed_voidage.fix(0.44)
    m.fs.unit.velocity_sup.fix(5 / 3600)
    # design spec
    m.fs.unit.conc_ratio_replace.fix(0.80)
    # parameters
    m.fs.unit.ds.fix(1.24e-14)
    m.fs.unit.kf.fix(3.73e-05)
    m.fs.unit.a0.fix(0.8)
    m.fs.unit.a1.fix(0)
    m.fs.unit.b0.fix(0.023)
    m.fs.unit.b1.fix(0.793673)
    m.fs.unit.b2.fix(0.039324)
    m.fs.unit.b3.fix(0.009326)
    m.fs.unit.b4.fix(0.08275)

    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e-2, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e5, index=("Liq", "TCE")
    )
    iscale.calculate_scaling_factors(m)

    return m


class TestGACCrittenden(UnitTestHarness):
    def configure(self):
        m = build_crittenden()

        # arguments for UnitTestHarness
        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance

        # values calculated by hand and match those reported in Crittenden, 2012
        self.unit_solutions[m.fs.unit.equil_conc] = 0.02097
        self.unit_solutions[m.fs.unit.dg] = 42890
        self.unit_solutions[m.fs.unit.N_Bi] = 45.79
        self.unit_solutions[m.fs.unit.min_N_St] = 36.64
        self.unit_solutions[m.fs.unit.throughput] = 1.139
        self.unit_solutions[m.fs.unit.min_residence_time] = 395.9
        self.unit_solutions[m.fs.unit.residence_time] = 264.0
        self.unit_solutions[m.fs.unit.min_operational_time] = 19340000
        self.unit_solutions[m.fs.unit.operational_time] = 13690000
        self.unit_solutions[m.fs.unit.bed_volumes_treated] = 22810
        self.unit_solutions[m.fs.unit.velocity_int] = 0.003157
        self.unit_solutions[m.fs.unit.bed_length] = 0.8333
        self.unit_solutions[m.fs.unit.bed_area] = 10.68
        self.unit_solutions[m.fs.unit.bed_volume] = 8.900
        self.unit_solutions[m.fs.unit.bed_diameter] = 3.688
        self.unit_solutions[m.fs.unit.bed_mass_gac] = 4004
        self.unit_solutions[m.fs.unit.conc_ratio_avg] = 0.2287
        self.unit_solutions[m.fs.unit.ele_operational_time[1]] = 6462000

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
                    "Liq", "TCE"
                ],
                "out": m.fs.unit.process_flow.properties_out[0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_out[0].flow_mol_phase_comp[
                    "Liq", "TCE"
                ],
            },
        }

        return m


def build_crittenden_surrogate():
    # trial problem from Crittenden, 2012 for removal of TCE
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = MCASParameterBlock(
        solute_list=["TCE"],
        mw_data={"H2O": 0.018, "TCE": 0.1314},
        diffus_calculation="HaydukLaudie",
        molar_volume_data={("Liq", "TCE"): 9.81e-5},
        ignore_neutral_charge=True,
    )
    m.fs.properties.visc_d_phase["Liq"] = 1.3097e-3
    m.fs.properties.dens_mass_const = 999.7
    m.fs.unit = GAC(
        property_package=m.fs.properties,
        film_transfer_coefficient_type="fixed",
        surface_diffusion_coefficient_type="fixed",
        finite_elements_ss_approximation=9,
        cphsdm_calaculation_method="surrogate",
    )

    # feed specifications
    unit_feed = m.fs.unit.process_flow.properties_in[0]
    unit_feed.pressure.fix(101325)
    unit_feed.temperature.fix(273.15 + 25)
    unit_feed.flow_mol_phase_comp["Liq", "H2O"].fix(823.8)
    unit_feed.flow_mol_phase_comp["Liq", "TCE"].fix(5.6444e-05)

    # adsorption isotherm
    m.fs.unit.freund_k.fix(1062e-6 * (1e6**0.48))
    m.fs.unit.freund_ninv.fix(0.48)
    # gac particle specifications
    m.fs.unit.particle_dens_app.fix(803.4)
    m.fs.unit.particle_dia.fix(0.001026)
    # adsorber bed specifications
    m.fs.unit.ebct.fix(10 * 60)
    m.fs.unit.bed_voidage.fix(0.44)
    m.fs.unit.velocity_sup.fix(5 / 3600)
    # design spec
    m.fs.unit.conc_ratio_replace.fix(0.80)
    # parameters
    m.fs.unit.ds.fix(1.24e-14)
    m.fs.unit.kf.fix(3.73e-05)

    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e-2, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e5, index=("Liq", "TCE")
    )
    iscale.calculate_scaling_factors(m)

    return m


class TestGACCrittendenSurrogate(UnitTestHarness):
    def configure(self):
        m = build_crittenden_surrogate()

        # arguments for UnitTestHarness
        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance

        # values calculated by hand and match those reported in Crittenden, 2012
        self.unit_solutions[m.fs.unit.equil_conc] = 0.02097
        self.unit_solutions[m.fs.unit.dg] = 42886
        self.unit_solutions[m.fs.unit.N_Bi] = 45.79
        self.unit_solutions[m.fs.unit.min_N_St] = 32.09
        self.unit_solutions[m.fs.unit.throughput] = 1.084
        self.unit_solutions[m.fs.unit.min_residence_time] = 346.9
        self.unit_solutions[m.fs.unit.residence_time] = 264.0
        self.unit_solutions[m.fs.unit.min_operational_time] = 16127386
        self.unit_solutions[m.fs.unit.operational_time] = 12573768
        self.unit_solutions[m.fs.unit.bed_volumes_treated] = 20956
        self.unit_solutions[m.fs.unit.velocity_int] = 0.0031565
        self.unit_solutions[m.fs.unit.bed_length] = 0.8333
        self.unit_solutions[m.fs.unit.bed_area] = 10.68
        self.unit_solutions[m.fs.unit.bed_volume] = 8.900
        self.unit_solutions[m.fs.unit.bed_diameter] = 3.688
        self.unit_solutions[m.fs.unit.bed_mass_gac] = 4004
        self.unit_solutions[m.fs.unit.conc_ratio_avg] = 0.20453
        self.unit_solutions[m.fs.unit.ele_operational_time[1]] = 7830550

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
                    "Liq", "TCE"
                ],
                "out": m.fs.unit.process_flow.properties_out[0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_out[0].flow_mol_phase_comp[
                    "Liq", "TCE"
                ],
            },
        }

        return m


# -----------------------------------------------------------------------------
def build_multicomponent():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # inserting arbitrary BackGround Solutes, Cations, and Anions to check handling
    # arbitrary diffusivity data for non-target species
    m.fs.properties = MCASParameterBlock(
        solute_list=["TCE", "BGSOL", "BGCAT", "BGAN"],
        mw_data={
            "H2O": 0.018,
            "TCE": 0.1314,
            "BGSOL": 0.1,
            "BGCAT": 0.1,
            "BGAN": 0.1,
        },
        charge={"BGCAT": 1, "BGAN": -2},
        diffus_calculation="HaydukLaudie",
        molar_volume_data={("Liq", "TCE"): 9.81e-5},
        diffusivity_data={
            ("Liq", "BGSOL"): 1e-5,
            ("Liq", "BGCAT"): 1e-5,
            ("Liq", "BGAN"): 1e-5,
        },
        ignore_neutral_charge=True,
    )
    m.fs.properties.visc_d_phase["Liq"] = 1.3097e-3
    m.fs.properties.dens_mass_const = 1000
    # testing target_species arg
    m.fs.unit = GAC(
        property_package=m.fs.properties,
        film_transfer_coefficient_type="calculated",
        surface_diffusion_coefficient_type="calculated",
        target_species={"TCE"},
    )

    unit_feed = m.fs.unit.process_flow.properties_in[0]
    # feed specifications
    unit_feed.pressure.fix(101325)
    unit_feed.temperature.fix(273.15 + 25)
    unit_feed.flow_mol_phase_comp["Liq", "H2O"].fix(824.0736620370348)
    unit_feed.flow_mol_phase_comp["Liq", "TCE"].fix(5.644342973110135e-05)
    unit_feed.flow_mol_phase_comp["Liq", "BGSOL"].fix(5e-05)
    unit_feed.flow_mol_phase_comp["Liq", "BGCAT"].fix(2e-05)
    unit_feed.flow_mol_phase_comp["Liq", "BGAN"].fix(1e-05)

    # trial problem from Crittenden, 2012 for removal of TCE
    # adsorption isotherm
    m.fs.unit.freund_k.fix(1062e-6 * (1e6**0.48))
    m.fs.unit.freund_ninv.fix(0.48)
    # gac particle specifications
    m.fs.unit.particle_dens_app.fix(803.4)
    m.fs.unit.particle_dia.fix(0.001026)
    # adsorber bed specifications
    m.fs.unit.ebct.fix(10 * 60)
    m.fs.unit.bed_voidage.fix(0.44)
    m.fs.unit.velocity_sup.fix(5 / 3600)
    # design spec
    m.fs.unit.conc_ratio_replace.fix(0.80)
    # parameters
    m.fs.unit.particle_porosity.fix(0.641)
    m.fs.unit.tort.fix(1)
    m.fs.unit.spdfr.fix(1)
    m.fs.unit.shape_correction_factor.fix(1.5)
    m.fs.unit.a0.fix(0.8)
    m.fs.unit.a1.fix(0)
    m.fs.unit.b0.fix(0.023)
    m.fs.unit.b1.fix(0.793673)
    m.fs.unit.b2.fix(0.039324)
    m.fs.unit.b3.fix(0.009326)
    m.fs.unit.b4.fix(0.08275)

    # scaling
    prop = m.fs.properties
    prop.set_default_scaling("flow_mol_phase_comp", 1e-2, index=("Liq", "H2O"))
    for j in prop.ion_set | prop.solute_set:
        prop.set_default_scaling("flow_mol_phase_comp", 1e5, index=("Liq", j))

    iscale.calculate_scaling_factors(m)

    return m


class TestGACMultiComponent(UnitTestHarness):
    def configure(self):
        m = build_multicomponent()

        # arguments for UnitTestHarness
        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance

        self.unit_solutions[m.fs.unit.N_Re] = 2.473
        self.unit_solutions[m.fs.unit.N_Sc] = 2001
        self.unit_solutions[m.fs.unit.kf] = 2.600e-5
        self.unit_solutions[m.fs.unit.ds] = 1.245e-14

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
                    "Liq", "TCE"
                ]
                + m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
                    "Liq", "BGSOL"
                ]
                + m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
                    "Liq", "BGCAT"
                ]
                + m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
                    "Liq", "BGAN"
                ],
                "out": m.fs.unit.process_flow.properties_out[0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_out[0].flow_mol_phase_comp[
                    "Liq", "TCE"
                ]
                + m.fs.unit.process_flow.properties_out[0].flow_mol_phase_comp[
                    "Liq", "BGSOL"
                ]
                + m.fs.unit.process_flow.properties_out[0].flow_mol_phase_comp[
                    "Liq", "BGCAT"
                ]
                + m.fs.unit.process_flow.properties_out[0].flow_mol_phase_comp[
                    "Liq", "BGAN"
                ],
            },
        }

        return m


# -----------------------------------------------------------------------------
def build_multicomponent_surrogate():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # inserting arbitrary BackGround Solutes, Cations, and Anions to check handling
    # arbitrary diffusivity data for non-target species
    m.fs.properties = MCASParameterBlock(
        solute_list=["TCE", "BGSOL", "BGCAT", "BGAN"],
        mw_data={
            "H2O": 0.018,
            "TCE": 0.1314,
            "BGSOL": 0.1,
            "BGCAT": 0.1,
            "BGAN": 0.1,
        },
        charge={"BGCAT": 1, "BGAN": -2},
        diffus_calculation="HaydukLaudie",
        molar_volume_data={("Liq", "TCE"): 9.81e-5},
        diffusivity_data={
            ("Liq", "BGSOL"): 1e-5,
            ("Liq", "BGCAT"): 1e-5,
            ("Liq", "BGAN"): 1e-5,
        },
        ignore_neutral_charge=True,
    )
    m.fs.properties.visc_d_phase["Liq"] = 1.3097e-3
    m.fs.properties.dens_mass_const = 1000
    # testing target_species arg
    m.fs.unit = GAC(
        property_package=m.fs.properties,
        film_transfer_coefficient_type="calculated",
        surface_diffusion_coefficient_type="calculated",
        target_species={"TCE"},
        cphsdm_calaculation_method="surrogate",
    )

    unit_feed = m.fs.unit.process_flow.properties_in[0]
    # feed specifications
    unit_feed.pressure.fix(101325)
    unit_feed.temperature.fix(273.15 + 25)
    unit_feed.flow_mol_phase_comp["Liq", "H2O"].fix(824.0736620370348)
    unit_feed.flow_mol_phase_comp["Liq", "TCE"].fix(5.644342973110135e-05)
    unit_feed.flow_mol_phase_comp["Liq", "BGSOL"].fix(5e-05)
    unit_feed.flow_mol_phase_comp["Liq", "BGCAT"].fix(2e-05)
    unit_feed.flow_mol_phase_comp["Liq", "BGAN"].fix(1e-05)

    # trial problem from Crittenden, 2012 for removal of TCE
    # adsorption isotherm
    m.fs.unit.freund_k.fix(1062e-6 * (1e6**0.48))
    m.fs.unit.freund_ninv.fix(0.48)
    # gac particle specifications
    m.fs.unit.particle_dens_app.fix(803.4)
    m.fs.unit.particle_dia.fix(0.001026)
    # adsorber bed specifications
    m.fs.unit.ebct.fix(10 * 60)
    m.fs.unit.bed_voidage.fix(0.44)
    m.fs.unit.velocity_sup.fix(5 / 3600)
    # design spec
    m.fs.unit.conc_ratio_replace.fix(0.80)
    # parameters
    m.fs.unit.particle_porosity.fix(0.641)
    m.fs.unit.tort.fix(1)
    m.fs.unit.spdfr.fix(1)
    m.fs.unit.shape_correction_factor.fix(1.5)

    # scaling
    prop = m.fs.properties
    prop.set_default_scaling("flow_mol_phase_comp", 1e-2, index=("Liq", "H2O"))
    for j in prop.ion_set | prop.solute_set:
        prop.set_default_scaling("flow_mol_phase_comp", 1e5, index=("Liq", j))

    iscale.calculate_scaling_factors(m)

    return m


class TestGACMultiComponentSurrogate(UnitTestHarness):
    def configure(self):
        m = build_multicomponent_surrogate()

        # arguments for UnitTestHarness
        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance

        self.unit_solutions[m.fs.unit.freund_k] = 0.805
        self.unit_solutions[m.fs.unit.freund_ninv] = 0.48
        self.unit_solutions[m.fs.unit.ds] = 1.244e-14
        self.unit_solutions[m.fs.unit.kf] = 2.598e-05
        self.unit_solutions[m.fs.unit.equil_conc] = 0.02097
        self.unit_solutions[m.fs.unit.dg] = 42887
        self.unit_solutions[m.fs.unit.N_Bi] = 31.78
        self.unit_solutions[m.fs.unit.ebct] = 600
        self.unit_solutions[m.fs.unit.residence_time] = 264
        self.unit_solutions[m.fs.unit.bed_mass_gac] = 4004
        self.unit_solutions[m.fs.unit.particle_dens_app] = 803
        self.unit_solutions[m.fs.unit.particle_dens_bulk] = 449.9
        self.unit_solutions[m.fs.unit.particle_dia] = 0.001026
        self.unit_solutions[m.fs.unit.min_N_St] = 21.97
        self.unit_solutions[m.fs.unit.min_ebct] = 774
        self.unit_solutions[m.fs.unit.throughput] = 1.11796
        self.unit_solutions[m.fs.unit.min_residence_time] = 340.84
        self.unit_solutions[m.fs.unit.min_operational_time] = 16342446
        self.unit_solutions[m.fs.unit.operational_time] = 13046915
        self.unit_solutions[m.fs.unit.bed_volumes_treated] = 21744
        self.unit_solutions[m.fs.unit.ele_throughput[1]] = 0.7368
        self.unit_solutions[m.fs.unit.ele_throughput[2]] = 0.8259
        self.unit_solutions[m.fs.unit.ele_throughput[3]] = 0.8842
        self.unit_solutions[m.fs.unit.ele_throughput[4]] = 0.9609
        self.unit_solutions[m.fs.unit.ele_throughput[5]] = 1.118
        self.unit_solutions[m.fs.unit.ele_min_operational_time[1]] = 10771565
        self.unit_solutions[m.fs.unit.ele_min_operational_time[2]] = 12074403
        self.unit_solutions[m.fs.unit.ele_min_operational_time[3]] = 12926494
        self.unit_solutions[m.fs.unit.ele_min_operational_time[4]] = 14046498
        self.unit_solutions[m.fs.unit.ele_min_operational_time[5]] = 16342446
        self.unit_solutions[m.fs.unit.ele_operational_time[0]] = 0
        self.unit_solutions[m.fs.unit.ele_operational_time[1]] = 7476035
        self.unit_solutions[m.fs.unit.ele_operational_time[2]] = 8778872
        self.unit_solutions[m.fs.unit.ele_operational_time[3]] = 9630963
        self.unit_solutions[m.fs.unit.ele_operational_time[4]] = 10750968
        self.unit_solutions[m.fs.unit.ele_operational_time[5]] = 13046915
        self.unit_solutions[m.fs.unit.ele_conc_ratio_term[1]] = 0.0028650
        self.unit_solutions[m.fs.unit.ele_conc_ratio_term[2]] = 0.0108595
        self.unit_solutions[m.fs.unit.ele_conc_ratio_term[3]] = 0.0200011
        self.unit_solutions[m.fs.unit.ele_conc_ratio_term[4]] = 0.0432441
        self.unit_solutions[m.fs.unit.ele_conc_ratio_term[5]] = 0.1234034
        self.unit_solutions[m.fs.unit.conc_ratio_avg] = 0.2003
        self.unit_solutions[m.fs.unit.mass_adsorbed] = 77.37
        self.unit_solutions[m.fs.unit.gac_usage_rate] = 0.000307
        self.unit_solutions[m.fs.unit.N_Re] = 2.47280
        self.unit_solutions[m.fs.unit.N_Sc] = 2001

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
                    "Liq", "TCE"
                ]
                + m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
                    "Liq", "BGSOL"
                ]
                + m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
                    "Liq", "BGCAT"
                ]
                + m.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp[
                    "Liq", "BGAN"
                ],
                "out": m.fs.unit.process_flow.properties_out[0].flow_mol_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.process_flow.properties_out[0].flow_mol_phase_comp[
                    "Liq", "TCE"
                ]
                + m.fs.unit.process_flow.properties_out[0].flow_mol_phase_comp[
                    "Liq", "BGSOL"
                ]
                + m.fs.unit.process_flow.properties_out[0].flow_mol_phase_comp[
                    "Liq", "BGCAT"
                ]
                + m.fs.unit.process_flow.properties_out[0].flow_mol_phase_comp[
                    "Liq", "BGAN"
                ],
            },
        }

        return m


# -----------------------------------------------------------------------------
class TestGACErrorLog:
    @pytest.mark.unit
    def test_error(self):

        with pytest.raises(
            ConfigurationError,
            match="'target species' is not specified for the GAC unit model, "
            "either specify 'target species' argument or reduce solute set "
            "to a single component",
        ):
            m = pyo.ConcreteModel()
            m.fs = FlowsheetBlock(dynamic=False)

            # inserting arbitrary BackGround Solutes, Cations, and Anions to check handling
            # arbitrary diffusivity data for non-target species
            m.fs.properties = MCASParameterBlock(
                solute_list=["TCE", "BGSOL", "BGCAT", "BGAN"],
                mw_data={
                    "H2O": 0.018,
                    "TCE": 0.1314,
                    "BGSOL": 0.1,
                    "BGCAT": 0.1,
                    "BGAN": 0.1,
                },
                charge={"BGCAT": 1, "BGAN": -2},
                diffus_calculation="HaydukLaudie",
                molar_volume_data={("Liq", "TCE"): 9.81e-5},
                diffusivity_data={
                    ("Liq", "BGSOL"): 1e-5,
                    ("Liq", "BGCAT"): 1e-5,
                    ("Liq", "BGAN"): 1e-5,
                },
                ignore_neutral_charge=True,
            )
            m.fs.properties.visc_d_phase["Liq"] = 1.3097e-3
            m.fs.properties.dens_mass_const = 1000

            # testing target_species arg
            m.fs.unit = GAC(
                property_package=m.fs.properties,
                film_transfer_coefficient_type="calculated",
                surface_diffusion_coefficient_type="calculated",
            )

        with pytest.raises(
            ConfigurationError,
            match="fs.unit received invalid argument for contactor_type:"
            " vessel. Argument must be a member of the ContactorType Enum.",
        ):
            m = pyo.ConcreteModel()
            m.fs = FlowsheetBlock(dynamic=False)

            m.fs.properties = MCASParameterBlock(
                solute_list=["TCE"],
                mw_data={"H2O": 0.018, "TCE": 0.1314},
                diffus_calculation="HaydukLaudie",
                molar_volume_data={("Liq", "TCE"): 9.81e-5},
                charge={"TCE": 0},
            )
            m.fs.properties.visc_d_phase["Liq"] = 1.3097e-3
            m.fs.properties.dens_mass_const = 1000

            m.fs.unit = GAC(
                property_package=m.fs.properties,
                film_transfer_coefficient_type="calculated",
                surface_diffusion_coefficient_type="calculated",
            )

            m.fs.costing = WaterTAPCosting()
            m.fs.costing.base_currency = pyo.units.USD_2020

            m.fs.unit.costing = UnitModelCostingBlock(
                flowsheet_costing_block=m.fs.costing,
                costing_method_arguments={"contactor_type": "vessel"},
            )

        with pytest.raises(
            ConfigurationError,
            match="item 0 within 'target_species' list is not of data type str",
        ):
            m = pyo.ConcreteModel()
            m.fs = FlowsheetBlock(dynamic=False)

            # inserting arbitrary BackGround Solutes, Cations, and Anions to check handling
            # arbitrary diffusivity data for non-target species
            m.fs.properties = MCASParameterBlock(
                solute_list=["TCE", "BGSOL", "BGCAT", "BGAN"],
                mw_data={
                    "H2O": 0.018,
                    "TCE": 0.1314,
                    "BGSOL": 0.1,
                    "BGCAT": 0.1,
                    "BGAN": 0.1,
                },
                charge={"BGCAT": 1, "BGAN": -2},
                diffus_calculation="HaydukLaudie",
                molar_volume_data={("Liq", "TCE"): 9.81e-5},
                diffusivity_data={
                    ("Liq", "BGSOL"): 1e-5,
                    ("Liq", "BGCAT"): 1e-5,
                    ("Liq", "BGAN"): 1e-5,
                },
                ignore_neutral_charge=True,
            )
            m.fs.properties.visc_d_phase["Liq"] = 1.3097e-3
            m.fs.properties.dens_mass_const = 1000

            # testing target_species arg
            m.fs.unit = GAC(
                property_package=m.fs.properties,
                film_transfer_coefficient_type="calculated",
                surface_diffusion_coefficient_type="calculated",
                target_species=range(2),
            )

        with pytest.raises(
            ConfigurationError,
            match="item species within 'target_species' list is not in 'component_list",
        ):
            m = pyo.ConcreteModel()
            m.fs = FlowsheetBlock(dynamic=False)

            # inserting arbitrary BackGround Solutes, Cations, and Anions to check handling
            # arbitrary diffusivity data for non-target species
            m.fs.properties = MCASParameterBlock(
                solute_list=["TCE", "BGSOL", "BGCAT", "BGAN"],
                mw_data={
                    "H2O": 0.018,
                    "TCE": 0.1314,
                    "BGSOL": 0.1,
                    "BGCAT": 0.1,
                    "BGAN": 0.1,
                },
                charge={"BGCAT": 1, "BGAN": -2},
                diffus_calculation="HaydukLaudie",
                molar_volume_data={("Liq", "TCE"): 9.81e-5},
                diffusivity_data={
                    ("Liq", "BGSOL"): 1e-5,
                    ("Liq", "BGCAT"): 1e-5,
                    ("Liq", "BGAN"): 1e-5,
                },
                ignore_neutral_charge=True,
            )
            m.fs.properties.visc_d_phase["Liq"] = 1.3097e-3
            m.fs.properties.dens_mass_const = 1000

            # testing target_species arg
            m.fs.unit = GAC(
                property_package=m.fs.properties,
                film_transfer_coefficient_type="calculated",
                surface_diffusion_coefficient_type="calculated",
                target_species={"species": "TCE"},
            )
