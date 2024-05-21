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

from pyomo.environ import ConcreteModel

from idaes.core import FlowsheetBlock
import idaes.core.util.scaling as iscale

from watertap.core.solvers import get_solver
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.unit_models.ion_exchange_0D import IonExchange0D
from watertap.unit_models.tests.unit_test_harness import UnitTestHarness

__author__ = "Kurban Sitterley"

solver = get_solver()
zero = 1e-8
relative_tolerance = 1e-3


def build_langmuir():

    target_ion = "Ca_2+"
    ion_props = {
        "solute_list": [target_ion],
        "diffusivity_data": {("Liq", target_ion): 9.2e-10},
        "mw_data": {"H2O": 0.018, target_ion: 0.04},
        "charge": {target_ion: 2},
    }
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(**ion_props)
    ix_config = {
        "property_package": m.fs.properties,
        "target_ion": target_ion,
    }
    m.fs.unit = ix = IonExchange0D(**ix_config)
    ix.process_flow.properties_in.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 0.5,
            ("conc_mass_phase_comp", ("Liq", target_ion)): 0.1,
            ("pressure", None): 101325,
            ("temperature", None): 298,
        },
        hold_state=True,
    )

    ix.process_flow.properties_in[0].flow_mass_phase_comp[...]
    ix.process_flow.properties_out[0].flow_mass_phase_comp[...]
    ix.regeneration_stream[0].flow_mass_phase_comp[...]

    ix.service_flow_rate.fix(15)
    ix.langmuir[target_ion].fix(0.9)
    ix.resin_max_capacity.fix(3)
    ix.bed_depth.fix(1.7)
    ix.dimensionless_time.fix()
    ix.number_columns.fix(8)
    ix.resin_diam.fix()
    ix.resin_bulk_dens.fix()
    ix.bed_porosity.fix()

    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e-4, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 10, index=("Liq", "Ca_2+")
    )

    iscale.calculate_scaling_factors(m)

    return m


class TestIXLangmuir(UnitTestHarness):
    def configure(self):
        m = build_langmuir()

        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance

        self.unit_solutions[m.fs.unit.c_norm["Ca_2+"]] = 0.4919
        self.unit_solutions[m.fs.unit.col_height] = 3.488
        self.unit_solutions[m.fs.unit.col_diam] = 3.351
        self.unit_solutions[m.fs.unit.col_height_to_diam_ratio] = 1.0408
        self.unit_solutions[m.fs.unit.number_columns] = 8
        self.unit_solutions[m.fs.unit.t_breakthru] = 52360
        self.unit_solutions[m.fs.unit.vel_bed] = 0.007083
        self.unit_solutions[m.fs.unit.N_Re] = 4.958
        self.unit_solutions[m.fs.unit.N_Sc["Ca_2+"]] = 1086
        self.unit_solutions[m.fs.unit.N_Sh["Ca_2+"]] = 26.296
        self.unit_solutions[m.fs.unit.N_Pe_particle] = 0.10782
        self.unit_solutions[m.fs.unit.N_Pe_bed] = 261.8
        self.unit_solutions[m.fs.unit.resin_eq_capacity] = 1.554
        self.unit_solutions[m.fs.unit.resin_unused_capacity] = 1.445
        self.unit_solutions[m.fs.unit.langmuir["Ca_2+"]] = 0.9
        self.unit_solutions[m.fs.unit.mass_removed["Ca_2+"]] = 65300
        self.unit_solutions[m.fs.unit.num_transfer_units] = 35.5483
        self.unit_solutions[m.fs.unit.dimensionless_time] = 1
        self.unit_solutions[m.fs.unit.partition_ratio] = 217.669
        self.unit_solutions[m.fs.unit.fluid_mass_transfer_coeff["Ca_2+"]] = 3.4561e-05
        self.unit_solutions[m.fs.unit.bw_flow] = 0.09803
        self.unit_solutions[m.fs.unit.bed_expansion_frac] = 0.46395
        self.unit_solutions[m.fs.unit.col_vol_tot] = 246.26
        self.unit_solutions[m.fs.unit.lh] = 0.0
        self.unit_solutions[m.fs.unit.separation_factor["Ca_2+"]] = 1.111
        self.unit_solutions[m.fs.unit.rate_coeff["Ca_2+"]] = 0.00021159
        self.unit_solutions[m.fs.unit.HTU["Ca_2+"]] = 0.04782
        self.unit_solutions[
            m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 499.95
        self.unit_solutions[
            m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                "Liq", "Ca_2+"
            ]
        ] = 0.05
        self.unit_solutions[
            m.fs.unit.process_flow.properties_out[0.0].flow_mass_phase_comp[
                "Liq", "H2O"
            ]
        ] = 499.95
        self.unit_solutions[
            m.fs.unit.process_flow.properties_out[0.0].flow_mass_phase_comp[
                "Liq", "Ca_2+"
            ]
        ] = 1.1458e-4
        self.unit_solutions[
            m.fs.unit.regeneration_stream[0.0].flow_mass_phase_comp["Liq", "Ca_2+"]
        ] = 0.049885
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "Ca_2+"]
        ] = -1.24713

        return m


def build_freundlich():

    target_ion = "Cl_-"

    ion_props = {
        "solute_list": [target_ion],
        "diffusivity_data": {("Liq", target_ion): 1e-9},
        "mw_data": {"H2O": 0.018, target_ion: 35.45e-3},
        "charge": {target_ion: -1},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(**ion_props)
    ix_config = {
        "property_package": m.fs.properties,
        "target_ion": target_ion,
        "isotherm": "freundlich",
        "regenerant": "NaOH",
        "hazardous_waste": True,
    }
    m.fs.unit = ix = IonExchange0D(**ix_config)

    ix.process_flow.properties_in.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 0.5,
            ("conc_mass_phase_comp", ("Liq", target_ion)): 1e-6,
            ("pressure", None): 101325,
            ("temperature", None): 298,
        },
        hold_state=True,
    )

    ix.process_flow.properties_in[0].flow_mass_phase_comp[...]
    ix.process_flow.properties_out[0].flow_mass_phase_comp[...]
    ix.regeneration_stream[0].flow_mass_phase_comp[...]

    ix.freundlich_n.fix(1.2)
    ix.bv_50.fix(20000)
    ix.bv.fix(18000)
    ix.resin_bulk_dens.fix(0.72)
    ix.bed_porosity.fix()
    ix.vel_bed.fix(6.15e-3)
    ix.resin_diam.fix(6.75e-4)
    ix.number_columns.fix(16)
    ix.service_flow_rate.fix(15)
    ix.c_norm.fix(0.25)

    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e-4, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e6, index=("Liq", "Cl_-")
    )
    iscale.calculate_scaling_factors(m)

    return m


class TestIXFreundlich(UnitTestHarness):
    def configure(self):
        m = build_freundlich()

        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance
        self.unit_solutions[m.fs.unit.col_height] = 3.160
        self.unit_solutions[m.fs.unit.col_diam] = 2.5435
        self.unit_solutions[m.fs.unit.N_Sh["Cl_-"]] = 24.083
        self.unit_solutions[m.fs.unit.N_Pe_particle] = 0.09901
        self.unit_solutions[m.fs.unit.N_Pe_bed] = 216.51
        self.unit_solutions[m.fs.unit.tb_traps[0]] = 0
        self.unit_solutions[m.fs.unit.tb_traps[1]] = 3344557
        self.unit_solutions[m.fs.unit.tb_traps[2]] = 3825939
        self.unit_solutions[m.fs.unit.tb_traps[3]] = 4034117
        self.unit_solutions[m.fs.unit.tb_traps[4]] = 4188551
        self.unit_solutions[m.fs.unit.tb_traps[5]] = 4320000
        self.unit_solutions[m.fs.unit.traps[1]] = 0.0038710157
        self.unit_solutions[m.fs.unit.traps[2]] = 0.0044572367
        self.unit_solutions[m.fs.unit.traps[3]] = 0.0048189406
        self.unit_solutions[m.fs.unit.traps[4]] = 0.0057197676
        self.unit_solutions[m.fs.unit.traps[5]] = 0.0066941563
        self.unit_solutions[m.fs.unit.c_norm_avg["Cl_-"]] = 0.02556
        self.unit_solutions[m.fs.unit.mass_transfer_coeff] = 0.159346
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "Cl_-"]
        ] = -1.3744e-05
        self.unit_solutions[
            m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 500
        self.unit_solutions[
            m.fs.unit.process_flow.properties_out[0.0].flow_mass_phase_comp[
                "Liq", "H2O"
            ]
        ] = 500
        self.unit_solutions[
            m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                "Liq", "Cl_-"
            ]
        ] = 5e-07
        self.unit_solutions[
            m.fs.unit.process_flow.properties_out[0.0].flow_mass_phase_comp[
                "Liq", "Cl_-"
            ]
        ] = 1.278e-08
        self.unit_solutions[
            m.fs.unit.regeneration_stream[0.0].flow_mass_phase_comp["Liq", "Cl_-"]
        ] = 4.872e-07

        return m


def build_inert():
    target_ion = "Cl_-"
    inert = "Ca_2+"

    ion_props = {
        "solute_list": [target_ion, inert],
        "diffusivity_data": {("Liq", target_ion): 1e-9, ("Liq", inert): 9.2e-10},
        "mw_data": {"H2O": 0.018, target_ion: 35.45e-3, inert: 40e-3},
        "charge": {target_ion: -1, inert: 2},
    }

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(**ion_props)
    ix_config = {
        "property_package": m.fs.properties,
        "target_ion": target_ion,
        "regenerant": "single_use",
        "isotherm": "freundlich",
    }
    m.fs.unit = ix = IonExchange0D(**ix_config)

    ix.process_flow.properties_in.calculate_state(
        var_args={
            ("flow_vol_phase", "Liq"): 0.5,
            ("conc_mass_phase_comp", ("Liq", target_ion)): 1e-6,
            ("conc_mass_phase_comp", ("Liq", inert)): 1e-7,
            ("pressure", None): 101325,
            ("temperature", None): 298,
        },
        hold_state=True,
    )

    ix.process_flow.properties_in[0].flow_mass_phase_comp[...]
    ix.process_flow.properties_out[0].flow_mass_phase_comp[...]
    ix.regeneration_stream[0].flow_mass_phase_comp[...]

    ix.freundlich_n.fix(1.2)
    ix.bv_50.fix(20000)
    ix.bv.fix(18000)
    ix.resin_bulk_dens.fix(0.72)
    ix.bed_porosity.fix()
    ix.vel_bed.fix(6.15e-3)
    ix.resin_diam.fix(6.75e-4)
    ix.number_columns.fix(16)
    ix.service_flow_rate.fix(15)
    ix.c_norm.fix(0.25)
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e-4, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e5, index=("Liq", "Cl_-")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e5, index=("Liq", "Ca_2+")
    )
    iscale.calculate_scaling_factors(m)

    return m


class TestIXInert(UnitTestHarness):
    def configure(self):
        m = build_inert()

        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance

        self.unit_solutions[m.fs.unit.number_columns] = 16
        self.unit_solutions[m.fs.unit.col_height] = 3.160
        self.unit_solutions[m.fs.unit.col_diam] = 2.5435
        self.unit_solutions[m.fs.unit.N_Sh["Cl_-"]] = 24.083
        self.unit_solutions[m.fs.unit.N_Pe_particle] = 0.09901
        self.unit_solutions[m.fs.unit.N_Pe_bed] = 216.51
        self.unit_solutions[m.fs.unit.tb_traps[1]] = 3344557
        self.unit_solutions[m.fs.unit.tb_traps[2]] = 3825939
        self.unit_solutions[m.fs.unit.tb_traps[3]] = 4034117
        self.unit_solutions[m.fs.unit.tb_traps[4]] = 4188551
        self.unit_solutions[m.fs.unit.tb_traps[5]] = 4320000
        self.unit_solutions[m.fs.unit.c_traps[2]] = 0.07
        self.unit_solutions[m.fs.unit.c_traps[3]] = 0.13
        self.unit_solutions[m.fs.unit.c_traps[4]] = 0.19
        self.unit_solutions[m.fs.unit.c_traps[5]] = 0.25
        self.unit_solutions[m.fs.unit.traps[4]] = 0.0057197676
        self.unit_solutions[m.fs.unit.traps[5]] = 0.0066941563
        self.unit_solutions[m.fs.unit.c_norm_avg["Cl_-"]] = 0.02556
        self.unit_solutions[m.fs.unit.mass_transfer_coeff] = 0.159346
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "Cl_-"]
        ] = -1.37435e-05
        self.unit_solutions[
            m.fs.unit.process_flow.mass_transfer_term[0.0, "Liq", "Ca_2+"]
        ] = 0
        self.unit_solutions[
            m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 500
        self.unit_solutions[
            m.fs.unit.process_flow.properties_out[0.0].flow_mass_phase_comp[
                "Liq", "H2O"
            ]
        ] = 500
        self.unit_solutions[
            m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                "Liq", "Cl_-"
            ]
        ] = 5e-07
        self.unit_solutions[
            m.fs.unit.process_flow.properties_in[0.0].flow_mass_phase_comp[
                "Liq", "Ca_2+"
            ]
        ] = 5e-08
        self.unit_solutions[
            m.fs.unit.process_flow.properties_out[0.0].flow_mass_phase_comp[
                "Liq", "Cl_-"
            ]
        ] = 1.2781e-08
        self.unit_solutions[
            m.fs.unit.process_flow.properties_out[0.0].flow_mass_phase_comp[
                "Liq", "Ca_2+"
            ]
        ] = 5e-08
        self.unit_solutions[
            m.fs.unit.regeneration_stream[0.0].flow_mass_phase_comp["Liq", "Cl_-"]
        ] = 4.872e-07

        return m
