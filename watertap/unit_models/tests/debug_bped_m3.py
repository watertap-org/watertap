"""
Standalone debug script for bped_m[3]:
  - Constant Current, 20A, n=50 cell triplets, salt_calculation=True, has_catalyst=True
  - Run:         python debug_bped_m3.py
  - Interactive: python -i debug_bped_m3.py
"""

import idaes.core.util.scaling as iscale
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.environ import ConcreteModel, value

from watertap.unit_models.Biploar_and_Electrodialysis_1D_nmsu import (
    Bipolar_and_Electrodialysis1D,
    ElectricalOperationMode,
    LimitingCurrentDensitybpemMethod,
)
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock

# ── 1. Build ──────────────────────────────────────────────────────────────────
print("=" * 60)
print("Building bped_m[3]: Constant Current, salt_calc=True")
print("=" * 60)

m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
m.fs.properties = MCASParameterBlock(
    solute_list=["Na_+", "Cl_-", "H_+", "OH_-"],
    mw_data={"Na_+": 23e-3, "Cl_-": 35.5e-3, "H_+": 1e-3, "OH_-": 17.0e-3},
    elec_mobility_data={
        ("Liq", "Na_+"): 5.19e-8,
        ("Liq", "Cl_-"): 7.92e-8,
        ("Liq", "H_+"): 36.23e-8,
        ("Liq", "OH_-"): 20.64e-8,
    },
    charge={"Na_+": 1, "Cl_-": -1, "H_+": 1, "OH_-": -1},
    diffusivity_data={
        ("Liq", "Na_+"): 1.33e-9,
        ("Liq", "Cl_-"): 2.03e-9,
        ("Liq", "H_+"): 9.31e-9,
        ("Liq", "OH_-"): 5.27e-9,
    },
)
m.fs.unit = Bipolar_and_Electrodialysis1D(
    property_package=m.fs.properties,
    operation_mode=ElectricalOperationMode.Constant_Current,
    limiting_current_density_method_bpem=LimitingCurrentDensitybpemMethod.Empirical,
    has_catalyst=True,
    salt_calculation=True,
)
u = m.fs.unit

# ── 2. Shared fixture settings ────────────────────────────────────────────────
print("Applying shared fixture settings...")
u.diffus_mass.fix((2.03 + 1.96) * 10**-9 / 2)
u.membrane_fixed_charge.fix(5e3)
u.salt_conc_aem_ref.fix(2000)
u.salt_conc_cem_ref.fix(2000)
u.salt_conc_dilu_ref.fix(2000)
u.membrane_areal_resistance_coef_0.set_value((1.89e-4 + 1.77e-4) / 2)
u.membrane_areal_resistance_coef_1.set_value(0)
u.conc_water.fix(50 * 1e3)
u.k2_zero.fix(2 * 10**-6)
u.relative_permittivity.fix(30)
u.membrane_fixed_catalyst_cem.fix(5e3)
u.membrane_fixed_catalyst_aem.fix(5e3)
u.k_a.fix(447)
u.k_b.fix(5e4)

u.ion_trans_number_membrane["bpem", "Na_+"].fix(0)
u.ion_trans_number_membrane["bpem", "Cl_-"].fix(0)
u.ion_trans_number_membrane["bpem", "H_+"].fix(1)
u.ion_trans_number_membrane["bpem", "OH_-"].fix(1)
u.ion_trans_number_membrane["cem", "Na_+"].fix(0.94)
u.ion_trans_number_membrane["cem", "Cl_-"].fix(0)
u.ion_trans_number_membrane["cem", "H_+"].fix(0.03)
u.ion_trans_number_membrane["cem", "OH_-"].fix(0.03)
u.ion_trans_number_membrane["aem", "Na_+"].fix(0)
u.ion_trans_number_membrane["aem", "Cl_-"].fix(0.94)
u.ion_trans_number_membrane["aem", "H_+"].fix(0.03)
u.ion_trans_number_membrane["aem", "OH_-"].fix(0.03)

u.solute_diffusivity_membrane["bpem", "Na_+"].fix(0)
u.solute_diffusivity_membrane["bpem", "Cl_-"].fix(0)
u.solute_diffusivity_membrane["bpem", "H_+"].fix(0)
u.solute_diffusivity_membrane["bpem", "OH_-"].fix(0)
u.solute_diffusivity_membrane["cem", "Na_+"].fix((1.8e-10 + 1.25e-10) / 2)
u.solute_diffusivity_membrane["cem", "Cl_-"].fix((1.8e-10 + 1.25e-10) / 2)
u.solute_diffusivity_membrane["cem", "H_+"].fix(0)
u.solute_diffusivity_membrane["cem", "OH_-"].fix(0)
u.solute_diffusivity_membrane["aem", "Na_+"].fix((1.8e-10 + 1.25e-10) / 2)
u.solute_diffusivity_membrane["aem", "Cl_-"].fix((1.8e-10 + 1.25e-10) / 2)
u.solute_diffusivity_membrane["aem", "H_+"].fix(0)
u.solute_diffusivity_membrane["aem", "OH_-"].fix(0)

u.inlet_basate.pressure.fix(101325)
u.inlet_basate.temperature.fix(298.15)
u.inlet_acidate.pressure.fix(101325)
u.inlet_acidate.temperature.fix(298.15)
u.inlet_diluate.pressure.fix(101325)
u.inlet_diluate.temperature.fix(298.15)
u.spacer_porosity.fix(1)

u.inlet_basate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-2)
u.inlet_basate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(0)
u.inlet_basate.flow_mol_phase_comp[0, "Liq", "H_+"].fix(0)
u.inlet_basate.flow_mol_phase_comp[0, "Liq", "OH_-"].fix(7.38e-2)
u.inlet_acidate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(0)
u.inlet_acidate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-2)
u.inlet_acidate.flow_mol_phase_comp[0, "Liq", "H_+"].fix(7.38e-2)
u.inlet_acidate.flow_mol_phase_comp[0, "Liq", "OH_-"].fix(0)
u.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-2)
u.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-2)
u.inlet_diluate.flow_mol_phase_comp[0, "Liq", "H_+"].fix(0)
u.inlet_diluate.flow_mol_phase_comp[0, "Liq", "OH_-"].fix(0)
u.inlet_basate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
u.inlet_acidate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
u.inlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)

u.shadow_factor.fix(1)
u.water_trans_number_membrane["cem"].fix((5.8 + 4.3) / 2)
u.water_permeability_membrane["cem"].fix((2.16e-14 + 1.75e-14) / 2)
u.water_trans_number_membrane["aem"].fix((5.8 + 4.3) / 2)
u.water_permeability_membrane["aem"].fix((2.16e-14 + 1.75e-14) / 2)
u.water_trans_number_membrane["bpem"].fix((5.8 + 4.3) / 2)
u.water_permeability_membrane["bpem"].fix((2.16e-14 + 1.75e-14) / 2)
u.electrodes_resistance.fix(0)
u.current_utilization.fix(1)
u.channel_height["diluate"].fix(2.7e-4)
u.channel_height["basate"].fix(2.7e-4)
u.channel_height["acidate"].fix(2.7e-4)
u.cell_width.fix(0.1)
u.cell_length.fix(0.79)
u.membrane_thickness["bpem"].fix(8e-4)
u.membrane_thickness["aem"].fix(4e-4)
u.membrane_thickness["cem"].fix(4e-4)
u.electrical_stage_num.fix(1)

m.fs.properties.set_default_scaling("flow_mol_phase_comp", 1e3, index=("Liq", "Na_+"))
m.fs.properties.set_default_scaling("flow_mol_phase_comp", 1e3, index=("Liq", "Cl_-"))
m.fs.properties.set_default_scaling("flow_mol_phase_comp", 1e3, index=("Liq", "H_+"))
m.fs.properties.set_default_scaling("flow_mol_phase_comp", 1e3, index=("Liq", "OH_-"))
m.fs.properties.set_default_scaling("flow_mol_phase_comp", 1e2, index=("Liq", "H2O"))
u.cell_triplet_num.fix(10)
iscale.set_scaling_factor(u.k_a, 1e-2)
iscale.set_scaling_factor(u.k_b, 1e-4)
iscale.set_scaling_factor(u.voltage_x, 1e-1)
iscale.set_scaling_factor(u.flux_splitting, 1e4)
iscale.set_scaling_factor(u.current_density_x, 1e-3)

# ── 3. bped_m[3]-specific ─────────────────────────────────────────────────────
print("Applying bped_m[3]-specific settings...")
u.current_applied.fix(2e1)
u.cell_triplet_num.fix(50)
iscale.set_scaling_factor(u.voltage_x, 1e-0)
m.fs.properties.set_default_scaling("flow_mol_phase_comp", 1e4, index=("Liq", "H_+"))
m.fs.properties.set_default_scaling("flow_mol_phase_comp", 1e4, index=("Liq", "OH_-"))
iscale.calculate_scaling_factors(m.fs)

# ── 4. Pre-init inspection ────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("PRE-INITIALIZATION STATE")
print("=" * 60)

print("DOF =", degrees_of_freedom(m))
print("hasattr current_applied:", hasattr(u, "current_applied"))
print("hasattr voltage_applied:", hasattr(u, "voltage_applied"))

t0 = m.fs.time.first()
print("current_applied[t0].fixed:", u.current_applied[t0].fixed)
print("current_applied[t0].value:", value(u.current_applied[t0]))

i_init = value(u.current_density_x[t0, 0.5])
ca_val = value(u.current_applied[t0])
sf_val = value(u.shadow_factor)
cw_val = value(u.cell_width)
cl_val = value(u.cell_length)
exp_i = ca_val / (cw_val * sf_val * cl_val)
print("current_density_x before init:", round(i_init, 4), "A/m2")
print("Expected current_density_x:   ", round(exp_i, 2), "A/m2")

na_inlet = value(u.diluate.properties[(t0, 0.0)].flow_mol_phase_comp["Liq", "Na_+"])
na_09 = value(u.diluate.properties[(t0, 0.9)].flow_mol_phase_comp["Liq", "Na_+"])
print("Na+ flow x=0.0 (inlet):", round(na_inlet, 6), "mol/s")
print("Na+ flow x=0.9 before: ", round(na_09, 6), "mol/s")

# ── 5. Run initialize ─────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("RUNNING initialize() -- watch for DEBUG prints")
print("=" * 60)
u.initialize(fail_on_warning=False)

# ── 6. Post-init inspection ───────────────────────────────────────────────────
print("\n" + "=" * 60)
print("POST-INITIALIZATION STATE")
print("=" * 60)

print("DOF =", degrees_of_freedom(m))
print(
    "current_density_x after init:",
    round(value(u.current_density_x[t0, 0.5]), 4),
    "A/m2",
)

na_09_after = value(u.diluate.properties[(t0, 0.9)].flow_mol_phase_comp["Liq", "Na_+"])
conc_na_09 = value(u.diluate.properties[(t0, 0.9)].conc_mol_phase_comp["Liq", "Na_+"])
print("Na+ flow x=0.9 after:  ", round(na_09_after, 6), "mol/s")
print("Na+ conc x=0.9 after:  ", round(conc_na_09, 4), "mol/m3")

print("\nResiduals for eq_conc_mol_phase_comp[Na+] at outlet nodes:")
for x in [0.8, 0.9, 1.0]:
    prop = u.diluate.properties[(t0, x)]
    conc = value(prop.conc_mol_phase_comp["Liq", "Na_+"])
    flow = value(prop.flow_mol_phase_comp["Liq", "Na_+"])
    vol = value(prop.flow_vol_phase["Liq"])
    resid = abs(conc - flow / vol)
    print(
        "  x="
        + str(x)
        + "  conc="
        + str(round(conc, 4))
        + "  flow/vol="
        + str(round(flow / vol, 4))
        + "  residual="
        + "{:.2e}".format(resid)
    )

print("\nDone. 'm' and 'u' available for interactive inspection.")
