#!/usr/bin/env python
"""
Quick smoke tests for the current BSM2_genericNP_no_bioP flowsheet using the
correct costing-first initialization pattern. Runs a handful of representative
parameter combinations similar to the sweep.
"""
import logging

import pyomo.environ as pyo

import watertap.flowsheets.electroNP.BSM2_genericNP_no_bioP as genericNP_flowsheet


logging.getLogger("idaes.core.util.scaling").disabled = True


def _build_and_solve(nh4_removal, p_removal, energy_intensity, mgcl2_dosage=0.388):
    # Build
    m = genericNP_flowsheet.build_flowsheet(has_genericNP=True, basis="mass")

    # Add costing before setting conditions (matching build_model pattern)
    genericNP_flowsheet.add_costing(m)
    for c in m.fs.component_objects(pyo.Constraint, descend_into=True):
        if "capital_cost" in c.name:
            c.deactivate()

    n_to_p_ratio = nh4_removal / p_removal
    genericNP_flowsheet.set_operating_conditions(
        m,
        p_removal=p_removal,
        n_to_p_ratio=n_to_p_ratio,
        energy_intensity=energy_intensity,
        mgcl2_dosage=mgcl2_dosage,
    )

    # Deactivate pressure constraints before init
    for mx in m.fs.mixers:
        mx.pressure_equality_constraints[0.0, 2].deactivate()
    m.fs.MX3.pressure_equality_constraints[0.0, 2].deactivate()
    m.fs.MX3.pressure_equality_constraints[0.0, 3].deactivate()

    # Initialize with costing present
    genericNP_flowsheet.initialize_system(m, has_genericNP=True)
    # Re-deactivate pressure constraints after init
    for mx in m.fs.mixers:
        mx.pressure_equality_constraints[0.0, 2].deactivate()
    m.fs.MX3.pressure_equality_constraints[0.0, 2].deactivate()
    m.fs.MX3.pressure_equality_constraints[0.0, 3].deactivate()

    solver = genericNP_flowsheet.get_solver()
    solver.options["max_iter"] = 5000
    results = solver.solve(m, tee=False)
    pyo.assert_optimal_termination(results)

    return m, results


def test_single_case(case_name, nh4_removal, p_removal=0.95, energy_intensity=0.044):
    print(f"\n{'='*70}")
    print(f"Testing Case: {case_name}")
    print(f"  NH4 removal: {nh4_removal:.2f}")
    print(f"  P removal: {p_removal:.2f}")
    print(f"  Energy intensity: {energy_intensity:.4f} kWh/kg")
    print(f"{'='*70}")

    try:
        m, results = _build_and_solve(nh4_removal, p_removal, energy_intensity)

        lcow = pyo.value(m.fs.costing.LCOW)
        treated_nh4 = pyo.value(m.fs.genericNP.treated.conc_mass_comp[0, "S_NH4"])
        treated_po4 = pyo.value(m.fs.genericNP.treated.conc_mass_comp[0, "S_PO4"])

        print(f"\n  Results:")
        print(f"    LCOW: ${lcow:.4f}/m³")
        print(f"    Treated NH4: {treated_nh4:.6f} kg/m³")
        print(f"    Treated PO4: {treated_po4:.6f} kg/m³")

        return True

    except Exception as e:
        print(f"FAILED with exception: {e}")
        import traceback

        traceback.print_exc()
        return False


def main():

    # Representative sweep points (keep short for quick smoke tests)
    test_cases = [
        ("Baseline", 0.285, 0.95, 0.044),
        ("Lower NH4 removal", 0.15, 0.95, 0.044),
        ("Higher NH4 removal", 0.60, 0.95, 0.044),
    ]

    results = []
    for case_name, nh4_removal, p_removal, energy_intensity in test_cases:
        success = test_single_case(case_name, nh4_removal, p_removal, energy_intensity)
        results.append((case_name, success))

    successes = sum(1 for _, success in results if success)
    total = len(results)

    for case_name, success in results:
        status = "PASS" if success else "FAIL"
        print(f"{status}: {case_name}")

    return successes == total


if __name__ == "__main__":
    import sys

    success = main()
    sys.exit(0 if success else 1)
