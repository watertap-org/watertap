#!/usr/bin/env python
"""Diagnose degrees of freedom for parameter sweep in electroNP flowsheets."""

import logging

logging.getLogger("idaes.core.util.scaling").disabled = True
# logging.getLogger("idaes.init").disabled = True
# logging.getLogger("idaes.solve").disabled = True

from idaes.core.util.model_statistics import degrees_of_freedom
from watertap.flowsheets.electroNP import BSM2_genericNP_no_bioP
from watertap.flowsheets.electroNP import BSM2_electroNP_no_bioP


def diagnose_flowsheet_dof(flowsheet_module, flowsheet_name, has_nutrient_removal=True):
    """Analyze DOF of a flowsheet at different stages.

    Parameters
    ----------
    flowsheet_module : module
        The flowsheet module (e.g., BSM2_genericNP_no_bioP)
    flowsheet_name : str
        Name of the flowsheet for display
    has_nutrient_removal : bool
        Whether the flowsheet has removal_factors (genericNP) vs other intensity vars

    Returns
    -------
    dict
        Results containing DOF at each stage
    """
    results = {
        "name": flowsheet_name,
        "dof_after_build": None,
        "dof_after_set_ops": None,
        "dof_after_init": None,
        "dof_after_deactivate": None,
        "dof_after_costing": None,
        "variables_to_sweep": [],
    }

    # Build flowsheet
    if flowsheet_name == "BSM2_genericNP_no_bioP":
        m = flowsheet_module.build_flowsheet(
            has_genericNP=has_nutrient_removal, basis="mass"
        )
    elif flowsheet_name == "BSM2_electroNP_no_bioP":
        m = flowsheet_module.build_flowsheet(has_electroNP=has_nutrient_removal)

    results["dof_after_build"] = degrees_of_freedom(m)

    # Set operating conditions (fixes variables)
    flowsheet_module.set_operating_conditions(m)
    results["dof_after_set_ops"] = degrees_of_freedom(m)

    # Initialize system
    if flowsheet_name == "BSM2_genericNP_no_bioP":
        flowsheet_module.initialize_system(m, has_genericNP=has_nutrient_removal)
    elif flowsheet_name == "BSM2_electroNP_no_bioP":
        flowsheet_module.initialize_system(m, has_electroNP=has_nutrient_removal)
    results["dof_after_init"] = degrees_of_freedom(m)

    # Deactivate pressure equality constraints (BSM2-specific)
    for mx in m.fs.mixers:
        if hasattr(mx, "pressure_equality_constraints"):
            mx.pressure_equality_constraints[0.0, 2].deactivate()
    m.fs.MX3.pressure_equality_constraints[0.0, 2].deactivate()
    m.fs.MX3.pressure_equality_constraints[0.0, 3].deactivate()

    results["dof_after_deactivate"] = degrees_of_freedom(m)

    # Add costing
    if flowsheet_name == "BSM2_genericNP_no_bioP":
        flowsheet_module.add_costing(m)
        results["dof_after_costing"] = degrees_of_freedom(m)

    # Identify which variables could be swept
    if has_nutrient_removal and flowsheet_name == "BSM2_genericNP_no_bioP":
        # GenericNP has removal_factors
        unit = m.fs.genericNP
        for comp in unit.removal_factors:
            var_name = f"removal_factors[{comp}]"
            unit.removal_factors[comp].unfix()
            dof_after_unfix = degrees_of_freedom(m)
            results["variables_to_sweep"].append(
                {
                    "name": var_name,
                    "dof_change": dof_after_unfix - results["dof_after_costing"],
                }
            )
            unit.removal_factors[comp].fix()

    return results


if __name__ == "__main__":
    # Analyze genericNP
    print("\n1. GenericNP Flowsheet (BSM2_genericNP_no_bioP)")
    results_generic = diagnose_flowsheet_dof(
        BSM2_genericNP_no_bioP, "BSM2_genericNP_no_bioP", has_nutrient_removal=True
    )
    print(f"   DOF after build:            {results_generic['dof_after_build']:+3d}")
    print(f"   DOF after set_conditions:   {results_generic['dof_after_set_ops']:+3d}")
    print(f"   DOF after init:             {results_generic['dof_after_init']:+3d}")
    print(
        f"   DOF after deactivate:       {results_generic['dof_after_deactivate']:+3d}"
    )
    print(f"   DOF after add_costing:      {results_generic['dof_after_costing']:+3d}")
    if results_generic["variables_to_sweep"]:
        print("   Sweepable variables:")
        for var in results_generic["variables_to_sweep"]:
            print(f"     - {var['name']}: DOF change {var['dof_change']:+d}")

    # Analyze electroNP
    print("\n2. ElectroNP Flowsheet (BSM2_electroNP_no_bioP)")
    results_electro = diagnose_flowsheet_dof(
        BSM2_electroNP_no_bioP, "BSM2_electroNP_no_bioP", has_nutrient_removal=False
    )
    print(f"   DOF after build:            {results_electro['dof_after_build']:+3d}")
    print(f"   DOF after set_conditions:   {results_electro['dof_after_set_ops']:+3d}")
    print(f"   DOF after init:             {results_electro['dof_after_init']:+3d}")
    print(
        f"   DOF after deactivate:       {results_electro['dof_after_deactivate']:+3d}"
    )
    if results_electro["variables_to_sweep"]:
        print("   Sweepable variables:")
        for var in results_electro["variables_to_sweep"]:
            print(f"     - {var['name']}: DOF change {var['dof_change']:+d}")
