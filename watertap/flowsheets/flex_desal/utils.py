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

"""
This module contains some utility functions
"""
from pyomo.environ import SolverFactory
from pyomo.common.dependencies import attempt_import

try:
    from gurobipy import nlfunc

    gurobipy_available = True
except:
    gurobipy_available = False


def get_gurobi_solver_model(m, mip_gap=0.01, time_limit=3600, tee=True):
    """
    Returns a Pyomo SolverFactory object that is compatible with Gurobi.
    This function is needed only when the RO recovery is a variable.
    """
    if not gurobipy_available:
        pass
    else:
        solver = SolverFactory("gurobi_persistent")
        solver.options["MIPGap"] = mip_gap
        solver.options["TimeLimit"] = time_limit
        solver.options["OutputFlag"] = int(tee)

        if (
            not m.period[1, 1]
            .reverse_osmosis.ro_skid[1]
            .calculate_energy_intensity.active
        ):
            # If the nonlinear constraint is not active, then return the solver
            # object directly
            solver.set_instance(m)
            return solver

        if m.params.surrogate_type == "quadratic_surrogate":
            # Model is quadratic, so Pyomo's writer can handle it.
            solver.set_instance(m)
            return solver

        # Nonlinear constraint is present.
        # Step 1: Deactivate the nonlinear constraint
        for p in m.period:
            for skid in m.period[p].reverse_osmosis.set_ro_skids:
                m.period[p].reverse_osmosis.ro_skid[
                    skid
                ].calculate_energy_intensity.deactivate()

        # pylint: disable = protected-access
        # Step 2: Build the gurobipy model
        solver.set_instance(m)
        gm = solver._solver_model  # Gurobipy model
        pm_to_gm = solver._pyomo_var_to_solver_var_map

        # Step 3: Add the nonlinear constraint
        coeffs = m.period[1, 1].reverse_osmosis.ro_skid[1].coeffs
        for p in m.period:
            for skid in m.period[p].reverse_osmosis.set_ro_skids:
                ro_skid = m.period[p].reverse_osmosis.ro_skid[skid]
                recovery_var = pm_to_gm[ro_skid.recovery]
                gm.addConstr(
                    (
                        pm_to_gm[ro_skid.energy_intensity]
                        == coeffs["a"] * nlfunc.exp(-coeffs["b"] * recovery_var)
                        + coeffs["c"] * recovery_var * recovery_var
                        + coeffs["d"]
                    ),
                    name=f"ro_energy_intensity_{p[0]}_{p[1]}_{skid}",
                )

    return solver


def fix_recovery(m, recovery):
    """Modifies the model for the fixed recovery case"""
    # Compute the energy intensity
    ro_skid = m.period[1, 1].reverse_osmosis.ro_skid[1]
    energy_intensity = m.params.ro.get_energy_intensity(recovery)

    for p in m.period:
        for skid in m.period[p].reverse_osmosis.set_ro_skids:
            ro_skid = m.period[p].reverse_osmosis.ro_skid[skid]
            ro_skid.recovery.fix(recovery)
            ro_skid.energy_intensity.fix(energy_intensity)
            ro_skid.calculate_energy_intensity.deactivate()


def update_recovery_bounds(m, lb, ub):
    """Updates the bounds on the recovery variable"""
    ro_skid = m.period[1, 1].reverse_osmosis.ro_skid[1]
    ei_lb, ei_ub = m.params.ro.get_energy_intensity_bounds(lb, ub)

    for p in m.period:
        for skid in m.period[p].reverse_osmosis.set_ro_skids:
            ro_skid = m.period[p].reverse_osmosis.ro_skid[skid]
            ro_skid.recovery.setlb(lb)
            ro_skid.recovery.setub(ub)
            ro_skid.energy_intensity.setlb(ei_lb)
            ro_skid.energy_intensity.setub(ei_ub)


def get_baseline_model(m):
    """Returns a baseline model from the given model"""
    bm = m.clone()

    # Ensure that the pretreatment unit is always on, and no leakage from it
    bm.fix_operation_var("pretreatment.op_mode", 1)
    bm.fix_operation_var("pretreatment.recovery", 1)

    # Ensure that the first three skids all always on
    for skid in [1, 2, 3]:
        bm.fix_operation_var(f"reverse_osmosis.ro_skid[{skid}].op_mode", 1)

    if hasattr(m.period[1, 1], "battery"):
        # If the battery model exists, ensure that there is no charging
        # and discharging
        bm.fix_operation_var("battery.power_charge", 0)
        bm.fix_operation_var("battery.power_discharge", 0)

    return bm
