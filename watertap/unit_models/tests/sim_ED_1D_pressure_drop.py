import pytest
import re
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.unit_models.electrodialysis_1D import (
    ElectricalOperationMode,
    Electrodialysis1D,
    PressureDropMethod,
    FrictionFactorMethod,
    HydraulicDiameterMethod,
    LimitingCurrentDensityMethod,
)
from watertap.costing import WaterTAPCosting
from pyomo.environ import (
    ConcreteModel,
    assert_optimal_termination,
    value,
    Set,
    Var,
    Constraint,
)
from idaes.core import (
    FlowsheetBlock,
    EnergyBalanceType,
    MaterialBalanceType,
    MomentumBalanceType,
    UnitModelCostingBlock,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
import idaes.core.util.scaling as iscale
from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog

solver = get_solver()


def main():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    ion_dict = {
        "solute_list": ["Na_+", "Cl_-"],
        "mw_data": {"H2O": 18e-3, "Na_+": 23e-3, "Cl_-": 35.5e-3},
        "elec_mobility_data": {("Liq", "Na_+"): 5.19e-8, ("Liq", "Cl_-"): 7.92e-8},
        "charge": {"Na_+": 1, "Cl_-": -1},
    }
    m.fs.properties = MCASParameterBlock(**ion_dict)
    m.fs.unit = Electrodialysis1D(
        property_package=m.fs.properties,
        operation_mode=ElectricalOperationMode.Constant_Voltage,
        # has_pressure_change=False,
        # limiting_current_density_method=LimitingCurrentDensityMethod.Theoretical,
        has_pressure_change=True,
        has_nonohmic_potential_membrane=False,
        has_Nernst_diffusion_layer=False,
        pressure_drop_method=PressureDropMethod.Darcy_Weisbach,
        # pressure_drop_method=PressureDropMethod.experimental,
        # pressure_drop_method=PressureDropMethod.none,
        friction_factor_method=FrictionFactorMethod.Gurreri,
        # friction_factor_method=FrictionFactorMethod.Kuroda,
        # friction_factor_method=FrictionFactorMethod.fixed,
        # hydraulic_diameter_method=HydraulicDiameterMethod.conventional,
        # hydraulic_diameter_method=HydraulicDiameterMethod.spacer_specific_area_known,
        hydraulic_diameter_method=HydraulicDiameterMethod.fixed,
    )

    m.fs.unit.inlet_diluate.pressure.fix(501325)
    m.fs.unit.inlet_diluate.temperature.fix(298.15)
    m.fs.unit.inlet_concentrate.pressure.fix(501325)
    m.fs.unit.inlet_concentrate.temperature.fix(298.15)
    m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(17.875)
    m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(5.56e-2)
    m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(5.56e-2)
    m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(17.875)
    m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(5.56e-2)
    m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(5.56e-2)
    m.fs.unit.water_trans_number_membrane["cem"].fix(5.8)
    m.fs.unit.water_trans_number_membrane["aem"].fix(4.3)
    m.fs.unit.water_permeability_membrane["cem"].fix(2.16e-14)
    m.fs.unit.water_permeability_membrane["aem"].fix(1.75e-14)
    m.fs.unit.electrodes_resistance.fix(0)
    m.fs.unit.cell_pair_num.fix(56)
    m.fs.unit.current_utilization.fix(1)
    m.fs.unit.channel_height.fix(7.1e-4)
    m.fs.unit.membrane_areal_resistance["cem"].fix(1.89e-4)
    m.fs.unit.membrane_areal_resistance["aem"].fix(1.77e-4)
    m.fs.unit.cell_width.fix(0.197)
    m.fs.unit.cell_length.fix(1.68)
    m.fs.unit.membrane_thickness["aem"].fix(1.3e-4)
    m.fs.unit.membrane_thickness["cem"].fix(1.3e-4)
    m.fs.unit.solute_diffusivity_membrane["cem", "Na_+"].fix(1.8e-10)
    m.fs.unit.solute_diffusivity_membrane["aem", "Na_+"].fix(1.25e-10)
    m.fs.unit.solute_diffusivity_membrane["cem", "Cl_-"].fix(1.8e-10)
    m.fs.unit.solute_diffusivity_membrane["aem", "Cl_-"].fix(1.25e-10)
    m.fs.unit.ion_trans_number_membrane["cem", "Na_+"].fix(1)
    m.fs.unit.ion_trans_number_membrane["aem", "Na_+"].fix(0)
    m.fs.unit.ion_trans_number_membrane["cem", "Cl_-"].fix(0)
    m.fs.unit.ion_trans_number_membrane["aem", "Cl_-"].fix(1)
    m.fs.unit.voltage_applied.fix(40)
    m.fs.unit.spacer_porosity.fix(0.83)
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 0.1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e2, index=("Liq", "Na_+")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e2, index=("Liq", "Cl_-")
    )
    iscale.set_scaling_factor(m.fs.unit.cell_width, 5)
    iscale.set_scaling_factor(m.fs.unit.cell_length, 1)
    iscale.set_scaling_factor(m.fs.unit.cell_pair_num, 0.1)

    # fix additional variables for specific pressure drop methods
    m.fs.unit.diffus_mass.fix(1.6e-9) if hasattr(m.fs.unit, "diffus_mass") else 0
    # m.fs.unit.spacer_specific_area.fix(10700) if hasattr(
    #     m.fs.unit, "spacer_specific_area"
    # ) else 0
    # m.display()

    m.fs.unit.hydraulic_diameter.fix(1.5e-3) if hasattr(
        m.fs.unit, "hydraulic_diameter"
    ) else 0
    # m.fs.unit.friction_factor.fix(20) if hasattr(
    #     m.fs.unit, "friction_factor"
    # ) else 0
    # m.fs.unit.pressure_drop.fix(4e4)

    print("DOF IS", degrees_of_freedom(m))
    m.fs.unit.initialize(optarg=solver.options)
    m.fs.unit.report()

    m.fs.unit.visc_d.pprint()
    m.fs.unit.dens_mass.pprint()
    m.fs.unit.pressure_drop_total.pprint() if hasattr(
        m.fs.unit, "pressure_drop_total"
    ) else print("No it doesn't have pressure_drop total ")

    m.fs.unit.pressure_drop.pprint() if hasattr(m.fs.unit, "pressure_drop") else print(
        "No it doesn't have pressure_drop "
    )
    m.fs.unit.hydraulic_diameter.pprint() if hasattr(
        m.fs.unit, "hydraulic_diameter"
    ) else print("No it doesn't have hydraulic diameter")

    m.fs.unit.current_dens_lim_ioa.pprint() if hasattr(
        m.fs.unit, "current_dens_lim_ioa"
    ) else print("No it doesn't have current_dens_lim_ioa")
    # m.fs.unit.pressure_drop_total.fix(0)

    m.fs.unit.diffus_mass.pprint() if hasattr(m.fs.unit, "diffus_mass") else print(
        "No it doesn't have diffus_mass"
    )

    m.fs.unit.friction_factor.pprint() if hasattr(
        m.fs.unit, "friction_factor"
    ) else print("No it doesn't have friction_factor")
    m.fs.unit.N_Re.pprint() if hasattr(m.fs.unit, "N_Re") else print(
        "No it doesn't have friction_factor"
    )
    m.fs.unit.velocity_diluate.pprint()
    m.fs.unit.hydraulic_diameter.pprint() if hasattr(
        m.fs.unit, "hydraulic_diameter"
    ) else print("No it doesn't have hydraulic_diameter")

    return m


if __name__ == "__main__":
    m = main()
