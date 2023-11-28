from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.unit_models.electrodialysis_0D import (
    ElectricalOperationMode,
    Electrodialysis0D,
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


def main():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    ion_dict = {
        "solute_list": ["Na_+", "Cl_-"],
        "mw_data": {"H2O": 18e-3, "Na_+": 23e-3, "Cl_-": 35.5e-3},
        "elec_mobility_data": {("Liq", "Na_+"): 5.19e-8, ("Liq", "Cl_-"): 7.92e-8},
        "charge": {"Na_+": 1, "Cl_-": -1},
        "diffusivity_data": {("Liq", "Na_+"): 1.33e-9, ("Liq", "Cl_-"): 2.03e-9},
    }
    m.fs.properties = MCASParameterBlock(**ion_dict)
    m.fs.unit = Electrodialysis0D(
        property_package=m.fs.properties,
        operation_mode=ElectricalOperationMode.Constant_Voltage,
        has_Nernst_diffusion_layer=True,
        has_nonohmic_potential_membrane=True,
        # limiting_current_density_method = LimitingCurrentDensityMethod.Empirical,
        limiting_current_density_method=LimitingCurrentDensityMethod.Theoretical,
        # limiting_current_density_data=500,
    )

    m.fs.unit.water_trans_number_membrane["cem"].fix(5.8)
    m.fs.unit.water_trans_number_membrane["aem"].fix(4.3)
    m.fs.unit.water_permeability_membrane["cem"].fix(2.16e-14)
    m.fs.unit.water_permeability_membrane["aem"].fix(1.75e-14)
    m.fs.unit.voltage.fix(0.5)
    m.fs.unit.electrodes_resistance.fix(0)
    m.fs.unit.cell_pair_num.fix(10)
    m.fs.unit.current_utilization.fix(1)
    m.fs.unit.channel_height.fix(5e-4)
    m.fs.unit.membrane_areal_resistance["cem"].fix(1.89e-4)
    m.fs.unit.membrane_areal_resistance["aem"].fix(1.77e-4)
    m.fs.unit.cell_width.fix(0.1)
    m.fs.unit.cell_length.fix(0.79)
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

    # check ion transfer number requirements
    assert (
        sum(
            value(m.fs.unit.ion_trans_number_membrane["cem", j])
            for j in m.fs.properties.ion_set
        )
        == 1
    )
    assert (
        sum(
            value(m.fs.unit.ion_trans_number_membrane["aem", j])
            for j in m.fs.properties.ion_set
        )
        == 1
    )
    assert sum(
        value(m.fs.unit.ion_trans_number_membrane["cem", j])
        for j in m.fs.properties.cation_set
    ) == sum(
        value(m.fs.unit.ion_trans_number_membrane["aem", j])
        for j in m.fs.properties.anion_set
    )
    # set the inlet stream
    m.fs.unit.inlet_diluate.pressure.fix(101325)
    m.fs.unit.inlet_diluate.temperature.fix(298.15)
    m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
    m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-4)
    m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-4)
    m.fs.unit.inlet_concentrate.pressure.fix(101325)
    m.fs.unit.inlet_concentrate.temperature.fix(298.15)
    m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
    m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-4)
    m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-4)
    m.fs.unit.spacer_porosity.fix(1)

    # assert degrees_of_freedom(m) == 0
    # set default scaling for state vars
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e3, index=("Liq", "Na_+")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e3, index=("Liq", "Cl_-")
    )

    iscale.calculate_scaling_factors(m.fs)

    m.fs.unit.diffus_mass.fix(1.6e-9) if hasattr(m.fs.unit, "diffus_mass") else 0
    m.fs.unit.spacer_specific_area.fix(10700) if hasattr(
        m.fs.unit, "spacer_specific_area"
    ) else 0
    # m.fs.unit.hydraulic_diameter.fix(1.5e-3) if hasattr(
    #     m.fs.unit, "hydraulic_diameter"
    # ) else 0
    # m.fs.unit.friction_factor.fix(20) if hasattr(
    #     m.fs.unit, "friction_factor"
    # ) else 0

    print(m.fs.unit.config.limiting_current_density_data)
    print("DOF IS", degrees_of_freedom(m))
    # m.display()
    solver = get_solver()
    # initialization_tester(m, outlvl=idaeslog.DEBUG)

    m.fs.unit.initialize(optarg=solver.options)
    m.fs.unit.report()
    m.fs.unit.current_dens_lim_ioa.pprint()

    results = solver.solve(m)
    print(results.solver.termination_condition)
    return m


if __name__ == "__main__":
    m = main()
