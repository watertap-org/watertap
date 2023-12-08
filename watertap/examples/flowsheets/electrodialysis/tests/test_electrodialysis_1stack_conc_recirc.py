#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
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
from pyomo.environ import value
from idaes.core.solvers import get_solver
from watertap.core.util.initialization import check_dof
import watertap.examples.flowsheets.electrodialysis.electrodialysis_1stack_conc_recirc as edfs

__author__ = "Xiangyu Bi"


class TestElectrodialysis1StackFS:
    @pytest.fixture(scope="class")
    def ED1D1Stack_conc_recirc(self):
        m = edfs.build()
        return m

    @pytest.mark.component
    def test_specific_operating_conditions(self, ED1D1Stack_conc_recirc):
        m = ED1D1Stack_conc_recirc
        solver = get_solver()
        # Testing a feeding salinity of 2g/L.
        init_arg = {
            ("flow_vol_phase", ("Liq")): 5.2e-4,
            ("conc_mol_phase_comp", ("Liq", "Na_+")): 34.188,
            ("conc_mol_phase_comp", ("Liq", "Cl_-")): 34.188,
        }
        m.fs.feed.properties.calculate_state(
            init_arg,
            hold_state=True,
        )
        m.fs.EDstack.voltage_applied[0].fix(10)
        m.fs.recovery_vol_H2O.fix(0.7)
        edfs._condition_base(m)
        check_dof(m)

        # Initialize and solve the model
        edfs.initialize_system(m, solver=solver)
        edfs.solve(m, solver=solver)

        assert value(m.fs.feed.properties[0].flow_vol_phase["Liq"]) == pytest.approx(
            5.2e-4, rel=1e-3
        )
        assert value(m.fs.prod.properties[0].flow_vol_phase["Liq"]) == pytest.approx(
            3.64e-4, rel=1e-3
        )
        assert value(m.fs.disp.properties[0].flow_vol_phase["Liq"]) == pytest.approx(
            1.56e-4, rel=1e-3
        )
        assert value(
            sum(
                m.fs.feed.properties[0].conc_mass_phase_comp["Liq", j]
                for j in m.fs.properties.ion_set
            )
        ) == pytest.approx(2.00, rel=1e-3)
        assert value(m.fs.product_salinity) == pytest.approx(1.1333, rel=1e-3)
        assert value(m.fs.disposal_salinity) == pytest.approx(4.0223, rel=1e-3)
        assert value(m.fs.mem_area) == pytest.approx(18.5338, rel=1e-3)
        assert value(m.fs.costing.specific_energy_consumption) == pytest.approx(
            0.1348, abs=0.001
        )
        assert value(m.fs.costing.LCOW) == pytest.approx(0.3996, rel=1e-3)
        assert value(m.fs.EDstack.inlet_concentrate.pressure[0]) == pytest.approx(
            169278.127, rel=1e-3
        )
        assert value(m.fs.EDstack.outlet_concentrate.pressure[0]) == pytest.approx(
            101325.00, rel=1e-3
        )
        assert value(m.fs.EDstack.inlet_diluate.pressure[0]) == pytest.approx(
            169278.127, rel=1e-3
        )
        assert value(m.fs.EDstack.outlet_diluate.pressure[0]) == pytest.approx(
            101325.00, rel=1e-3
        )
        assert value(m.fs.EDstack.inlet_concentrate.temperature[0]) == pytest.approx(
            298.15, rel=1e-3
        )
        assert value(m.fs.EDstack.outlet_concentrate.temperature[0]) == pytest.approx(
            298.15, rel=1e-3
        )
        assert value(m.fs.EDstack.inlet_diluate.temperature[0]) == pytest.approx(
            298.15, rel=1e-3
        )
        assert value(m.fs.EDstack.outlet_diluate.temperature[0]) == pytest.approx(
            298.15, rel=1e-3
        )

    @pytest.mark.component
    def test_optimization(self):
        m = edfs.main()
        assert value(m.fs.feed.properties[0].flow_vol_phase["Liq"]) == pytest.approx(
            5.2e-4, rel=1e-3
        )
        assert value(m.fs.prod.properties[0].flow_vol_phase["Liq"]) == pytest.approx(
            3.64e-4, rel=1e-3
        )
        assert value(m.fs.disp.properties[0].flow_vol_phase["Liq"]) == pytest.approx(
            1.56e-4, rel=1e-3
        )
        assert value(
            sum(
                m.fs.feed.properties[0].conc_mass_phase_comp["Liq", j]
                for j in m.fs.properties.ion_set
            )
        ) == pytest.approx(2.00, rel=1e-3)
        assert value(m.fs.product_salinity) == pytest.approx(0.1000, rel=1e-3)
        assert value(m.fs.disposal_salinity) == pytest.approx(6.4333, rel=1e-3)
        assert value(m.fs.mem_area) == pytest.approx(15.1283, rel=1e-3)
        assert value(m.fs.EDstack.cell_pair_num) == pytest.approx(16, rel=1e-8)
        assert value(m.fs.EDstack.cell_length) == pytest.approx(4.800, rel=1e-3)
        assert value(m.fs.EDstack.voltage_applied[0]) == pytest.approx(18.785, rel=1e-3)
        assert value(m.fs.costing.specific_energy_consumption) == pytest.approx(
            1.7775, rel=1e-3
        )
        assert value(m.fs.costing.LCOW) == pytest.approx(0.6254, rel=1e-3)
        assert value(m.fs.EDstack.inlet_concentrate.pressure[0]) == pytest.approx(
            784785.332, rel=1e-3
        )
        assert value(m.fs.EDstack.outlet_concentrate.pressure[0]) == pytest.approx(
            101325.00, rel=1e-3
        )
        assert value(m.fs.EDstack.inlet_diluate.pressure[0]) == pytest.approx(
            784785.332, rel=1e-3
        )
        assert value(m.fs.EDstack.outlet_diluate.pressure[0]) == pytest.approx(
            101325.00, rel=1e-3
        )
        assert value(m.fs.EDstack.inlet_concentrate.temperature[0]) == pytest.approx(
            298.15, rel=1e-3
        )
        assert value(m.fs.EDstack.outlet_concentrate.temperature[0]) == pytest.approx(
            298.15, rel=1e-3
        )
        assert value(m.fs.EDstack.inlet_diluate.temperature[0]) == pytest.approx(
            298.15, rel=1e-3
        )
        assert value(m.fs.EDstack.outlet_diluate.temperature[0]) == pytest.approx(
            298.15, rel=1e-3
        )
