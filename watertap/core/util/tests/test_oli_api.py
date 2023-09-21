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

# imports from WaterTAP repo
from watertap.core.util.credential_manager import CredentialManager
from watertap.core.util.oli_api import OLIApi

from pyomo.environ import ConcreteModel, assert_optimal_termination
from idaes.core import FlowsheetBlock
from idaes.core.util.scaling import calculate_scaling_factors
from idaes.core.solvers import get_solver
import watertap.property_models.multicomp_aq_sol_prop_pack as props

__author__ = "Adam Atia, Paul Vecchiarelli"


# Set up solver
solver = get_solver()

# TODO: alter tests; login and dependent functions now in credential_manager.py
# TODO: write requests scripts that simulate login success/failure and return dummy objects in expected format
@pytest.mark.unit
def test_login():

    # Enter login credentials below, or leave empty and submit during user input prompt
    username = "dummy@dummy.edu"
    password = "dummy_pass"
    root_url = "https://dummy_root.com"
    auth_url = "https://dummy_url.com/dummy"
    config_file = None
    encryption_key = None
    
    oliapi = OLIApi(
        username=username, password=password, root_url=root_url, auth_url=auth_url,
        config_file=config_file, encryption_key=encryption_key, credential_manager_class=CredentialManager
    )

    # TODO: Want capability to test successful login (as well as testing desired exceptions to be raise upon intentional fail)
    # is desired from OLI.
    # login will fail due to dummy credentials right now, with unintentional exception.
    # oliapi.login(fail_flag=False)


# Start test class
# TODO: Consider using dummy metadata rather than importing property package
class TestOLIAPI_WaterTAP:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = props.MCASParameterBlock(
            solute_list=[
                "NAION",
                "CLION",
                "CAION",
                "SO4ION",
                "MGION",
                "KION",
                "HCO3ION",
            ],
            mw_data={
                "H2O": 18e-3,
                "NAION": 22.989770e-3,
                "CLION": 35.45e-3,
                "CAION": 40.08e-3,
                "SO4ION": 96.06e-3,
                "MGION": 24.305e-3,
                "KION": 39.10e-3,
                "HCO3ION": 61.02e-3,
            },
            charge={
                "NAION": 1,
                "CLION": -1,
                "CAION": 2,
                "SO4ION": -2,
                "MGION": 2,
                "KION": 1,
                "HCO3ION": -1,
            },
        )
        stream = m.fs.stream = m.fs.properties.build_state_block([0])
        stream[0].temperature.fix(298)
        stream[0].pressure.fix(101325)

        stream.calculate_state(
            var_args={
                ("mass_frac_phase_comp", ("Liq", "NAION")): 10556e-6,
                ("mass_frac_phase_comp", ("Liq", "CLION")): 18980e-6,
                ("mass_frac_phase_comp", ("Liq", "MGION")): 1262e-6,
                ("mass_frac_phase_comp", ("Liq", "CAION")): 400e-6,
                ("mass_frac_phase_comp", ("Liq", "SO4ION")): 2649e-6,
                ("mass_frac_phase_comp", ("Liq", "KION")): 380e-6,
                ("mass_frac_phase_comp", ("Liq", "HCO3ION")): 140e-6,
                ("flow_vol_phase", "Liq"): 1e-3,
            },
            hold_state=True,  # fixes the calculated component mass flow rates
        )

        stream[0].conc_mass_phase_comp
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1, index=("Liq", "NAION")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1, index=("Liq", "CLION")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1, index=("Liq", "MGION")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1, index=("Liq", "CAION")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1, index=("Liq", "SO4ION")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1, index=("Liq", "KION")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1, index=("Liq", "HCO3ION")
        )

        calculate_scaling_factors(m)
        stream.initialize()
        res = solver.solve(m)
        assert_optimal_termination(res)

        return m

    @pytest.mark.unit
    def test_example(self, model):
        m = model
        # Log-in to OLI API
        username = "dummy@dummy.edu"
        password = "dummy_pass"
        root_url = "https://dummy_root.com"
        auth_url = "https://dummy_url.com/dummy"
        config_file = None
        encryption_key = None
        
        try:
            oliapi = OLIApi(
                username=username,
                password=password,
                root_url=root_url,
                auth_url=auth_url,
                config_file=config_file,
                encryption_key=encryption_key,
                credential_manager_class=CredentialManager,
            )
        except ConnectionError:
            pass
        # Continue if login  successful
        # if oliapi.login():

        # All code below would be wrapped in the conditional above to run only if oli login was successful

        # TODO: Need OLI testing for automatic dbs file generation to run line below:
        # Create chemistry file ID using OLI's chem-builder and the solute set defined in WaterTAP's MCAS
        # chemistry_file_id = oliapi.get_dbs_file_id(ions=m.fs.properties.solute_set, model_name="Test_Seawater")

        # Use MCAS stateblock to create dictionary with input concentration data required for OLI call function.
        # Will add check to see in dict is as expected
        brine_input_clone = oliapi.create_input_dict(m.fs.stream)

        # TODO: Need OLI testing to call wateranalysis function in line below:
        # call OLI's wateranalysis function and save results
        # result = oliapi.call("wateranalysis", chemistry_file_id, brine_input_clone)

        # TODO: Need OLI testing to get the job id in line below:
        # save jobid for troubleshooting with OLI support team
        # jobid = oliapi.get_job_id(chemistry_file_id)

        # TODO: Need OLI testing to get scaling tendency example results to verify a single value for gypsum scaling,
        # specific to feed comp specified earlier:
        # Extract scaling tendency for gypsum as an example and verify value is as expected within tolerance
        # SI = result['result']['additionalProperties']['scalingTendencies']['values']['CASO4.2H2O']
        # assert SI == pytest.approx(1, rel=1e-2)

        # TODO: Need OLI testing to test flash history method
        # The function below just gives us a record
        # flash_hist = oliapi.get_flash_history(chemistry_file_id)

    @pytest.mark.unit
    def test_create_dbs_dict(self, model):
        m = model
        username = "dummy@dummy.edu"
        password = "dummy_pass"
        root_url = "https://dummy_root.com"
        auth_url = "https://dummy_url.com/dummy"
        config_file = None
        encryption_key = None
        
        try:
            oliapi = OLIApi(
                username=username,
                password=password,
                root_url=root_url,
                auth_url=auth_url,
                config_file=config_file,
                encryption_key=encryption_key,
                credential_manager_class=CredentialManager,
            )
        except ConnectionError:
            pass
        brine_input_clone = oliapi.create_input_dict(m.fs.stream)
        test_clone = {
            "params": {
                "waterAnalysisInputs": [
                    {
                        "group": "Cations",
                        "name": "NAION",
                        "unit": "mg/L",
                        "value": 10555.99999998432,
                        "charge": 1,
                    },
                    {
                        "group": "Anions",
                        "name": "CLION",
                        "unit": "mg/L",
                        "value": 18979.999999969295,
                        "charge": -1,
                    },
                    {
                        "group": "Cations",
                        "name": "CAION",
                        "unit": "mg/L",
                        "value": 399.9999999402999,
                        "charge": 2,
                    },
                    {
                        "group": "Anions",
                        "name": "SO4ION",
                        "unit": "mg/L",
                        "value": 2648.9999999917354,
                        "charge": -2,
                    },
                    {
                        "group": "Cations",
                        "name": "MGION",
                        "unit": "mg/L",
                        "value": 1261.9999999828099,
                        "charge": 2,
                    },
                    {
                        "group": "Cations",
                        "name": "KION",
                        "unit": "mg/L",
                        "value": 379.9999999812004,
                        "charge": 1,
                    },
                    {
                        "group": "Anions",
                        "name": "HCO3ION",
                        "unit": "mg/L",
                        "value": 140.00000186069485,
                        "charge": -1,
                    },
                    {
                        "group": "Properties",
                        "name": "Temperature",
                        "unit": "°C",
                        "value": 24.850000000000023,
                    },
                    {
                        "group": "Properties",
                        "name": "Pressure",
                        "unit": "Pa",
                        "value": 101325,
                    },
                    {
                        "group": "Electroneutrality Options",
                        "name": "ElectroNeutralityBalanceType",
                        "value": "DominantIon",
                    },
                    {
                        "group": "Calculation Options",
                        "name": "CalcType",
                        "value": "EquilCalcOnly",
                    },
                    {
                        "group": "Calculation Options",
                        "name": "CalcAlkalnity",
                        "value": False,
                    },
                    {
                        "group": "Calculation Options",
                        "name": "AllowSolidsToForm",
                        "value": False,
                    },
                ],
                "optionalProperties": {
                    "scalingIndex": True,
                    "scalingTendencies": True,
                    "kValuesMBased": True,
                },
                "unitSetInfo": {
                    "tds": "mg/L",
                    "solid_phs_comp": "g/g",
                    "liquid_phs_comp": "mg/L",
                },
            }
        }
        assert isinstance(brine_input_clone, dict)
        assert brine_input_clone.keys() == test_clone.keys()
        assert brine_input_clone["params"].keys() == test_clone["params"].keys()
        for k in brine_input_clone["params"].keys():
            test_water_inputs = test_clone["params"][k]
            brine_inputs = brine_input_clone["params"][k]
            if k == "waterAnalysisInputs":
                assert isinstance(brine_inputs, list)
                for i, e in enumerate(brine_inputs):
                    for k, v in e.items():
                        vtest = test_water_inputs[i][k]
                        if isinstance(vtest, (float, int)) and isinstance(
                            v, (float, int)
                        ):
                            assert v == pytest.approx(vtest, rel=1e-3)
                        else:
                            assert v == vtest
            else:
                assert brine_inputs == test_water_inputs
