"""
Test for api module, with real workflow
"""
import pprint
import pytest
from watertap.ui.api import AnalysisWorkflow, Build, Initialize, Optimize, Steps
from watertap.examples.flowsheets.case_studies.seawater_RO_desalination import (
    seawater_RO_desalination as srd,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.environ import units as pyunits
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


class BuildRO(Build):
    def build_model(self, input_data):
        kwargs = self.flowsheet_data
        model = srd.build(**kwargs)
        set_operating_conditions_p(model, input_data)
        srd.assert_degrees_of_freedom(model, 0)
        return model


class InitializeRO(Initialize):
    def initialize_model(self, data):
        model = self.workflow.model
        srd.initialize_system(model)
        srd.assert_degrees_of_freedom(model, 0)


class Solve(Optimize):
    def solve(self, data):
        return srd.solve(self.workflow.model)


class AddCosting(Build):
    def build_model(self, data):
        # assume model already built
        model = self.workflow.model
        srd.add_costing(model)
        return model

class InitCosting(Initialize):
    def initialize_model(self, data):
        model = self.workflow.model
        srd.initialize_costing(model)
        srd.assert_degrees_of_freedom(model, 0)


perf_vars = dict(
    feed=dict(
        flow_vol={
            "value": 0.3092,
            "units": pyunits.m**3 / pyunits.s,
            "desc": "flow volume",
        },
        conc_mass_tds={
            "value": 35,
            "units": pyunits.kg / pyunits.m**3,
            "desc": "concentration of total dissolved solids",
        },
        conc_mass_tss={
            "value": 0.03,
            "units": pyunits.kg / pyunits.m**3,
            "desc": "concentration of total SS",
        },
        temperature={"value": 298, "units": pyunits.K, "desc": "temperature"},
        pressure={"value": 1e5, "units": pyunits.Pa, "desc": "pressure"},
    )
)

cost_vars = dict()


@pytest.mark.component
def test_workflow():
    wf = AnalysisWorkflow()
    wf.set_flowsheet_data({"erd_type": "pressure_exchanger"})
    # TODO: set strategy and data (if any) for each step
    wf.set_strategy(Steps.build, BuildRO)
    wf.set_input(Steps.build, perf_vars)
    wf.set_strategy(Steps.init, InitializeRO)
    wf.set_strategy(Steps.optimize, Solve)
    wf.set_strategy(Steps.build_costing, AddCosting)
    wf.set_strategy(Steps.init_costing, InitCosting)
    wf.set_strategy(Steps.optimize_costing, Solve)   # re-use solve class
    wf.set_standard_workflow(costing=True)
    print("Workflow is created")
    print("Performance variables:")
    pprint.pprint(wf.get_input(Steps.build))
    print("Costing variables:")
    pprint.pprint(wf.get_input(Steps.build_costing))
    wf.run_all()
    # TODO: Check results
    print(wf.optimize_result)


## Some modified methods


def uv(entry):
    return entry["value"] * entry["units"]


def set_operating_conditions_p(m, p):
    _log.debug(f"set_operating_conditions_p: parameters={p}")

    prtrt = m.fs.pretreatment
    desal = m.fs.desalination
    psttrt = m.fs.posttreatment

    # ---specifications---
    # feed
    # flow_vol = 0.3092 * pyunits.m**3 / pyunits.s
    # conc_mass_tds = 35 * pyunits.kg / pyunits.m**3
    # conc_mass_tss = 0.03 * pyunits.kg / pyunits.m**3
    # temperature = 298 * pyunits.K
    # pressure = 1e5 * pyunits.Pa
    #

    # Set feed parameters
    feed = p["feed"]
    flow_vol = uv(feed["flow_vol"])
    conc_mass_tds = uv(feed["conc_mass_tds"])
    conc_mass_tss = uv(feed["conc_mass_tss"])
    temperature = uv(feed["temperature"])
    pressure = uv(feed["pressure"])

    m.fs.feed.flow_vol[0].fix(flow_vol)
    m.fs.feed.conc_mass_comp[0, "tds"].fix(conc_mass_tds)
    m.fs.feed.conc_mass_comp[0, "tss"].fix(conc_mass_tss)
    srd.solve(m.fs.feed)

    m.fs.tb_prtrt_desal.properties_out[0].temperature.fix(temperature)
    m.fs.tb_prtrt_desal.properties_out[0].pressure.fix(pressure)

    # ---pretreatment---
    # intake

    # ferric chloride
    m.db.get_unit_operation_parameters("chemical_addition")
    prtrt.ferric_chloride_addition.load_parameters_from_database()
    prtrt.ferric_chloride_addition.chemical_dosage.fix(20)

    # chlorination
    m.db.get_unit_operation_parameters("chlorination")
    prtrt.chlorination.load_parameters_from_database(use_default_removal=True)

    # static mixer
    m.db.get_unit_operation_parameters("static_mixer")
    prtrt.static_mixer.load_parameters_from_database(use_default_removal=True)

    # storage tank
    m.db.get_unit_operation_parameters("storage_tank")
    prtrt.storage_tank_1.load_parameters_from_database(use_default_removal=True)
    prtrt.storage_tank_1.storage_time.fix(2)

    # media filtration
    m.db.get_unit_operation_parameters("media_filtration")
    prtrt.media_filtration.load_parameters_from_database(use_default_removal=True)

    # backwash handling
    m.db.get_unit_operation_parameters("backwash_solids_handling")
    prtrt.backwash_handling.load_parameters_from_database(use_default_removal=True)

    # anti-scalant
    prtrt.anti_scalant_addition.load_parameters_from_database()
    prtrt.anti_scalant_addition.chemical_dosage.fix(5)

    # cartridge filtration
    m.db.get_unit_operation_parameters("cartridge_filtration")
    prtrt.cartridge_filtration.load_parameters_from_database(use_default_removal=True)

    # ---desalination---
    # pump 1, high pressure pump, 2 degrees of freedom (efficiency and outlet pressure)
    desal.P1.efficiency_pump.fix(0.80)  # pump efficiency [-]
    operating_pressure = 70e5 * pyunits.Pa
    desal.P1.control_volume.properties_out[0].pressure.fix(operating_pressure)

    # RO unit
    desal.RO.A_comp.fix(4.2e-12)  # membrane water permeability coefficient [m/s-Pa]
    desal.RO.B_comp.fix(3.5e-8)  # membrane salt permeability coefficient [m/s]
    desal.RO.channel_height.fix(1e-3)  # channel height in membrane stage [m]
    desal.RO.spacer_porosity.fix(0.97)  # spacer porosity in membrane stage [-]
    desal.RO.permeate.pressure[0].fix(101325)  # atmospheric pressure [Pa]
    desal.RO.width.fix(1000)  # stage width [m]
    desal.RO.area.fix(
        flow_vol * 4.5e4 * pyunits.s / pyunits.m
    )  # stage area [m2] TODO: replace with actual area

    if m.erd_type == "pressure_exchanger":
        # splitter (no degrees of freedom)

        # pressure exchanger, 1 degree of freedom (efficiency)
        desal.PXR.efficiency_pressure_exchanger.fix(0.95)

        # pump 2, booster pump, 1 degree of freedom (efficiency, pressure must match high pressure pump)
        desal.P2.efficiency_pump.fix(0.80)

        # mixer, no degrees of freedom
    elif m.erd_type == "pump_as_turbine":
        # ERD, 2 degrees of freedom (efficiency, outlet pressure)
        desal.ERD.efficiency_pump.fix(0.95)
        desal.ERD.control_volume.properties_out[0].pressure.fix(
            101325
        )  # atmospheric pressure [Pa]

    # ---posttreatment---
    # storage tank 2
    psttrt.storage_tank_2.load_parameters_from_database(use_default_removal=True)
    psttrt.storage_tank_2.storage_time.fix(1)

    # uv aop
    m.db.get_unit_operation_parameters("uv_aop")
    psttrt.uv_aop.load_parameters_from_database(use_default_removal=True)
    psttrt.uv_aop.uv_reduced_equivalent_dose.fix(
        350
    )  # TODO: check this was the right thing to fix
    psttrt.uv_aop.uv_transmittance_in.fix(
        0.95
    )  # TODO: check this was the right thing to fix

    # co2 addition
    m.db.get_unit_operation_parameters("co2_addition")
    psttrt.co2_addition.load_parameters_from_database(use_default_removal=True)

    # lime
    psttrt.lime_addition.load_parameters_from_database()
    psttrt.lime_addition.chemical_dosage.fix(2.3)

    # storage tank 3
    psttrt.storage_tank_3.load_parameters_from_database(use_default_removal=True)
    psttrt.storage_tank_3.storage_time.fix(1)

    # ---product and disposal---
    m.db.get_unit_operation_parameters("municipal_drinking")
    m.fs.municipal.load_parameters_from_database()

    m.db.get_unit_operation_parameters("landfill")
    m.fs.landfill.load_parameters_from_database()
