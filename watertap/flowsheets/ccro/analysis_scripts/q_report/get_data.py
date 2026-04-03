import os
from psPlotKit.data_manager.ps_data_manager import PsDataManager
from psPlotKit.data_manager.ps_costing import (
    PsCostingGroup,
    PsCostingManager,
)
from psPlotKit.data_manager.costing_packages.watertap_costing import (
    WaterTapCostingPackage,
)

here = os.path.dirname(os.path.abspath(__file__))
par_dir = os.path.dirname(here)


def get_ss_data():

    dm = PsDataManager()
    dm.register_data_file(
        f"{par_dir}/output/RO_w_ERD_analysisType_BW_recovery_sweep.h5",
        directory="Brackish water",
    )
    dm.register_data_file(
        f"{par_dir}/output/RO_w_ERD_analysisType_PW_recovery_sweep.h5",
        directory="Produced water",
    )
    dm.register_data_file(
        f"{par_dir}/output/RO_w_ERD_analysisType_SW_recovery_sweep.h5",
        directory="Seawater",
    )
    for stage in [1, 2]:
        dm.register_data_key(
            f"fs.stage[{stage}].RO.flux_mass_phase_comp_avg[0.0,Liq,H2O]",
            (f"RO {stage}", "Flux"),
            assign_units="L/(m^2*hr)",
            conversion_factor=3600,
        )
        dm.register_data_key(
            f"fs.stage[{stage}].pump.control_volume.properties_out[0.0].pressure",
            (f"RO {stage}", "Pressure"),
            "bar",
        )
        dm.register_data_key(
            f"fs.stage[{stage}].RO.area",
            (f"RO {stage}", "Area"),
        )
        # dm.register_data_key(
        #     f"fs.stage[{stage}].pump.control_volume.work[0.0]",
        #     (f"RO {stage}", "pump work"),
        #     "kW",
        # )
        dm.register_data_key(
            f"fs.stage[{stage}].RO.feed_side.properties[0.0,0.0].conc_mass_phase_comp[Liq,NaCl]",
            (f"RO {stage}", "Inlet concentration"),
            "g/L",
        )
        dm.register_data_key(
            f"fs.stage[{stage}].pump.control_volume.work[0.0]",
            (f"RO {stage}", "Pump work"),
            "kW",
        )
        # dm.register_data_key(
        #     f"fs.stage[{stage}].pump.work_fluid[0.0]",
        #     (f"RO {stage}", "Pump size"),
        #     "kW",
        # )
    # dm.register_data_key("fs.ERD.control_volume.work[0.0]", "ERD work", "kW")
    # dm.register_data_key("fs.system_recovery", "Water recovery", "%")

    dm.register_data_key("fs.system_recovery", "Water recovery", "%")
    dm.register_data_key("fs.costing.LCOW", "LCOW")
    dm.register_data_key("fs.costing.SEC", "SEC")
    dm.register_data_key(
        "fs.costing.total_capital_cost",
        "CAPEX",
        assign_units="USD",
        conversion_factor=1e-3,
    )
    dm.register_data_key(
        "fs.costing.total_operating_cost",
        "OPEX",
        assign_units="USD/yr",
        conversion_factor=1e-3,
    )

    dm.register_data_key(
        "fs.product.properties[0.0].flow_vol_phase[Liq]",
        "Permeate flow rate",
        "L/hr",
    )
    dm.register_data_key(
        "fs.product.properties[0.0].conc_mass_phase_comp[Liq,NaCl]",
        "Permeate concentration",
        "g/L",
    )
    dm.register_data_key(
        "fs.disposal.properties[0.0].conc_mass_phase_comp[Liq,NaCl]",
        "Brine concentration",
        "g/L",
    )
    dm.register_data_key(
        "fs.total_area",
        "Membrane Area",
        assign_units="m^2",
    )
    dm.register_data_key(
        "fs.specific_area",
        "Specific area",
        assign_units="m^2/m^3/s",
    )
    dm.register_data_key(
        "fs.system_flux",
        "Flux",
        assign_units="L/(m^2*hr)",
    )
    dm.register_data_key(
        "fs.total_power",
        "Pump work",
        assign_units="kW",
    )

    dm.load_data()
    ek = dm.get_expression_keys()
    dm.register_expression(
        ek.CAPEX / ek.Permeate_flow_rate,
        "Specific CAPEX",
        assign_units="USD/(m^3/s)",
    )
    dm.register_expression(
        ek.OPEX / ek.Permeate_flow_rate,
        "Specific OPEX",
        assign_units="USD/(m^3/s)",
    )
    # for stage in [1, 2]:
    #     dm.register_expression(
    #         ek[f"RO_{stage}_Pump_work"] / ek.Permeate_flow_rate,
    #         (f"RO {stage}", "Specific pump work"),
    #         assign_units="kW/(m^3/s)",
    #         # conversion_factor=3600,
    #     )
    # dm.register_expression(
    #     ek.RO_2_Pump_size / ek.Permeate_flow_rate,
    #     (f"RO {stage}", "Specific pump size"),
    #     assign_units="kW/(m^3/s)",
    #     # conversion_factor=3600,
    # )
    # dm.register_expression(
    #     (ek.RO_1_Pump_work + ek.RO_2_Pump_work),
    #     "Pump work",
    #     assign_units="kW",
    # )
    # dm.register_expression(
    #     (ek.RO_1_Pump_work + ek.RO_2_Pump_work) / ek.Permeate_flow_rate,
    #     "Specific pump work",
    #     assign_units="kW/(m^3/s)",
    # )
    dm.register_expression(
        (ek.RO_1_Pump_work + ek.RO_2_Pump_work) / ek.Permeate_flow_rate,
        "Specific pump work",
        assign_units="kW/(m^3/s)",
    )
    # dm.register_expression(
    #     (ek.RO_1_Area + ek.RO_2_Area) / ek.Permeate_flow_rate,
    #     "Specific area",
    #     assign_units="m^2/L/hr",
    # )
    # dm.register_expression(
    #     ek.RO_1_pump_work / ek.Permeate_flow_rate,
    #     "RO 1 Specific pump work",
    #     "kW/(m^3/s)",
    # )
    # dm.register_expression(
    #     ek.RO_2_pump_work / ek.Permeate_flow_rate,
    #     "RO 2 Specific pump work",
    #     "kW/(m^3/s)",
    # )
    dm.register_expression(ek.OPEX / ek.CAPEX, "OPEX/CAPEX Ratio")
    dm.evaluate_expressions()
    dm.display()

    # dm.reduce_data(
    #     stack_keys=("Produced water", "stage_sim_cases"),
    #     data_key="LCOW",
    #     reduction_type="min",
    #     directory=("Produced water", "optimal_design"),
    # )

    # dm.display()
    # dm.reduce_data(
    #     stack_keys=("Brackish water", "stage_sim_cases"),
    #     data_key="LCOW",
    #     reduction_type="min",
    #     directory=("Brackish water", "optimal_design"),
    # )
    # dm.reduce_data(
    #     stack_keys="stage_sim_cases",
    #     data_key="LCOW",
    #     reduction_type="min",
    #     directory="optimal_design",
    # )

    # dm.display()

    package_ss = WaterTapCostingPackage(
        # costing_block="fs.costing", validation_key="fs.costing.LCOW"
    )
    # package_ss.add_flow_cost("electricity", "electricity_cost")
    package_ss.register_product_flow()

    RO1 = PsCostingGroup("Stage 1")
    RO1.add_unit(
        "stage[1].RO",
        capex_keys="capital_cost",
        fixed_opex_keys="fixed_operating_cost",
    )
    RO2 = PsCostingGroup("Stage 2")
    RO2.add_unit(
        "stage[2].RO",
        capex_keys="capital_cost",
        fixed_opex_keys="fixed_operating_cost",
    )
    pumps = PsCostingGroup("Pumps + ERD")
    pumps.add_unit(
        "pump",
        capex_keys="capital_cost",
        flow_keys={"electricity": "control_volume.work"},
    )

    pumps.add_unit(
        "ERD",
        capex_keys="capital_cost",
        flow_keys={"electricity": "control_volume.work"},
    )
    cm_ss = PsCostingManager(
        dm,
        package_ss,
        [RO1, RO2, pumps],
    )
    cm_ss.build()
    dm.register_expression(
        ek[("costing", "total", "LCOW_opex")] / ek[("costing", "total", "LCOW")],
        "OPEX LCOW Fraction",
    )
    dm.register_expression(
        ek[("costing", "total", "LCOW_capex")] / ek[("costing", "total", "LCOW")],
        "CAPEX LCOW Fraction",
    )
    dm.register_expression(
        1 / ek.Flux,
        "Specific Membrane Area",
    )

    return dm


def get_ccro_recov_data(last_block=19, n_flushing_points=5):

    last_filt_blk = int(last_block - n_flushing_points)

    dm = PsDataManager()

    dm.register_data_file(
        f"{par_dir}/output/ccro_analysisType_BW_recovery_sweep.h5",
        directory="Brackish water",
    )
    dm.register_data_file(
        f"{par_dir}/output/ccro_analysisType_SW_recovery_sweep.h5",
        directory="Seawater",
    )
    dm.register_data_file(
        f"{par_dir}/output/ccro_analysisType_PW_recovery_sweep.h5",
        directory="Produced water",
    )

    # COSTING
    dm.register_data_key("overall_recovery", "Water recovery", "%")
    dm.register_data_key("costing.LCOW", "LCOW", assign_units="USD/m^3")
    dm.register_data_key("costing.SEC", "SEC", assign_units="kWh/m^3")

    dm.register_data_key(
        "costing.total_capital_cost",
        "CAPEX",
        assign_units="USD",
        conversion_factor=1e-3,
    )
    dm.register_data_key(
        "costing.total_operating_cost",
        "OPEX",
        assign_units="USD/year",
        conversion_factor=1e-3,
    )
    # SYSTEM

    dm.register_data_key("flushing.flushing_efficiency", "Flushing efficiency", "%")
    dm.register_data_key(
        "recycle_loop_concentration", "Recycle loop concentration", "g/L"
    )
    dm.register_data_key(
        "flushing.pre_flushing_concentration", "Pre Flushing Concentration", "g/L"
    )
    # dm.register_data_key("recycle_flowrate", "Recycle rate", "L/s")
    dm.register_data_key(
        "blocks[0].process.fs.P2.control_volume.properties_out[0.0].flow_vol_phase[Liq]",
        "Recycle rate",
        "L/s",
    )
    dm.register_data_key(
        "blocks[0].process.fs.raw_feed.properties[0.0].pressure_osm_phase[Liq]",
        "Feed osmotic pressure",
        "bar",
    )
    dm.register_data_key(
        "blocks[0].process.fs.raw_feed.properties[0.0].flow_vol_phase[Liq]",
        "Feed flow rate",
        "L/s",
    )
    dm.register_data_key("avg_feed_flow_rate", "Average feed flow rate")
    dm.register_data_key("avg_product_flow_rate", "Average product flow rate", "m^3/s")
    # dm.register_data_key("avg_product_flow_rate", "Permeate flow rate", "m^3/s")
    dm.register_data_key("filtration_ramp_rate", "Ramp rate", "bar/min")
    dm.register_data_key("total_cycle_time", "Total cycle time", "min")
    dm.register_data_key("overall_rejection", "Overall rejection", "%")
    dm.register_data_key("cycle_time_ratio", "Cycle time ratio", "%")
    dm.register_data_key("permeate_concentration", "Permeate concentration", "g/L")
    dm.register_data_key("total_permeate_vol", "Total permeate volume", "m^3")
    # dm.register_data_key(
    #     "blocks[0].process.fs.P2.control_volume.properties_out[0.0].pressure",
    #     "Recycle Pump Pressure",
    #     "bar",
    # )

    # RO
    dm.register_data_key("blocks[0].process.fs.RO.area", "Area", "m^2")
    dm.register_data_key("blocks[0].process.fs.RO.area", "Membrane Area", "m^2")
    dm.register_data_key("blocks[0].process.fs.RO.area", "Total Area", "m^2")
    dm.register_data_key("blocks[0].process.fs.RO.length", "Membrane Length", "m")
    dm.register_data_key("blocks[0].process.fs.RO.width", "Membrane Width", "m")
    dm.register_data_key("blocks[0].process.fs.RO.feed_side.velocity[0.0,0.0]", "Inlet Velocity", "m/s")
    dm.register_data_key(
        "blocks[0].process.fs.RO.recovery_vol_phase[0.0,Liq]",
        "Single Pass Recovery",
        "%",
    )

    dm.register_data_key(
        f"blocks[0].process.fs.RO.feed_side.properties[0.0,0.0].conc_mass_phase_comp[Liq,NaCl]",
        "Inlet concentration",
        "g/L",
    )
    dm.register_data_key(
        f"blocks[0].process.fs.RO.feed_side.properties[0.0,0.0].conc_mass_phase_comp[Liq,NaCl]",
        "Initial concentration",
        "g/L",
    )
    dm.register_data_key(
        f"blocks[0].process.fs.P1.control_volume.properties_out[0.0].pressure",
        "Initial Pressure",
        "bar",
    )
    dm.register_data_key(
        f"blocks[{last_block}].process.fs.P1.control_volume.properties_out[0.0].pressure",
        "Final Pressure",
        "bar",
    )
    dm.register_data_key(
        f"blocks[{last_block}].process.fs.P1.control_volume.properties_out[0.0].pressure",
        "Pressure",
        "bar",
    )
    # dm.register_data_key(
    #     f"blocks[{last_block}].process.fs.P1.control_volume.pressure",
    #     "Pressure",
    #     "bar",
    # )
    dm.register_data_key(
        f"blocks[{last_block}].process.fs.P1.control_volume.work[0.0]",
        "Pump size",
        "kW",
    )
    dm.register_data_key(
        f"blocks[{last_block}].process.fs.P1.total_power",
        "Pump work",
        "kW",
    )
    dm.register_data_key(
        f"blocks[{last_block}].process.fs.RO.feed_side.properties[0.0,0.0].conc_mass_phase_comp[Liq,NaCl]",
        "Final concentration",
        "g/L",
    )

    dm.load_data()

    ek = dm.get_expression_keys()
    dm.register_expression(
        ek.Pump_size / ek.Average_product_flow_rate, "Specific pump size", "kW/(m^3/s)"
    )
    dm.register_expression(
        ek.Pump_work / ek.Average_product_flow_rate, "Specific pump work", "kW/(m^3/s)"
    )
    dm.register_expression(
        ek.Average_product_flow_rate / ek.Area,
        "Flux",
        "L/(m^2*hr)",
    )
    dm.register_expression(
        ek.Area / ek.Average_product_flow_rate,
        "Specific area",
        "m^2/(L/hr)",
    )
    dm.register_expression(ek.OPEX / ek.CAPEX, "OPEX/CAPEX Ratio")
    dm.register_expression(
        ek.OPEX / ek.Average_product_flow_rate, "Specific OPEX", "USD/m^3"
    )
    dm.register_expression(
        ek.CAPEX / ek.Average_product_flow_rate, "Specific CAPEX", "USD/(m^3/s)"
    )

    dm.evaluate_expressions()
    ek = dm.get_expression_keys()
    dm.register_expression(
        ek.Recycle_rate / ek.Feed_flow_rate,
        "Recycle ratio",
    )

    # lets create costing pacakage
    package_ccro = WaterTapCostingPackage(
        costing_block="costing", validation_key="costing.LCOW"
    )
    package_ccro.register_product_flow("avg_product_flow_rate")

    # Lets create our groups
    RO = PsCostingGroup("RO")
    RO.add_unit(
        f"blocks[0].process.fs.RO",  # only adding "block[0] to specify wher capex is, normally acn just say "RO
        capex_keys="capital_cost",
        fixed_opex_keys="fixed_operating_cost",
    )
    feed_pump = PsCostingGroup("Feed pump")
    feed_pump.add_unit(
        f"blocks[{last_block}].process.fs.P1",
        capex_keys="capital_cost",
        flow_keys={"electricity": "total_power"},
    )
    recycle_pump = PsCostingGroup("Recycle pump")
    recycle_pump.add_unit(
        f"blocks[{last_block}].process.fs.P2",
        capex_keys="capital_cost",
        flow_keys={"electricity": "total_power"},
    )
    conduit = PsCostingGroup("Conduit")
    conduit.add_unit(
        "conduit",
        capex_keys="capital_cost",
    )

    cm_recov = PsCostingManager(
        dm, package_ccro, [RO, feed_pump, recycle_pump, conduit]
    )
    cm_recov.build()
    dm.register_expression(
        ek[("costing", "total", "LCOW_opex")] / ek[("costing", "total", "LCOW")],
        "OPEX LCOW Fraction",
    )
    dm.register_expression(
        ek[("costing", "total", "LCOW_capex")] / ek[("costing", "total", "LCOW")],
        "CAPEX LCOW Fraction",
    )
    dm.register_expression(
        1 / ek.Flux,
        "Specific Membrane Area",
    )
    return dm


def get_ccro_SEC_data(last_block=19, n_flushing_points=5):

    last_filt_blk = int(last_block - n_flushing_points)

    dm = PsDataManager()

    dm.register_data_file(
        f"{par_dir}/output/ccro_analysisType_BW_recovery_sweep_SEC_obj.h5",
        directory="Brackish water",
    )
    dm.register_data_file(
        f"{par_dir}/output/ccro_analysisType_SW_recovery_sweep_SEC_obj.h5",
        directory="Seawater",
    )
    dm.register_data_file(
        f"{par_dir}/output/ccro_analysisType_PW_recovery_sweep_SEC_obj.h5",
        directory="Produced water",
    )

    # COSTING
    dm.register_data_key("overall_recovery", "Water recovery", "%")
    dm.register_data_key("costing.LCOW", "LCOW", assign_units="USD/m^3")
    dm.register_data_key("costing.SEC", "SEC", assign_units="kWh/m^3")

    dm.register_data_key(
        "costing.total_capital_cost",
        "CAPEX",
        assign_units="USD",
        conversion_factor=1e-3,
    )
    dm.register_data_key(
        "costing.total_operating_cost",
        "OPEX",
        assign_units="USD/year",
        conversion_factor=1e-3,
    )
    # SYSTEM

    dm.register_data_key("flushing.flushing_efficiency", "Flushing efficiency", "%")
    dm.register_data_key(
        "recycle_loop_concentration", "Recycle loop concentration", "g/L"
    )
    dm.register_data_key(
        "flushing.pre_flushing_concentration", "Pre Flushing Concentration", "g/L"
    )
    # dm.register_data_key("recycle_flowrate", "Recycle rate", "L/s")
    dm.register_data_key(
        "blocks[0].process.fs.P2.control_volume.properties_out[0.0].flow_vol_phase[Liq]",
        "Recycle rate",
        "L/s",
    )
    dm.register_data_key(
        "blocks[0].process.fs.raw_feed.properties[0.0].pressure_osm_phase[Liq]",
        "Feed osmotic pressure",
        "bar",
    )
    dm.register_data_key(
        "blocks[0].process.fs.raw_feed.properties[0.0].flow_vol_phase[Liq]",
        "Feed flow rate",
        "L/s",
    )
    dm.register_data_key("avg_feed_flow_rate", "Average feed flow rate")
    dm.register_data_key("avg_product_flow_rate", "Average product flow rate", "m^3/s")
    # dm.register_data_key("avg_product_flow_rate", "Permeate flow rate", "m^3/s")
    dm.register_data_key("filtration_ramp_rate", "Ramp rate", "bar/min")
    dm.register_data_key("total_cycle_time", "Total cycle time", "min")
    dm.register_data_key("overall_rejection", "Overall rejection", "%")
    dm.register_data_key("cycle_time_ratio", "Cycle time ratio", "%")
    dm.register_data_key("permeate_concentration", "Permeate concentration", "g/L")
    dm.register_data_key("total_permeate_vol", "Total permeate volume", "m^3")
    # dm.register_data_key(
    #     "blocks[0].process.fs.P2.control_volume.properties_out[0.0].pressure",
    #     "Recycle Pump Pressure",
    #     "bar",
    # )

    # RO
    dm.register_data_key("blocks[0].process.fs.RO.area", "Area", "m^2")
    dm.register_data_key("blocks[0].process.fs.RO.area", "Membrane Area", "m^2")
    dm.register_data_key("blocks[0].process.fs.RO.area", "Total Area", "m^2")
    dm.register_data_key("blocks[0].process.fs.RO.length", "Membrane Length", "m")
    dm.register_data_key("blocks[0].process.fs.RO.width", "Membrane Width", "m")
    dm.register_data_key("blocks[0].process.fs.RO.feed_side.velocity[0.0,0.0]", "Inlet Velocity", "m/s")
    dm.register_data_key(
        "blocks[0].process.fs.RO.recovery_vol_phase[0.0,Liq]",
        "Single Pass Recovery",
        "%",
    )

    dm.register_data_key(
        f"blocks[{last_block + n_flushing_points}].process.fs.RO.feed_side.properties[0.0,0.0].conc_mass_phase_comp[Liq,NaCl]",
        "Inlet concentration",
        "g/L",
    )
    dm.register_data_key(
        f"blocks[{last_block + n_flushing_points}].process.fs.RO.feed_side.properties[0.0,0.0].conc_mass_phase_comp[Liq,NaCl]",
        "Initial concentration",
        "g/L",
    )
    # dm.register_data_key(
    #     f"blocks[0].process.fs.RO.feed_side.properties[0.0,0.0].conc_mass_phase_comp[Liq,NaCl]",
    #     "Initial concentration",
    #     "g/L",
    # )
    dm.register_data_key(
        f"blocks[0].process.fs.P1.control_volume.properties_out[0.0].pressure",
        "Initial Pressure",
        "bar",
    )
    dm.register_data_key(
        f"blocks[{last_block}].process.fs.P1.control_volume.properties_out[0.0].pressure",
        "Final Pressure",
        "bar",
    )

    dm.register_data_key(
        f"blocks[{last_block}].process.fs.P1.control_volume.work[0.0]",
        "Pump size",
        "kW",
    )
    dm.register_data_key(
        f"blocks[{last_block}].process.fs.P1.total_power",
        "Pump work",
        "kW",
    )
    dm.register_data_key(
        f"blocks[{last_block}].process.fs.RO.feed_side.properties[0.0,0.0].conc_mass_phase_comp[Liq,NaCl]",
        "Final concentration",
        "g/L",
    )

    dm.load_data()

    ek = dm.get_expression_keys()
    dm.register_expression(
        ek.Recycle_rate / ek.Feed_flow_rate,
        "Recycle ratio",
    )
    dm.register_expression(
        ek.Pump_size / ek.Average_product_flow_rate, "Specific pump size", "kW/(m^3/s)"
    )
    dm.register_expression(
        ek.Pump_work / ek.Average_product_flow_rate, "Specific pump work", "kW/(m^3/s)"
    )
    dm.register_expression(
        ek.Average_product_flow_rate / ek.Area,
        "Flux",
        "L/(m^2*hr)",
    )
    dm.register_expression(
        ek.Area / ek.Average_product_flow_rate,
        "Specific area",
        "m^2/(L/hr)",
    )
    dm.register_expression(ek.OPEX / ek.CAPEX, "OPEX/CAPEX Ratio")
    dm.register_expression(
        ek.OPEX / ek.Average_product_flow_rate, "Specific OPEX", "USD/m^3"
    )
    dm.register_expression(
        ek.CAPEX / ek.Average_product_flow_rate, "Specific CAPEX", "USD/(m^3/s)"
    )

    dm.evaluate_expressions()

    # lets create costing pacakage
    package_ccro = WaterTapCostingPackage(
        costing_block="costing", validation_key="costing.LCOW"
    )
    package_ccro.register_product_flow("avg_product_flow_rate")

    # Lets create our groups
    RO = PsCostingGroup("RO")
    RO.add_unit(
        f"blocks[0].process.fs.RO",  # only adding "block[0] to specify wher capex is, normally acn just say "RO
        capex_keys="capital_cost",
        fixed_opex_keys="fixed_operating_cost",
    )
    feed_pump = PsCostingGroup("Feed pump")
    feed_pump.add_unit(
        f"blocks[{last_block}].process.fs.P1",
        capex_keys="capital_cost",
        flow_keys={"electricity": "total_power"},
    )
    recycle_pump = PsCostingGroup("Recycle pump")
    recycle_pump.add_unit(
        f"blocks[{last_block}].process.fs.P2",
        capex_keys="capital_cost",
        flow_keys={"electricity": "total_power"},
    )
    conduit = PsCostingGroup("Conduit")
    conduit.add_unit(
        "conduit",
        capex_keys="capital_cost",
    )

    cm_recov = PsCostingManager(
        dm, package_ccro, [RO, feed_pump, recycle_pump, conduit]
    )
    cm_recov.build()
    dm.register_expression(
        ek[("costing", "total", "LCOW_opex")] / ek[("costing", "total", "LCOW")],
        "OPEX LCOW Fraction",
    )
    dm.register_expression(
        ek[("costing", "total", "LCOW_capex")] / ek[("costing", "total", "LCOW")],
        "CAPEX LCOW Fraction",
    )
    dm.register_expression(
        1 / ek.Flux,
        "Specific Membrane Area",
    )
    return dm



def get_ccro_rr_data(last_block=19, n_flushing_points=5):

    last_filt_blk = int(last_block - n_flushing_points)

    dm = PsDataManager()

    dm.register_data_file(
        f"{par_dir}/output/ccro_analysisType_BW_recycle_rate_sweep.h5",
        directory="Brackish water",
    )
    dm.register_data_file(
        f"{par_dir}/output/ccro_analysisType_SW_recycle_rate_sweep.h5",
        directory="Seawater",
    )
    dm.register_data_file(
        f"{par_dir}/output/ccro_analysisType_PW_recycle_rate_sweep.h5",
        directory="Produced water",
    )

    # COSTING
    dm.register_data_key("overall_recovery", "Water recovery", "%")
    dm.register_data_key("costing.LCOW", "LCOW", assign_units="USD/m^3")
    dm.register_data_key("costing.SEC", "SEC", assign_units="kWh/m^3")

    dm.register_data_key(
        "costing.total_capital_cost",
        "CAPEX",
        assign_units="USD",
        conversion_factor=1e-3,
    )
    dm.register_data_key(
        "costing.total_operating_cost",
        "OPEX",
        assign_units="USD/year",
        conversion_factor=1e-3,
    )
    # SYSTEM

    dm.register_data_key("flushing.flushing_efficiency", "Flushing efficiency", "%")
    dm.register_data_key(
        "recycle_loop_concentration", "Recycle loop concentration", "g/L"
    )
    dm.register_data_key(
        "flushing.pre_flushing_concentration", "Pre Flushing Concentration", "g/L"
    )
    dm.register_data_key("recycle_flowrate", "Recycle rate", "L/s")
    # dm.register_data_key(
        # "blocks[0].process.fs.P2.control_volume.properties_out[0.0].flow_vol_phase[Liq]",
    #     "Recycle rate",
    #     "L/s",
    # )
    dm.register_data_key(
        "blocks[0].process.fs.raw_feed.properties[0.0].pressure_osm_phase[Liq]",
        "Feed osmotic pressure",
        "bar",
    )
    dm.register_data_key(
        "blocks[0].process.fs.raw_feed.properties[0.0].flow_vol_phase[Liq]",
        "Feed flow rate",
        "L/s",
    )
    dm.register_data_key("avg_feed_flow_rate", "Average feed flow rate")
    dm.register_data_key("avg_product_flow_rate", "Average product flow rate", "m^3/s")
    # dm.register_data_key("avg_product_flow_rate", "Permeate flow rate", "m^3/s")
    dm.register_data_key("filtration_ramp_rate", "Ramp rate", "bar/min")
    dm.register_data_key("total_cycle_time", "Total cycle time", "min")
    dm.register_data_key("overall_rejection", "Overall rejection", "%")
    dm.register_data_key("cycle_time_ratio", "Cycle time ratio", "%")
    dm.register_data_key("permeate_concentration", "Permeate concentration", "g/L")
    dm.register_data_key("total_permeate_vol", "Total permeate volume", "m^3")
    # dm.register_data_key(
    #     "blocks[0].process.fs.P2.control_volume.properties_out[0.0].pressure",
    #     "Recycle Pump Pressure",
    #     "bar",
    # )

    # RO
    dm.register_data_key("blocks[0].process.fs.RO.area", "Area", "m^2")
    dm.register_data_key("blocks[0].process.fs.RO.area", "Membrane Area", "m^2")
    dm.register_data_key("blocks[0].process.fs.RO.area", "Total Area", "m^2")
    dm.register_data_key("blocks[0].process.fs.RO.length", "Membrane Length", "m")
    dm.register_data_key("blocks[0].process.fs.RO.width", "Membrane Width", "m")
    dm.register_data_key("blocks[0].process.fs.RO.feed_side.velocity[0.0,0.0]", "Inlet Velocity", "cm/s")
    dm.register_data_key(
        "blocks[0].process.fs.RO.recovery_vol_phase[0.0,Liq]",
        "Single Pass Recovery",
        "%",
    )

    dm.register_data_key(
        f"blocks[0].process.fs.RO.feed_side.properties[0.0,0.0].conc_mass_phase_comp[Liq,NaCl]",
        "Inlet concentration",
        "g/L",
    )
    dm.register_data_key(
        f"blocks[0].process.fs.RO.feed_side.properties[0.0,0.0].conc_mass_phase_comp[Liq,NaCl]",
        "Initial concentration",
        "g/L",
    )
    dm.register_data_key(
        f"blocks[0].process.fs.P1.control_volume.properties_out[0.0].pressure",
        "Initial Pressure",
        "bar",
    )
    dm.register_data_key(
        f"blocks[{last_block}].process.fs.P1.control_volume.properties_out[0.0].pressure",
        "Final Pressure",
        "bar",
    )
    dm.register_data_key(
        f"blocks[{last_block}].process.fs.P1.control_volume.properties_out[0.0].pressure",
        "Pressure",
        "bar",
    )
    # dm.register_data_key(
    #     f"blocks[{last_block}].process.fs.P1.control_volume.pressure",
    #     "Pressure",
    #     "bar",
    # )
    dm.register_data_key(
        f"blocks[{last_block}].process.fs.P1.control_volume.work[0.0]",
        "Pump size",
        "kW",
    )
    dm.register_data_key(
        f"blocks[{last_block}].process.fs.P1.total_power",
        "Pump work",
        "kW",
    )
    dm.register_data_key(
        f"blocks[{last_block}].process.fs.RO.feed_side.properties[0.0,0.0].conc_mass_phase_comp[Liq,NaCl]",
        "Final concentration",
        "g/L",
    )

    dm.load_data()

    ek = dm.get_expression_keys()
    dm.register_expression(
        ek.Pump_size / ek.Average_product_flow_rate, "Specific pump size", "kW/(m^3/s)"
    )
    dm.register_expression(
        ek.Pump_work / ek.Average_product_flow_rate, "Specific pump work", "kW/(m^3/s)"
    )
    dm.register_expression(
        ek.Average_product_flow_rate / ek.Area,
        "Flux",
        "L/(m^2*hr)",
    )
    dm.register_expression(
        ek.Area / ek.Average_product_flow_rate,
        "Specific area",
        "m^2/(L/hr)",
    )
    dm.register_expression(ek.OPEX / ek.CAPEX, "OPEX/CAPEX Ratio")
    dm.register_expression(
        ek.OPEX / ek.Average_product_flow_rate, "Specific OPEX", "USD/m^3"
    )
    dm.register_expression(
        ek.CAPEX / ek.Average_product_flow_rate, "Specific CAPEX", "USD/(m^3/s)"
    )

    dm.evaluate_expressions()
    ek = dm.get_expression_keys()
    dm.register_expression(
        ek.Recycle_rate / ek.Feed_flow_rate,
        "Recycle ratio",
    )

    # lets create costing pacakage
    package_ccro = WaterTapCostingPackage(
        costing_block="costing", validation_key="costing.LCOW"
    )
    package_ccro.register_product_flow("avg_product_flow_rate")

    # Lets create our groups
    RO = PsCostingGroup("RO")
    RO.add_unit(
        f"blocks[0].process.fs.RO",  # only adding "block[0] to specify wher capex is, normally acn just say "RO
        capex_keys="capital_cost",
        fixed_opex_keys="fixed_operating_cost",
    )
    feed_pump = PsCostingGroup("Feed pump")
    feed_pump.add_unit(
        f"blocks[{last_block}].process.fs.P1",
        capex_keys="capital_cost",
        flow_keys={"electricity": "total_power"},
    )
    recycle_pump = PsCostingGroup("Recycle pump")
    recycle_pump.add_unit(
        f"blocks[{last_block}].process.fs.P2",
        capex_keys="capital_cost",
        flow_keys={"electricity": "total_power"},
    )
    conduit = PsCostingGroup("Conduit")
    conduit.add_unit(
        "conduit",
        capex_keys="capital_cost",
    )

    cm_recov = PsCostingManager(
        dm, package_ccro, [RO, feed_pump, recycle_pump, conduit]
    )
    cm_recov.build()
    dm.register_expression(
        ek[("costing", "total", "LCOW_opex")] / ek[("costing", "total", "LCOW")],
        "OPEX LCOW Fraction",
    )
    dm.register_expression(
        ek[("costing", "total", "LCOW_capex")] / ek[("costing", "total", "LCOW")],
        "CAPEX LCOW Fraction",
    )
    dm.register_expression(
        1 / ek.Flux,
        "Specific Membrane Area",
    )
    return dm


if __name__ == "__main__":
    # dm = get_ss_data()
    # dm = get_ccro_data()
    # dm = get_ccro_SEC_data()
    dm = get_ccro_rr_data()
    # dm[('costing', 'total', 'LCOW_capex')].data
    # dm.display()
    # ek = dm.get_expression_keys()
    # dm.register_expression(
    #     ek[("costing", "total", "LCOW_opex")]
    #     / ek[("costing", "total", "LCOW")],
    #     "OPEX LCOW Fraction",
    # )
    # dm.register_expression(
    #     ek[("costing", "total", "LCOW_capex")]
    #     / ek[("costing", "total", "LCOW")],
    #     "CAPEX LCOW Fraction",
    # )
    # dm.evaluate_expressions()
    # for k in ek:
    #     print(k)
    # print(dm["Seawater","OPEX LCOW Fraction"].data)
    # print(
    #     dm[
    #         (
    #             "Brackish water",
    #             ("add_erd", "False"),
    #             ("stage_sim_cases", "1_stage_1_pump"),
    #             ("costing", "total", "LCOW_opex"),
    #         )
    #     ].data
    # )
