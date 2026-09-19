.. _0D_index:

Zero-Order Unit Models
======================

The zero-order models rely on default model parameter values specified in YAML files, but these data are generally only meant to help initiate model usage for new users. Users should supply their own values, if possible, instead of relying on the default parameter values. The YAML database files for the zero-order models are located in :code:`watertap/data/techno_economic/`. The name of the YAML file should match the associated model name without the ZO letters and use snake case style. For example, the YAML file for the :code:`DualMediaFiltrationZO()` zero-order model is named :code:`dual_media_filtration.yaml`.

.. toctree::
   :hidden:
   :maxdepth: 1

   aeration_basin_zo
   anaerobic_digestion_oxidation_zo
   anaerobic_digestion_reactive_zo
   anaerobic_mbr_mec_zo
   autothermal_hydrothermal_liquefaction_zo
   backwash_solids_handling_zo
   bio_active_filtration_zo
   bioreactor_zo
   blending_reservoir_zo
   brine_concentrator_zo
   buffer_tank_zo
   CANDOP_zo
   cartridge_filtration_zo
   centrifuge_zo
   chemical_addition_zo
   chlorination_zo
   clarifier_zo
   cloth_media_filtration_zo
   co2_addition_zo
   coag_and_floc_zo
   cofermentation_zo
   constructed_wetlands_zo
   conventional_activated_sludge_zo
   cooling_supply_zo
   cooling_tower_zo
   crystallizer_zo
   decarbonator_zo
   deep_well_injection_zo
   dissolved_air_flotation_zo
   dmbr_zo
   dual_media_filtration_zo
   electrochemical_nutrient_removal_zo
   electrocoagulation_zo
   electrodialysis_reversal_zo
   energy_recovery_zo
   evaporation_pond_zo
   feed_water_tank_zo
   feed_zo
   filter_press_zo
   fixed_bed_zo
   gac_zo
   gas_sparged_membrane_zo
   hrcs_zo
   hydrothermal_gasification_zo
   injection_well_disposal_zo
   intrusion_mitigation_zo
   ion_exchange_zo
   iron_and_manganese_removal_zo
   landfill_zo
   mabr_zo
   magprex_zo
   mbr_zo
   media_filtration_zo
   membrane_evaporator_zo
   metab_zo
   microbial_battery_zo
   microfiltration_zo
   municipal_drinking_zo
   municipal_wwtp_zo
   nanofiltration_zo
   ozone_aop_zo
   ozone_zo
   peracetic_acid_disinfection_zo
   photothermal_membrane_zo
   primary_separator_zo
   pump_electricity_zo
   pump_zo
   screen_zo
   secondary_treatment_wwtp_zo
   sedimentation_zo
   settling_pond_zo
   sludge_tank_zo
   smp_zo
   static_mixer_zo
   storage_tank_zo
   struvite_classifier_zo
   suboxic_activated_sludge_process_zo
   supercritical_salt_precipitation_zo
   surface_discharge_zo
   sw_onshore_intake_zo
   tramp_oil_tank_zo
   tri_media_filtration_zo
   ultra_filtration_zo
   uv_aop_zo
   uv_zo
   vfa_recovery_zo
   waiv_zo
   walnut_shell_filter_zo
   water_pumping_station_zo
   well_field_zo

.. csv-table::
   :header: "Unit Model", "Compatible Property Packages", "GitHub Code"
   :widths: 35, 30, 35

   ":ref:`Aeration Basin <aeration_basin_zo>`", "Any", "`AerationBasinZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/aeration_basin_zo.py>`_"
   ":ref:`Anaerobic Digestion Oxidation <anaerobic_digestion_oxidation_zo>`", "Any", "`AnaerobicDigestionOxidationZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/anaerobic_digestion_oxidation_zo.py>`_"
   ":ref:`Anaerobic Digestion Reactive <anaerobic_digestion_reactive_zo>`", "Any", "`AnaerobicDigestionReactiveZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/anaerobic_digestion_reactive_zo.py>`_"
   ":ref:`Anaerobic MBR MEC <anaerobic_mbr_mec_zo>`", "Any", "`AnaerobicMBRMECZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/anaerobic_mbr_mec_zo.py>`_"
   ":ref:`Autothermal Hydrothermal Liquefaction <autothermal_hydrothermal_liquefaction_zo>`", "Any", "`AutothermalHydrothermalLiquefactionZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/autothermal_hydrothermal_liquefaction_zo.py>`_"
   ":ref:`Backwash Solids Handling <backwash_solids_handling_zo>`", "Any", "`BackwashSolidsHandlingZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/backwash_solids_handling_zo.py>`_"
   ":ref:`Bio Active Filtration <bio_active_filtration_zo>`", "Any", "`BioActiveFiltrationZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/bio_active_filtration_zo.py>`_"
   ":ref:`Bioreactor <bioreactor_zo>`", "Any", "`BioreactorZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/bioreactor_zo.py>`_"
   ":ref:`Blending Reservoir <blending_reservoir_zo>`", "Any", "`BlendingReservoirZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/blending_reservoir_zo.py>`_"
   ":ref:`Brine Concentrator <brine_concentrator_zo>`", "Any", "`BrineConcentratorZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/brine_concentrator_zo.py>`_"
   ":ref:`Buffer Tank <buffer_tank_zo>`", "Any", "`BufferTankZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/buffer_tank_zo.py>`_"
   ":ref:`CANDOP <CANDOP_zo>`", "Any", "`CANDOPZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/CANDOP_zo.py>`_"
   ":ref:`Cartridge Filtration <cartridge_filtration_zo>`", "Any", "`CartridgeFiltrationZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/cartridge_filtration_zo.py>`_"
   ":ref:`Centrifuge <centrifuge_zo>`", "Any", "`CentrifugeZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/centrifuge_zo.py>`_"
   ":ref:`Chemical Addition <chemical_addition_zo>`", "Any", "`ChemicalAdditionZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/chemical_addition_zo.py>`_"
   ":ref:`Chlorination <chlorination_zo>`", "Any", "`ChlorinationZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/chlorination_zo.py>`_"
   ":ref:`Clarifier <clarifier_zo>`", "Any", "`ClarifierZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/clarifier_zo.py>`_"
   ":ref:`Cloth Media Filtration <cloth_media_filtration_zo>`", "Any", "`ClothMediaFiltrationZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/cloth_media_filtration_zo.py>`_"
   ":ref:`CO2 Addition <co2_addition_zo>`", "Any", "`CO2AdditionZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/co2_addition_zo.py>`_"
   ":ref:`Coagulation and Flocculation <coag_and_floc_zo>`", "Any", "`CoagulationFlocculationZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/coag_and_floc_zo.py>`_"
   ":ref:`Cofermentation <cofermentation_zo>`", "Any", "`CofermentationZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/cofermentation_zo.py>`_"
   ":ref:`Constructed Wetlands <constructed_wetlands_zo>`", "Any", "`ConstructedWetlandsZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/constructed_wetlands_zo.py>`_"
   ":ref:`Conventional Activated Sludge <conventional_activated_sludge_zo>`", "Any", "`ConventionalActivatedSludgeZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/conventional_activated_sludge_zo.py>`_"
   ":ref:`Cooling Supply <cooling_supply_zo>`", "Any", "`CoolingSupplyZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/cooling_supply_zo.py>`_"
   ":ref:`Cooling Tower <cooling_tower_zo>`", "Any", "`CoolingTowerZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/cooling_tower_zo.py>`_"
   ":ref:`Crystallizer <crystallizer_zo>`", "Any", "`CrystallizerZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/crystallizer_zo.py>`_"
   ":ref:`Decarbonator <decarbonator_zo>`", "Any", "`DecarbonatorZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/decarbonator_zo.py>`_"
   ":ref:`Deep Well Injection <deep_well_injection_zo>`", "Any", "`DeepWellInjectionZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/deep_well_injection_zo.py>`_"
   ":ref:`Dissolved Air Flotation <dissolved_air_flotation_zo>`", "Any", "`DissolvedAirFlotationZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/dissolved_air_flotation_zo.py>`_"
   ":ref:`DMBR <dmbr_zo>`", "Any", "`DMBRZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/dmbr_zo.py>`_"
   ":ref:`Dual Media Filtration <dual_media_filtration_zo>`", "Any", "`DualMediaFiltrationZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/dual_media_filtration_zo.py>`_"
   ":ref:`Electrochemical Nutrient Removal <electrochemical_nutrient_removal_zo>`", "Any", "`ElectrochemicalNutrientRemovalZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/electrochemical_nutrient_removal_zo.py>`_"
   ":ref:`Electrocoagulation <electrocoagulation_zo>`", "Any", "`ElectrocoagulationZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/electrocoagulation_zo.py>`_"
   ":ref:`Electrodialysis Reversal <electrodialysis_reversal_zo>`", "Any", "`ElectrodialysisReversalZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/electrodialysis_reversal_zo.py>`_"
   ":ref:`Energy Recovery <energy_recovery_zo>`", "Any", "`EnergyRecoveryZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/energy_recovery_zo.py>`_"
   ":ref:`Evaporation Pond <evaporation_pond_zo>`", "Any", "`EvaporationPondZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/evaporation_pond_zo.py>`_"
   ":ref:`Feed Water Tank <feed_water_tank_zo>`", "Any", "`FeedWaterTankZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/feed_water_tank_zo.py>`_"
   ":ref:`Feed <feed_zo>`", "Any", "`FeedZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/feed_zo.py>`_"
   ":ref:`Filter Press <filter_press_zo>`", "Any", "`FilterPressZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/filter_press_zo.py>`_"
   ":ref:`Fixed Bed <fixed_bed_zo>`", "Any", "`FixedBedZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/fixed_bed_zo.py>`_"
   ":ref:`GAC <gac_zo>`", "Any", "`GACZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/gac_zo.py>`_"
   ":ref:`Gas-Sparged Membrane <gas_sparged_membrane_zo>`", "Any", "`GasSpargedMembraneZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/gas_sparged_membrane_zo.py>`_"
   ":ref:`HRCS <hrcs_zo>`", "Any", "`HRCSZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/hrcs_zo.py>`_"
   ":ref:`Hydrothermal Gasification <hydrothermal_gasification_zo>`", "Any", "`HydrothermalGasificationZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/hydrothermal_gasification_zo.py>`_"
   ":ref:`Injection Well Disposal <injection_well_disposal_zo>`", "Any", "`InjectionWellDisposalZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/injection_well_disposal_zo.py>`_"
   ":ref:`Intrusion Mitigation <intrusion_mitigation_zo>`", "Any", "`IntrusionMitigationZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/intrusion_mitigation_zo.py>`_"
   ":ref:`Ion Exchange <ion_exchange_zo>`", "Any", "`IonExchangeZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/ion_exchange_zo.py>`_"
   ":ref:`Iron and Manganese Removal <iron_and_manganese_removal_zo>`", "Any", "`IronAndManganeseRemovalZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/iron_and_manganese_removal_zo.py>`_"
   ":ref:`Landfill <landfill_zo>`", "Any", "`LandfillZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/landfill_zo.py>`_"
   ":ref:`MABR <mabr_zo>`", "Any", "`MABRZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/mabr_zo.py>`_"
   ":ref:`MAGPREX <magprex_zo>`", "Any", "`MAGPREXZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/magprex_zo.py>`_"
   ":ref:`MBR <mbr_zo>`", "Any", "`MBRZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/mbr_zo.py>`_"
   ":ref:`Media Filtration <media_filtration_zo>`", "Any", "`MediaFiltrationZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/media_filtration_zo.py>`_"
   ":ref:`Membrane Evaporator <membrane_evaporator_zo>`", "Any", "`MembraneEvaporatorZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/membrane_evaporator_zo.py>`_"
   ":ref:`Metab <metab_zo>`", "Any", "`MetabZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/metab_zo.py>`_"
   ":ref:`Microbial Battery <microbial_battery_zo>`", "Any", "`MicrobialBatteryZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/microbial_battery_zo.py>`_"
   ":ref:`Microfiltration <microfiltration_zo>`", "Any", "`MicrofiltrationZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/microfiltration_zo.py>`_"
   ":ref:`Municipal Drinking Water <municipal_drinking_zo>`", "Any", "`MunicipalDrinkingZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/municipal_drinking_zo.py>`_"
   ":ref:`Municipal WWTP <municipal_wwtp_zo>`", "Any", "`MunicipalWWTPZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/municipal_wwtp_zo.py>`_"
   ":ref:`Nanofiltration <nanofiltration_zo>`", "Any", "`NanofiltrationZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/nanofiltration_zo.py>`_"
   ":ref:`Ozone AOP <ozone_aop_zo>`", "Any", "`OzoneAOPZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/ozone_aop_zo.py>`_"
   ":ref:`Ozone <ozone_zo>`", "Any", "`OzoneZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/ozone_zo.py>`_"
   ":ref:`Peracetic Acid Disinfection <peracetic_acid_disinfection_zo>`", "Any", "`PeraceticAcidDisinfectionZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/peracetic_acid_disinfection_zo.py>`_"
   ":ref:`Photothermal Membrane <photothermal_membrane_zo>`", "Any", "`PhotothermalMembraneZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/photothermal_membrane_zo.py>`_"
   ":ref:`Primary Separator <primary_separator_zo>`", "Any", "`PrimarySeparatorZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/primary_separator_zo.py>`_"
   ":ref:`Pump Electricity <pump_electricity_zo>`", "Any", "`PumpElectricityZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/pump_electricity_zo.py>`_"
   ":ref:`Pump <pump_zo>`", "Any", "`PumpZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/pump_zo.py>`_"
   ":ref:`Screen <screen_zo>`", "Any", "`ScreenZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/screen_zo.py>`_"
   ":ref:`Secondary Treatment WWTP <secondary_treatment_wwtp_zo>`", "Any", "`SecondaryTreatmentWWTPZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/secondary_treatment_wwtp_zo.py>`_"
   ":ref:`Sedimentation <sedimentation_zo>`", "Any", "`SedimentationZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/sedimentation_zo.py>`_"
   ":ref:`Settling Pond <settling_pond_zo>`", "Any", "`SettlingPondZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/settling_pond_zo.py>`_"
   ":ref:`Sludge Tank <sludge_tank_zo>`", "Any", "`SludgeTankZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/sludge_tank_zo.py>`_"
   ":ref:`SMP <smp_zo>`", "Any", "`SMPZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/smp_zo.py>`_"
   ":ref:`Static Mixer <static_mixer_zo>`", "Any", "`StaticMixerZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/static_mixer_zo.py>`_"
   ":ref:`Storage Tank <storage_tank_zo>`", "Any", "`StorageTankZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/storage_tank_zo.py>`_"
   ":ref:`Struvite Classifier <struvite_classifier_zo>`", "Any", "`StruviteClassifierZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/struvite_classifier_zo.py>`_"
   ":ref:`Suboxic Activated Sludge Process <suboxic_activated_sludge_process_zo>`", "Any", "`SuboxicActivatedSludgeProcessZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/suboxic_activated_sludge_process_zo.py>`_"
   ":ref:`Supercritical Salt Precipitation <supercritical_salt_precipitation_zo>`", "Any", "`SupercriticalSaltPrecipitationZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/supercritical_salt_precipitation_zo.py>`_"
   ":ref:`Surface Discharge <surface_discharge_zo>`", "Any", "`SurfaceDischargeZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/surface_discharge_zo.py>`_"
   ":ref:`SW Onshore Intake <sw_onshore_intake_zo>`", "Any", "`SWOnshoreIntakeZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/sw_onshore_intake_zo.py>`_"
   ":ref:`Tramp Oil Tank <tramp_oil_tank_zo>`", "Any", "`TrampOilTankZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/tramp_oil_tank_zo.py>`_"
   ":ref:`Tri-Media Filtration <tri_media_filtration_zo>`", "Any", "`TriMediaFiltrationZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/tri_media_filtration_zo.py>`_"
   ":ref:`Ultra Filtration <ultra_filtration_zo>`", "Any", "`UltraFiltrationZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/ultra_filtration_zo.py>`_"
   ":ref:`UV AOP <uv_aop_zo>`", "Any", "`UVAOPZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/uv_aop_zo.py>`_"
   ":ref:`UV <uv_zo>`", "Any", "`UVZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/uv_zo.py>`_"
   ":ref:`VFA Recovery <vfa_recovery_zo>`", "Any", "`VFARecoveryZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/vfa_recovery_zo.py>`_"
   ":ref:`Wind-Aided Intensified Evaporation Unit <waiv_zo>`", "Any", "`WAIVZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/waiv_zo.py>`_"
   ":ref:`Walnut Shell Filter <walnut_shell_filter_zo>`", "Any", "`WalnutShellFilterZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/walnut_shell_filter_zo.py>`_"
   ":ref:`Water Pumping Station <water_pumping_station_zo>`", "Any", "`WaterPumpingStationZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/water_pumping_station_zo.py>`_"
   ":ref:`Well Field <well_field_zo>`", "Any", "`WellFieldZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/zero_order/well_field_zo.py>`_"
