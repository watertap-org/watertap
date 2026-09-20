Unit Models
===========

.. toctree::
   :hidden:
   :maxdepth: 1

   anaerobic_digester
   aeration_tank
   boron_removal
   clarifier
   coag_floc_model
   crystallizer_0D
   cstr
   cstr_injection
   dewatering
   electrocoagulation
   electrodialysis_0D
   electrodialysis_1D
   electrodialysis_bipolar_1D
   electrolyzer
   electroNP_ZO
   energy_recovery_device
   gac
   generic_desalter
   generic_separator
   ion_exchange_0D
   membrane_distillation_0D
   membrane_distillation_1D
   mvc
   nanofiltration_ZO
   nanofiltration_0D
   nanofiltration_dspmde_0D
   osmotically_assisted_reverse_osmosis_0D
   osmotically_assisted_reverse_osmosis_1D
   pressure_exchanger
   pump
   reverse_osmosis_0D
   reverse_osmosis_1D
   steam_ejector
   steam_heater_0D
   stoichiometric_reactor
   thickener
   translators/index
   uv_aop
   zero_order_unit_models/index
   unit_model_utilities

.. csv-table::
   :header: "Unit Model", "Compatible Property Packages", "GitHub Code"
   :widths: 35, 30, 35

   ":ref:`Aeration Tank <aeration_tank>`", ":ref:`ASM1 <ASM1>`, :ref:`ASM2d <ASM2d>`, :ref:`Modified ASM2d <modified_ASM2d>`", "`AerationTank <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/aeration_tank.py>`_"
   ":ref:`Anaerobic Digester <anaerobic_digester>`", ":ref:`ADM1 <ADM1>`, :ref:`Modified ADM1 <modified_ADM1>`", "`AD <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/anaerobic_digester.py>`_"
   ":ref:`Boron Removal <boron_removal>`", ":ref:`MCAS <mcas_tech_ref>`", "`BoronRemoval <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/boron_removal.py>`_"
   ":ref:`Clarifier <clarifier>`", "Any", "`Clarifier <https://github.com/watertap-org/watertap/unit_models/clarifier.py>`_"
   ":ref:`Coagulation/Flocculation Model <coagulation_flocculation>`", ":ref:`Coagulation <coagulation>`", "`CoagulationFlocculation <https://github.com/watertap-org/watertap/unit_models/coag_floc_model.py>`_"
   ":ref:`0D Crystallizer <crystallizer_0D>`", "Crystallization (TBD)", "`Crystallization <https://github.com/watertap-org/watertap/unit_models/crystallizer_0D.py>`_"
   ":ref:`CSTR <CSTR>`", "Any", "`CSTR <https://github.com/watertap-org/watertap/unit_models/cstr.py>`_"
   ":ref:`CSTR Injection <CSTR_injection>`", ":ref:`ASM1 <ASM1>`, :ref:`ASM2d <ASM2d>`, :ref:`Modified ASM2d <modified_ASM2d>`", "`CSTR_Injection <https://github.com/watertap-org/watertap/unit_models/cstr_injection.py>`_"
   ":ref:`Dewatering Unit <dewatering>`", ":ref:`ASM1 <ASM1>`, :ref:`ASM2d <ASM2d>`, :ref:`Modified ASM2d <modified_ASM2d>`", "`DewateringUnit <https://github.com/watertap-org/watertap/unit_models/dewatering.py>`_"
   ":ref:`Electrocoagulation <EC_0D>`", ":ref:`MCAS <mcas_tech_ref>`", "`Electrocoagulation <https://github.com/watertap-org/watertap/unit_models/electrocoagulation.py>`_"
   ":ref:`0D Electrodialysis <ED_0D>`", ":ref:`MCAS <mcas_tech_ref>`", "`Electrodialysis0D <https://github.com/watertap-org/watertap/unit_models/electrodialysis_0D.py>`_"
   ":ref:`1D Electrodialysis <ED_1D>`", ":ref:`MCAS <mcas_tech_ref>`", "`Electrodialysis1D <https://github.com/watertap-org/watertap/unit_models/electrodialysis_1D.py>`_"
   ":ref:`1D Bipolar Electrodialysis <ED_bipolar_1D>`", ":ref:`MCAS <mcas_tech_ref>`", "`Electrodialysis_Bipolar_1D <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/electrodialysis_bipolar_1D.py>`_"
   ":ref:`Electrolyzer <electrolyzer>`", ":ref:`MCAS <mcas_tech_ref>`", "`Electrolyzer <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/electrolyzer.py>`_"
   ":ref:`ElectroNP ZO <electroNP>`", ":ref:`Modified ASM2d <modified_ASM2d>`", "`ElectroNPZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/electroNP_ZO.py>`_"
   ":ref:`Energy Recovery Device <ERD>`", "Any", "`EnergyRecoveryDevice <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/energy_recovery_device.py>`_"
   ":ref:`GAC <GAC>`", ":ref:`MCAS <mcas_tech_ref>`", "`GAC <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/gac.py>`_"
   ":ref:`Generic Desalter <generic_desalter>`", ":ref:`MCAS <mcas_tech_ref>`", "`GenericDesalter <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/generic_desalter.py>`_"
   ":ref:`Generic Separator <generic_separator>`", ":ref:`MCAS <mcas_tech_ref>`", "`GenericSeparator <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/generic_separator.py>`_"
   ":ref:`0D Ion Exchange <IX_0D>`", ":ref:`MCAS <mcas_tech_ref>`", "`IonExchange0D <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/ion_exchange_0D.py>`_"
   ":ref:`0D Membrane Distillation <MD_0D>`", "Any", "`MembraneDistillation0D <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/membrane_distillation_0D.py>`_"
   ":ref:`1D Membrane Distillation <MD_1D>`", "Any", "`MembraneDistillation1D <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/membrane_distillation_1D.py>`_"
   ":ref:`MVC <MVC>`", "Any", "`MVC <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/mvc.py>`_"
   ":ref:`ZO Nanofiltration <NF_ZO>`", ":ref:`MCAS <mcas_tech_ref>`", "`NanofiltrationZO <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/nanofiltration_ZO.py>`_"
   ":ref:`0D Nanofiltration <nanofiltration_0D>`", ":ref:`MCAS <mcas_tech_ref>`", "`Nanofiltration0D <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/nanofiltration_0D.py>`_"
   ":ref:`0D DSPM-DE Nanofiltration <nanofiltration_DSPMDE>`", ":ref:`MCAS <mcas_tech_ref>`", "`NanofiltrationDSPMDE0D <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/nanofiltration_dspmde_0D.py>`_"
   ":ref:`0D Osmotically Assisted Reverse Osmosis <OARO_0D>`", "Any", "`OsmoticallyAssistedReverseOsmosis0D <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/osmotically_assisted_reverse_osmosis_0D.py>`_"
   ":ref:`1D Osmotically Assisted Reverse Osmosis <OARO_1D>`", "Any", "`OsmoticallyAssistedReverseOsmosis1D <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/osmotically_assisted_reverse_osmosis_1D.py>`_"
   ":ref:`Pressure Exchanger <pressure_exchanger>`", "Any", "`PressureExchanger <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/pressure_exchanger.py>`_"
   ":ref:`Pump <pump>`", "Any", "`Pump <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/pump.py>`_"
   ":ref:`0D Reverse Osmosis <RO_0D>`", "Any", "`ReverseOsmosis0D <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/reverse_osmosis_0D.py>`_"
   ":ref:`1D Reverse Osmosis <RO_1D>`", "Any", "`ReverseOsmosis1D <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/reverse_osmosis_1D.py>`_"
   ":ref:`Steam Ejector <steam_ejector>`", ":ref:`Water <water>`", "`SteamEjector <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/steam_ejector.py>`_"
   ":ref:`0D Steam Heater <steam_heater_0D>`", ":ref:`Water <water>`, :ref:`Seawater <seawater>`", "`SteamHeater0D <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/steam_heater_0D.py>`_"
   ":ref:`Stoichiometric Reactor <stoichiometric_reactor>`", ":ref:`MCAS <mcas_tech_ref>`", "`StoichiometricReactor <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/stoichiometric_reactor.py>`_"
   ":ref:`Thickener <thickener>`", ":ref:`ASM1 <ASM1>`, :ref:`ASM2d <ASM2d>`, :ref:`Modified ASM2d <modified_ASM2d>`", "`Thickener <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/thickener.py>`_"
   ":doc:`Translators <translators/index>`", "See :ref:`Translator Index <translator_index>`", "`Translators <https://github.com/watertap-org/watertap/tree/main/watertap/unit_models/translators>`_"
   ":ref:`UV AOP <UV_AOP>`", ":ref:`MCAS <mcas_tech_ref>`, NDMA (TBD)", "`Ultraviolet0D <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/uv_aop.py>`_"
   ":doc:`Zero-Order Unit Models <zero_order_unit_models/index>`", "See :ref:`Zero-Order Unit Model Index <0D_index>`", "`0D Unit Models <https://github.com/watertap-org/watertap/tree/main/watertap/unit_models/zero_order>`_"
   ":doc:`Unit Model Utility Functions <unit_model_utilities>`", "`-`", "`Unit Model Utilities <https://github.com/watertap-org/watertap/blob/main/watertap/core/util/unit_models.py>`_"