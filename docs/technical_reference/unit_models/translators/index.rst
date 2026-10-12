.. _translator_index:

Translators
===========

Translator blocks connects two unit models that use different property or
reaction packages by translating the upstream stream variables (e.g., component flows/concentrations)
into the corresponding variables expected by the downstream property package.

.. toctree::
   :hidden:
   :maxdepth: 1

   translator_adm1_asm1
   translator_adm1_asm2d
   translator_asm1_adm1
   translator_asm2d_adm1

.. csv-table::
   :header: "Unit Model", "Inlet Property Package", "Outlet Property Package", "GitHub Code"
   :widths: 35, 20, 20, 35

   ":ref:`ADM1 to ASM1 Translator <ADM1_ASM1_translator>`", ":ref:`ADM1 <ADM1>`", ":ref:`ASM1 <ASM1>`", "`Translator_ADM1_ASM1 <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/translators/translator_adm1_asm1.py>`_"
   ":ref:`ADM1 to ASM2d Translator <ADM1_ASM2d_translator>`", ":ref:`Modified ADM1 <modified_ADM1>`", ":ref:`Modified ASM2d <modified_ASM2d>`", "`Translator_ADM1_ASM2d <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/translators/translator_adm1_asm2d.py>`_"
   ":ref:`ASM1 to ADM1 Translator <ASM1_ADM1_translator>`", ":ref:`ASM1 <ASM1>`", ":ref:`ADM1 <ADM1>`", "`Translator_ASM1_ADM1 <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/translators/translator_asm1_adm1.py>`_"
   ":ref:`ASM2d to ADM1 Translator <ASM2d_ADM1_translator>`", ":ref:`Modified ASM2d <modified_ASM2d>`", ":ref:`Modified ADM1 <modified_ADM1>`", "`Translator_ASM2d_ASM1 <https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/translators/translator_asm2d_adm1.py>`_"