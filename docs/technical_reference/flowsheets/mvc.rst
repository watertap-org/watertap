Mechanical Vapor Compression
============================

Introduction
------------

Mechanical Vapor Compression (MVC) desalination is an evaporative, compressor-driven process where vapor generated
from saline feedwater is compressed to raise its temperature and pressure. This compressed vapor then condenses,
releasing latent heat to evaporate the feedwater, thus producing distillate without external heating or cooling.
The system primarily consists of an evaporator/condenser and a compressor, with optional heat exchangers for enhanced
heat integration, and may utilize external steam to initiate the process.

Implementation
--------------

Figure 1 illustrates the process flow diagram for a Mechanical Vapor Compression (MVC) desalination system.
In this configuration, feedwater is pumped and divided between the cold inlets of a distillate heat
exchanger (HEX) and a brine HEX. The preheated feed exiting these heat exchangers combines and enters the evaporator,
where it absorbs heat from the condensing steam and vaporizes, leaving behind concentrated brine. The generated vapor
is then directed to a compressor, which elevates its temperature and pressure. This compressed vapor is recycled back
to the evaporator, where it condenses and provides the necessary heat for evaporating the preheated feed. The hot
concentrated brine and distillate pass through heat exchangers to preheat the incoming feedwater before exiting the
system. The flowsheet relies on the following key assumptions:

   * supports steady-state only
   * property package(s) supporting liquid and vapors is provided

.. figure:: ../../_static/flowsheets/mvc.png
    :width: 600
    :align: center

    Figure 1. MVC flowsheet

Documentation for each of the WaterTAP unit models can be found below.
    * `Pressure Changer <https://idaes-pse.readthedocs.io/en/latest/reference_guides/model_libraries/generic/unit_models/pressure_changer.html>`_
    * `Evaporator <https://watertap.readthedocs.io/en/latest/technical_reference/unit_models/mvc.html>`_
    * `Compressor <https://watertap.readthedocs.io/en/latest/technical_reference/unit_models/mvc.html>`_
    * `Condenser <https://watertap.readthedocs.io/en/latest/technical_reference/unit_models/mvc.html>`_

Documentation for each of the IDAES unit models can be found below.
    * `Feed <https://idaes-pse.readthedocs.io/en/latest/reference_guides/model_libraries/generic/unit_models/feed.html>`_
    * `Heat Exchanger <https://idaes-pse.readthedocs.io/en/latest/reference_guides/model_libraries/generic/unit_models/heat_exchanger.html>`_
    * `Separator <https://idaes-pse.readthedocs.io/en/latest/reference_guides/model_libraries/generic/unit_models/separator.html>`_
    * `Product <https://idaes-pse.readthedocs.io/en/latest/reference_guides/model_libraries/generic/unit_models/product.html>`_
    * `Mixer <https://idaes-pse.readthedocs.io/en/latest/reference_guides/model_libraries/generic/unit_models/mixer.html>`_
    * `Translator <https://idaes-pse.readthedocs.io/en/latest/reference_guides/model_libraries/generic/unit_models/translator.html>`_

Documentation for each of the property models can be found below.
    * `Water <https://watertap.readthedocs.io/en/latest/technical_reference/property_models/water.html>`_
    * `Seawater <https://watertap.readthedocs.io/en/latest/technical_reference/property_models/seawater.html>`_

Documentation for the costing relationships can be found below.
    * `WaterTAP Costing Package <https://watertap.readthedocs.io/en/latest/technical_reference/costing/watertap_costing.html>`_

The objective function is to minimize the levelized cost of water, which can be represented by the following equation
where :math:`Q` represents volumetric flow, :math:`f_{crf}` represents capital recovery factor
:math:`C_{cap,tot}` represents total capital cost, :math:`C_{op,tot}` represents total operating cost, and
:math:`f_{util}` represents the utilization factor:

    .. math::

        LCOW_{Q} = \frac{f_{crf}   C_{cap,tot} + C_{op,tot}}{f_{util} Q}

Degrees of Freedom
------------------
The following variables are specified for MVC flowsheet based on the default settings:
    * feed water conditions (mass flow, mass fractions, temperature, and pressure)
    * feed pump efficiency and pressure change (ΔP)
    * distillate HEX heat transfer coefficient, cold-side ΔP, and hot-side ΔP
    * brine HEX heat transfer coefficient, cold-side ΔP, and hot-side ΔP
    * evaporator overall heat transfer coefficient
    * compressor efficiency
    * brine pump efficiency and ΔP
    * distillate pump efficiency and ΔP
    * translator block outlet TDS concentration

Flowsheet Specifications
------------------------

.. csv-table::
   :header: "Description", "Value", "Units"

   "**Feed Water**"
   "Water mass flow","40", ":math:`\text{kg/s}`"
   "TDS mass fraction", "0.1", ":math:`\text{dimensionless}`"
   "Temperature", "298.15", ":math:`\text{K}`"
   "Pressure", "101325", ":math:`\text{Pa}`"

   "**Feed Pump**"
   "Pump efficiency", "0.8", ":math:`\text{dimensionless}`"
   "Pressure change", "7000", ":math:`\text{Pa}`"

   "**Separator**"
   "Total flow split fraction to distillate HEX", "0.5", ":math:`\text{dimensionless}`"

   "**Distillate HEX**"
   "Overall heat transfer coefficient", "2000", ":math:`W/\left(m^2K\right)`"
   "Area", "125", ":math:`\text{m}^2`"
   "Cold-side pressure change", "7000", ":math:`\text{Pa}`"
   "Hot-side pressure change", "7000", ":math:`\text{Pa}`"

   "**Brine HEX**"
   "Overall heat transfer coefficient", "2000", ":math:`W/\left(m^2K\right)`"
   "Area", "125", ":math:`\text{m}^2`"
   "Cold-side pressure change", "7000", ":math:`\text{Pa}`"
   "Hot-side pressure change", "7000", ":math:`\text{Pa}`"

   "**Evaporator**"
   "Outlet brine temperature", "343.15", ":math:`\text{K}`"
   "Overall heat transfer coefficient", "3000", ":math:`W/\left(m^2K\right)`"

   "**Compressor**"
   "Compressor efficiency", "0.8", ":math:`\text{dimensionless}`"
   "Pressure ratio", "1.6", ":math:`\text{dimensionless}`"

   "**Brine Pump**"
   "Pump efficiency", "0.8", ":math:`\text{dimensionless}`"
   "Pressure change", "40000", ":math:`\text{Pa}`"

   "**Distillate Pump**"
   "Pump efficiency", "0.8", ":math:`\text{dimensionless}`"
   "Pressure change", "40000", ":math:`\text{Pa}`"

   "**Translator Block**"
   "Outlet TDS mass flow", "1e-5", ":math:`\text{kg/s}`"
