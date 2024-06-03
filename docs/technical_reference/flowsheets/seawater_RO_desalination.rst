Seawater RO Desalination
========================

Introduction
------------

This flowsheet represents a full-scale seawater reverse osmosis treatment facility.
The flowsheet includes pretreatment, desalination, post-treatment, and waste handling unit processes available in WaterTAP.
Unit models used on this flowsheet span the spectrum from simple (zero order models *add link*) to complex (reverse osmosis).
This flowsheet includes examples of several different types of modeling features available in WaterTAP, including:

* Zero order models
* Zero order property package
* Reverse osmosis and pump models
* Sodium chloride property package
* Translator blocks
* WaterTAP costing package
* Zero order costing package
* Unit model costing packages
* WaterTAP database
* Several base IDAES models

**add inlet conditions.**

Implementation
--------------

The demonstration file itself contains several core functions that are used to build, specify, initialize, and solve the model, as well as
some helper functions that group these core functions together for convenience. Building and solving the flowsheet proceeds in six steps:

1. Creating and instantiating the model using ``build()``:

    This function will create the core components and structure of the flowsheet. 
    First, the ``FlowsheetBlock``, ``Database``, property models, are created. The zero order property models require the user
    to provide a ``solute_list`` (e.g., TDS and TSS), while the seawater property model is pre-populated with TDS as the only solute.
    Separate ``Block`` are created to contain all the unit models required to model the pretreatment, desalination, and post-treatment
    parts of the treatment train:

        * Pre-treatment (``m.fs.pretreatment``): comprised entirely of zero order models, this block contains the intake, chemical addition, media, and cartridge filtration unit models.
        * Desalination (``m.fs.desalination``): contains all the unit models needed to represent the pumping, reverse osmosis (RO), and energy recovery device (ERD) processes.
        * Post-treatment (``m.fs.posttreatment``): includes post desalination disinfection, remineralization, and storage unit models.

    ``Translator`` blocks are added with appropriate constraints and ``Arc`` are used to connect the unit processes in the proper order. 
    Finally, default scaling factors are set and scaling factors are calculated for all variables.



The majority of the unit models on the flowsheet are split up on to separate ``Blocks`` based on the purpose the treatment process serves.
On the main flowsheet block, there are three separate blocks:



Outside of these blocks, the flowsheet contains a feed, disposal, and product models, as well as ``Translator`` blocks that are used to 
allow for models that use different property models to exchange information. Each of these blocks will be describe further in sections below. 

Pre-Treatment
^^^^^^^^^^^^^


.. figure:: ../../_static/flowsheets/sw_fs_pretreat.png
    :width: 800
    :align: center


Desalination
^^^^^^^^^^^^

Desalination!

Post-Treatment
^^^^^^^^^^^^^^

Post treatment!