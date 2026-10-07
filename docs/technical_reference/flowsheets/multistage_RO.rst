.. _multistage_RO_flowsheet:

Multi-Stage Reverse Osmosis
===========================

Introduction
------------

This flowsheet can be used to create a N-stage reverse osmosis (RO) system with different configuration options available to the user.

Implementation
--------------

This flowsheet will build, solve, and optimize an N-stage RO system with the following property models:

Property models:

- :ref:`Seawater <seawater>`
- :ref:`NaCl with temperature dependence <nacl_t_dependent>`

Unit models:

- :doc:`Feed model <idaes:reference_guides/model_libraries/generic/unit_models/product>`
- :ref:`Pump model <pump>`
- :ref:`1D RO model <RO_1D>`
- :ref:`ERD model <ERD>`

The functions contained in this flowsheet are designed to allow the user to specify the build of their flowsheet, including:

- The number of stages
- Whether or not to include booster pumps between stages
- Where in the system to place booster pumps
- To include an ERD or not

