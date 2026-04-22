.. _ref_ui-fsapi:

WaterTAP User Interface API
===========================
This is the Python API that will be called by the UI backend to process
changes from the UI front-end.

.. include:: <isonum.txt>

.. py:currentmodule:: idaes_flowsheet_processor.api

.. contents:: Contents
    :depth: 1
    :local:

Introduction
------------

This page describes an application programming interface (API) that is designed to help communicate information about a flowsheet and its variables to a user interface layer (or possibly several different kinds of user interface layers at the same time).
The API also provides the ability, intended for the UI developer, to update the variable values and run "actions" such as building and solving the flowsheet.

For how-to guides on specific tasks, see the :ref:`howto_ui-api`.

The entire external-facing API is in :class:`FlowsheetInterface`.
