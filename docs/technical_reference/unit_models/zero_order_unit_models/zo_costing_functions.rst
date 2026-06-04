ZO Costing Functions
===================

Costing functions for zero-order unit models that do not follow power law..

ozone_zo
--------
.. start_ozone_zo_costing

Ozone capital cost is calculated using the following equation:

.. math:: 
   C_{capital} = Cost\ Factor * (A + B * dose_{ozone} + C * ln(Q) + D * ln(Q) * dose_{ozone})

where:

- :math:`C_{capital}` is the capital cost of the ozone contactor, in dollars

- :math:`Cost\ Factor` is a scaling factor to adjust the cost based on the size of the system, in dollars.

- :math:`A`, :math:`B`, :math:`C`, and :math:`D` are coefficients determined from cost data, dimensionless.

- :math:`dose_{ozone}` is the amount of ozone consumption, in mg/L.

- :math:`Q` is the volumetric flow rate of the influent water, in m^{3}/h.

.. end_ozone_zo_costing
