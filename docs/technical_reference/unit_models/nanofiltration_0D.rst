Nanofiltration (0D)
===================

Introduction
------------

Nanofiltration (NF) is an important membrane-based water filtration technology. Generally, the pores found in nanofiltration membranes are smaller than those in micro- or ultrafiltration. Unlike reverse osmosis, however, nanofiltration membranes do have actual physical pores (1-10 nanometers in dimension). 

One of the intriguing features of NF membranes is that there are three modes of solute transports that need to be considered. Firstly, the pores in NF membranes allow for convective transport of sufficiently small solutes, just as in micro- and ultrafiltration. Moreover, molecules are able to pass through the membranes through diffusion due to concentration gradients, just as in reverse osmosis. Finally, solute transport across NF membranes is governed by electromigration, i.e., the attraction and/or repulsion of charged solutes based on charges within the membrane itself. 

As such, NF sits at the phenomenological intersection of micro-/ultrafiltration and reverse osmosis. NF membranes are particularly effective at rejecting charged solutions which explains why they are often deployed to lower the concentration of scale-forming (and multivalent) elements such as Calcium or Magnesium. Generally speaking, NF is most often used to treat solutions with fairly low TDS concentrations (i.e. <30k mg/L TDS). 

Modeling the rejection of NF is challenging since it is not always clear how the individual modes of transport described above contribute to the overall membrane performance. That being said, several NF transport models have been proposed over the years. For this implementation, we choose the Kedem-Katchalsky model. 

The main assumptions of the implemented model are as follows:

1) Nanofiltration is modeled as cross-flow filtration as seen in Figure 1
2) Model dimensionality is limited to 0D
3) Model dynamics are restricted to steady-state only
4) Single liquid phase only
5) Single solute only (e.g. NaCl)
6) Isothermal operation

.. figure:: ../../_static/unit_models/nanofiltration.png
    :width: 400
    :align: center
    
    Figure 1. Schematic representation of a nanofiltration unit modeled in IDAES

Ports
---------

The model provides three ports (Pyomo notation in parenthesis):

* Inlet port (inlet)
* Permeate port (permeate)
* Brine port (retentate)

Parameters
----------

The following model parameters need to be specified (Pyomo notation in parenthesis):

* Water permeability coefficient (A)
* Salt permeability coefficient (B)
* Water density (dens_H2O)
* Reflection coefficient (sigma)
* Membrane area (area)
* Pressure drop across the unit (deltaP)

Variables
----------

The model inputs including the following (Pyomo notation in parenthesis):

* Feed components mass flow  (inlet: flow_mass_comp)
* Feed stream temperature (inlet: temperature)
* Feed stream pressure (inlet: pressure)
* Permeate pressure (permeate: pressure)

The model outputs are as follows (Pyomo notation in parenthesis):

* Permeate components mass flow (permeate: flow_mass_comp)
* Permeate stream temperature (permeate: temperature)
* Permeate stream pressure (permeate: pressure)
* Brine stream components mass flow (retentate: flow_mass_comp)
* Brine stream temperature (retentate: temperature)
* Brine stream pressure (retentate: pressure) 

Performance Equations
---------------------

The following equations govern the performance of the nanofiltration unit:

Inlet water and salt fluxes:

.. math::
  J_{in,t,j}=A_t\cdot\rho^{H2O}\cdot(P_{F,t}-P_{P,t}\ )-\sigma\cdot(\pi_{F,t}-\pi_{P,t}\ )
.. math::
  J_{in,t,j}=B_t\cdot(C_{F,t,j}-C_{P,t,j}\ )+(1-\sigma)\cdot J_{in,t,j}\cdot 1/\rho^{H2O}\ \cdot \tilde{c}_{in,t}

Outlet water and salt fluxes:

.. math::
  J_{out,j}=A_t\cdot\rho^{H2O}\cdot(P_{B,t}-P_{P,t})-\sigma\cdot(\pi_{B,t}-\pi_{P,t})
.. math::
  J_{in,t,j}=B_t\cdot(C_{F,t,j}-C_{P,t,j})+(1-\sigma)\cdot J_{out,t,j}\cdot 1/\rho^{H2O} \cdot \tilde{c}_{out,t}

Average flux:

.. math::
  J_{Avg,t,j}=0.5\cdot(J_{in,t,j}+J_{out,t,j})

Average inlet and outlet concentrations:

.. math::
  \tilde{c}_{in,t}=(C_{F,t,j}\cdot C_{P,t,j}\cdot(C_{F,t,j}+C_{P,t,j})/2)^{1/3}
.. math::
  \tilde{c}_{out,t}=(C_{B,t,j}\cdot C_{P,t,j}\cdot(C_{B,t,j}+C_{P,t,j})/2)^{1/3}

Permeate mass flow:

.. math::
  M_{P,t,j}=x_A \cdot J_{Avg,t,j}

