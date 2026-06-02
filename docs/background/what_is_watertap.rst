What is WaterTAP
------------------

The Water treatment Technoeconomic Assessment Platform (WaterTAP) is a Python-based, open-source library of water treatment models than can be used to assess water treatment trains through simulation, optimization, and other advanced methods.
WaterTAP development is funded by the `National Alliance for Water Innovation (NAWI) <https://www.nawihub.org/>`_, the U.S. Department of Energy’s Energy-Water Desalination Hub.

Motivation
^^^^^^^^^^

NAWI’s objective is to advance early-stage water treatment technologies and enable the treatment of non-traditional water sources. 
As part of these efforts, NAWI is funding the development of WaterTAP to predict the performance and cost of next-generation desalination systems and aid in their design and operation. 
Through the development of WaterTAP, NAWI is seeking to provide the broader water research community an integrated, 
open source modeling and simulation capability to evaluate water treatment options and identify high impact opportunities for innovation including novel materials, processes, and networks.

IDAES
^^^^^

WaterTAP is based on the `Institute of Design of Advanced Energy Systems (IDAES) Integrated Platform <https://idaes.org/>`_, 
which is an advanced process systems engineering tool developed by the U.S. Department of Energy. 
This platform provides significant benefits over traditional modeling approaches because it supports equation-oriented algebraic modeling, 
open-source and commercial solvers, libraries of process units and property packages, and a diverse set of analysis tools,
all within a single fully-featured programming language (Python).

For more information on the IDAES platform and related projects, now called *IDAES+*,
see the `IDAES+ web pages <https://idaesplus.readthedocs.io/latest/>`_.

Advantages
^^^^^^^^^^

WaterTAP has significant advantages over typical modeling approaches including:

* **Modular** – WaterTAP’s modular approach allows users to rapidly assemble components to represent a treatment train and more fully assess the impact of a technology or innovation
* **Multi-hierarchical** – WaterTAP provides models with multiple levels of detail thereby allowing a user to select the appropriate relationships and computational demand for their application
* **Customizable** – WaterTAP allows users to modify the standard models or create custom models to suit their needs
* **IDAES Capabilities** – WaterTAP includes the advantages of the `IDAES Platform <https://idaes-pse.readthedocs.io/en/stable/explanations/why_idaes.html>`_. These advantages include: 1) an equation-oriented approach, which greatly benefits simulation and optimization-based analyses by supporting linear, non-linear, and mixed-integer problems as well as providing access to highly efficient derivative-based solvers; and 2) support for advanced capabilities like dynamics, parameter estimation, conceptual design, surrogate modeling, and uncertainty quantification
* **Open-source** – all WaterTAP code is made freely available for use, modification, and redistribution. WaterTAP's license is located :ref:`here<license:License Agreement>`.
