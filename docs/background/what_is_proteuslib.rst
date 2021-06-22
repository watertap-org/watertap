What is ProteusLib
------------------

ProteusLib is a Python based, open source library of water treatment models than can be used to assess water treatment trains through simulation, optimization, and other advanced methods. 
ProteusLib development is funded by the `National Alliance for Water Innovation (NAWI) <https://www.nawihub.org/>`_, the U.S. Department of Energy’s Energy-Water Desalination Hub.

Motivation
^^^^^^^^^^

NAWI’s objective is to advance early-stage water treatment technologies and enable the treatment of non-traditional water sources. 
As part of these efforts, NAWI is funding the development of Proteuslib to predict the performance and cost of next-generation desalination systems and aid in their design and operation. 
Through the development of Proteuslib, NAWI is seeking to provide the broader water research community an integrated, 
open source modeling and simulation capability to evaluate water treatment options and identify high impact opportunities for innovation including novel materials, processes, and networks.

IDAES
^^^^^

Proteuslib is based on the `Institute of Design of Advanced Energy Systems (IDAES) Integrated Platform <https://idaes.org/>`_, 
which is an advanced process systems engineering tool developed by the U.S. Department of Energy. 
This platform provides significant benefits over traditional modeling approaches because it supports equation-oriented algebraic modeling, 
open source and commercial solvers, libraries of process units and property packages, and a diverse set of analysis tools, 
all within a single fully featured programming language (Python).

Advantages
^^^^^^^^^^

ProteusLib has significant advantages over typical modeling approaches including:

* **Modular** – Proteuslib’s modular approach allows users to rapidly assemble components to represent a treatment train and more fully assess the impact of a technology or innovation
* **Multi-hierarchical** – ProteusLib provides models with multiple levels of detail thereby allowing a user to select the appropriate relationships and computational demand for their application
* **Customizable** – ProteusLib allows users to modify the standard models or create custom models to suit their needs
* **IDAES Capabilities** – ProteusLib includes the advantages of the `IDAES Platform <https://idaes-pse.readthedocs.io/en/stable/user_guide/why_idaes.html>`_. These advantages include: 1) an equation-oriented approach, which greatly benefits simulation and optimization based analyses by supporting linear, non-linear, and mixed-integer problems as well as providing access to highly efficient derivative-based solvers; and 2) support for advanced capabilities like dynamics, parameter estimation, conceptual design, surrogate modeling, and uncertainty quantification
* **Open Source** – all ProteusLib code is made freely available for use, modification, and redistribution. ProteusLib's license is located :ref:`here<license:License Agreement>`.
