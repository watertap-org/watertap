import pyomo.environ as pyo
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP

m = pyo.ConcreteModel()
m.x = pyo.Var()
m.c = pyo.Constraint(expr=(0, m.x, 1))
m.o = pyo.Objective(expr=m.x)

nlp = PyomoNLP(m)
