from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock

import ion_prop_pack as props

# create model, flowsheet
m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
# attach property package
m.fs.properties = props.PropParameterBlock(default={'ion_list': ['Na_+', 'Cl_-', 'Mg_2+']})
# build a state block, must specify a time which by convention for steady state models is just 0
m.fs.stream = m.fs.properties.build_state_block([0], default={})

# display the state block, it only has the state variables and they are all unfixed
print('\n---first display---')
m.fs.stream[0].display()
