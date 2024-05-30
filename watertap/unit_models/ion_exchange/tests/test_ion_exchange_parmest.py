import os
import matplotlib.pyplot as plt
from watertap.unit_models import IXParmest

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
data_file = __location__ + "/ix_parmest_test_data.csv"

# Data required
resin_data = dict(diameter=675e-6, density=0.72, porosity=0.4)
compound_data = dict(name="pfos", mw_comp=0.41407, diffusivity=0.49e-9, charge=-1)
# Initial guess for thetas
initial_guess_dict = dict(mass_transfer_coeff=1, freundlich_n=1.1)

set_bounds_dict = {
    "mass_transfer_coeff": [0, None],
    "freundlich_n": [1.05, 100],
    "service_flow_rate": [0, None],
    "loading_rate": [0, None],
    "ebct": [0, None],
    "resin_diam": [0, None],
    "bv": [0, None],
    "bed_depth": [0, None],
    "bed_diameter": [0, None],
}


ixp = IXParmest(
    data_file=data_file,
    initial_guess_dict=initial_guess_dict,
    resin_data=resin_data,
    compound_data=compound_data,
    set_bounds_dict=set_bounds_dict,
)

ixp.test_initial_guess()
ixp.estimate_bv50()
ixp.plot_initial_guess()
ixp.run_parmest()
ixp.test_theta()
ixp.plot_theta()
plt.show()
