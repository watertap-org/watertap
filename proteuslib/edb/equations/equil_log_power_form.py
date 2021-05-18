'''
    This file creates a custom equilibrium constraint form for evaluating
    that constraint in terms of the log (base 10) form. This is necessary
    as the typical equilibrium reactions in solution rely on very small
    molar concentrations, which make convergence very difficult if we
    do not use this equation form to appropriately scale the constraint.
'''

# Import the 'get_concentration_term' function from the idaes framework
from idaes.generic_models.properties.core.generic.generic_reaction import \
    get_concentration_term

# Import log10 function from pyomo
from pyomo.environ import log10
from pyomo.environ import units as pyunits

# ----------------------------------------------------------------------------
class log_power_law():

    @staticmethod
    def build_parameters(rblock, config):
        pass

    @staticmethod
    def return_expression(b, rblock, r_idx, T):
        e = 0
        # Get reaction orders and construct log power law expression
        for p, j in b.state_ref.params._phase_component_set:
            o = rblock.reaction_order[p, j]
            e += o*log10(get_concentration_term(b, r_idx)[p, j] * (1.0/pyunits.get_units(get_concentration_term(b, r_idx)[p, j])))
        if pyunits.get_units(b.k_eq[r_idx]) is not None:
            return log10(b.k_eq[r_idx] * (1/pyunits.get_units(b.k_eq[r_idx])) ) == e
        else:
            return log10(b.k_eq[r_idx]) == e
