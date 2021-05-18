'''
    This file creates alternative forms of the van't Hoff equation that
    uses enthalpy and entropy to calculate the equilibrium coefficient, rather
    than a reference state coefficient and temperature. This was necessary as
    a majority of the equilibrium data in literature for aqueous systems uses
    this form of the van't Hoff expression.

    K = exp(-dH/R/T + dS/R)
'''

from pyomo.environ import exp, Var, Param, units as pyunits
from pyomo.environ import value

from idaes.generic_models.properties.core.generic.generic_reaction import \
    ConcentrationForm
from idaes.core.util.misc import set_param_from_config
from idaes.core.util.constants import Constants as c
from idaes.core.util.exceptions import BurntToast, ConfigurationError


# -----------------------------------------------------------------------------
# Function for alternative van't Hoff expression using enthalpy and entropy
class van_t_hoff_aqueous():

    @staticmethod
    def build_parameters(rblock, config):
        parent = rblock.parent_block()
        units = parent.get_metadata().derived_units

        c_form = config.concentration_form
        if c_form is None:
            raise ConfigurationError(
                "{} concentration_form configuration argument was not set. "
                "Please ensure that this argument is included in your "
                "configuration dict.".format(rblock.name))
        elif (c_form == ConcentrationForm.moleFraction or
              c_form == ConcentrationForm.massFraction):
            raise ConfigurationError(
                "{} concentration_form configuration argument must be either "
                " set as ConcentrationForm.molarity or ConcentrationForm.activity.")
        else:
            order = 0
            for p, j in parent.config.property_package._phase_component_set:
                order += rblock.reaction_order[p, j].value

            if (c_form == ConcentrationForm.molarity or
                    c_form == ConcentrationForm.activity):
                c_units = units["density_mole"]
            elif c_form == ConcentrationForm.molality:
                raise ConfigurationError(
                    "{} concentration_form configuration argument must be either "
                    " set as ConcentrationForm.molarity or ConcentrationForm.activity.")
            elif c_form == ConcentrationForm.partialPressure:
                raise ConfigurationError(
                    "{} concentration_form configuration argument must be either "
                    " set as ConcentrationForm.molarity or ConcentrationForm.activity.")
            else:
                raise BurntToast(
                    "{} get_concentration_term received unrecognised "
                    "ConcentrationForm ({}). This should not happen - please "
                    "contact the IDAES developers with this bug."
                    .format(rblock.name, c_form))

            e_units = c_units**order

        # Create a variable for entropy of reaction
        #   NOTE: We do not need to create a variable for reaction enthalpy
        #       because that variable already exists from another method
        rblock.ds_rxn_ref = Var(doc="Specific entropy of reaction at reference state",
                                units=pyunits.J/pyunits.mol/pyunits.K)
        set_param_from_config(rblock, param="ds_rxn_ref", config=config)

        # Add a param for reference temperature (not read in from file)
        rblock.T_eq_ref = Param(
                doc="Reference temperature for equilibrium constant",
                units=units["temperature"],
                initialize=293)

        # Convert the enthalpy value (if necessary)
        dh_val = pyunits.convert_value(rblock.dh_rxn_ref.value, from_units=pyunits.get_units(rblock.dh_rxn_ref), to_units=pyunits.J/pyunits.mol)

        # Calculate a k_eq at reference temperature in K
        k_eq = exp(-(dh_val/8.3145/rblock.T_eq_ref.value) + (rblock.ds_rxn_ref.value/8.3145))

        # Add parameter for equilibrium coefficient at reference temperature
        #   and initialize its value to the above calculated value, which is
        #   done on a (mol/L) basis for aqueous chemistry
        rblock.k_eq_ref = Param(
                doc="Equilibrium constant at reference state",
                units=e_units,
                initialize=k_eq)

        # Convert units of k value to what IDAES wants
        rblock.k_eq_ref.value = pyunits.convert_value(rblock.k_eq_ref.value, from_units=(pyunits.mol/pyunits.L)**order, to_units=e_units)

        # Lastly, convert units for reference temperature back to the expected Units
        rblock.T_eq_ref.value = pyunits.convert_value(rblock.T_eq_ref.value, from_units=pyunits.K, to_units=units["temperature"])

    @staticmethod
    def return_expression(b, rblock, r_idx, T):
        units = rblock.parent_block().get_metadata().derived_units

        return rblock.k_eq_ref * exp(
            -(b.dh_rxn[r_idx] /
              pyunits.convert(c.gas_constant,
                              to_units=units["gas_constant"])) *
            (1/T - 1/rblock.T_eq_ref))
