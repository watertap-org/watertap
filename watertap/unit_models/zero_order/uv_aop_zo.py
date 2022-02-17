###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################
"""
This module contains a zero-order representation of a UV-A0P unit.
operation.
"""

from pyomo.environ import Constraint, units as pyunits, Var
from pyomo.common.config import ConfigBlock, ConfigValue, In
from idaes.core.util.exceptions import ConfigurationError
from idaes.core import declare_process_block_class
from watertap.core import build_siso, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Adam Atia"


@declare_process_block_class("UVAOPZO")
class UVAOPZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a UV-AOP unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    CONFIG.declare("has_aop", ConfigValue(
            default=False,
            domain=In([True, False]),
            description="Oxidant construction flag",
            doc="""Indicates whether oxidant dose and flow terms should be constructed or not.
            **default** - False."""))

    def build(self):
        super().build()

        self._tech_type = "uv_aop"

        if self.config.process_subtype is None:
            raise ConfigurationError(
                f"{self.name} - zero-order uv/aop operations "
                "require the process_subtype configuration argument to be set")

        build_siso(self)
        constant_intensity(self)

        self.uv_reduced_equivalent_dose = Var(self.flowsheet().time,
                                              units=pyunits.mJ/pyunits.cm**2,
                                              doc="Reduced equivalent dosage")
        self.uvt_in = Var(self.flowsheet().time,
                          units=pyunits.dimensionless,
                          doc="UV transmittance of solution at UV reactor inlet")

        self._fixed_perf_vars.append(self.uv_reduced_equivalent_dose)
        self._fixed_perf_vars.append(self.uvt_in)


        if self.config.has_aop:
            self.oxidant_dose = Var(self.flowsheet().time,
                                    units=pyunits.mg/pyunits.L,
                                    doc="Oxidant dosage")

            self.chemical_flow_mass = Var(
                self.flowsheet().time,
                units=pyunits.kg/pyunits.s,
                bounds=(0, None),
                doc="Mass flow rate of oxidant solution")


            self._fixed_perf_vars.append(self.oxidant_dose)

            @self.Constraint(self.flowsheet().time,
                             doc="Chemical mass flow constraint")
            def chemical_flow_mass_constraint(b, t):
                return (b.chemical_flow_mass[t] ==
                        pyunits.convert(b.oxidant_dose[t]
                        * b.properties_in[t].flow_vol, to_units=pyunits.kg/pyunits.s))

            self._perf_var_dict["Oxidant Dosage (mg/L)"] = self.oxidant_dose
            self._perf_var_dict["Oxidant Flow (kg/s)"] = self.chemical_flow_mass

    def load_parameters_from_database(self, use_default_removal=False):
        """
                Method to load parameters from database, with modification for UV-AOP.

                Args:
                    use_default_removal - (optional) indicate whether to use defined
                                          default removal fraction if no specific value
                                          defined in database

                Returns:
                    None
                """
        # Get parameter dict from database
        if self._tech_type is None:
            raise NotImplementedError(
                f"{self.name} derived zero order unit model has not "
                f"implemented the _tech_type attribute. This is required "
                f"to identify the database file to load parameters from.")

        # Get parameter dict from database
        pdict = self.config.database.get_unit_operation_parameters(
            self._tech_type, subtype=self.config.process_subtype)

        if not self.config.has_aop:
            # make sure that a subtype WITH advanced oxidation is not selected when has_aop is False
            if 'oxidant_dose' in pdict:
                raise KeyError(f'{self.name} - database provided contains an entry'
                               f' for oxidant_dose, but has_aop was set to False.')

        if self._has_recovery_removal:
            self.set_recovery_and_removal(pdict, use_default_removal)

        for v in self._fixed_perf_vars:
            self.set_param_from_data(v, pdict)
            # set_param_from_data will handle exception when a subtype WITHOUT advanced oxidation is
            # selected when has_aop is True




