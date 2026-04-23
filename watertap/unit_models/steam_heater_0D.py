#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################


from pyomo.environ import NonNegativeReals, Var, units as pyunits

from idaes.core import declare_process_block_class
from idaes.core.util.model_serializer import to_json, from_json, StoreSpec
from idaes.models.unit_models.heat_exchanger import HeatExchangerData
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.core.solvers import get_solver
from watertap.costing.unit_models.heat_exchanger import cost_heat_exchanger


_log = idaeslog.getLogger(__name__)


__author__ = "Elmira Shamlou"


@declare_process_block_class("SteamHeater0D")
class SteamHeater0DData(HeatExchangerData):
    CONFIG = HeatExchangerData.CONFIG()

    def build(self):
        super().build()

        self.pressure_deltaP = Var(
            self.flowsheet().time,
            initialize=0.0,
            bounds=(0.0, None),
            domain=NonNegativeReals,
            units=pyunits.Pa,
        )
        self.pressure_deltaP.fix(0)

        @self.Constraint(self.flowsheet().time)
        def outlet_pressure_sat(b, t):
            return (
                b.hot_side.properties_out[t].pressure
                == b.hot_side.properties_out[t].pressure_sat
                + b.pressure_deltaP[t]
            )

        self._apply_total_condensation()

    def _apply_total_condensation(self):
        for t in self.flowsheet().time:
            for j in self.hot_side.config.property_package.component_list:
                v = self.hot_side.properties_out[t].flow_mass_phase_comp[
                    "Vap", j
                ]
                v.fix(v.lb)

    def set_subcooling_margin(self, margin_Pa):
        for t in self.flowsheet().time:
            self.pressure_deltaP[t].fix(float(margin_Pa))

    def release_subcooling_margin(self):
        for t in self.flowsheet().time:
            self.pressure_deltaP[t].unfix()

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        for t in self.flowsheet().time:
            sf_p = iscale.get_scaling_factor(
                self.hot_side.properties_out[t].pressure, default=1e-5
            )
            if iscale.get_scaling_factor(self.pressure_deltaP[t]) is None:
                iscale.set_scaling_factor(self.pressure_deltaP[t], sf_p)
            if iscale.get_scaling_factor(self.outlet_pressure_sat[t]) is None:
                iscale.constraint_scaling_transform(
                    self.outlet_pressure_sat[t], sf_p
                )

    def initialize_build(self, *args, **kwargs):
        solver = kwargs.get("solver", None)
        optarg = kwargs.get("optarg", {})
        outlvl = kwargs.get("outlvl", idaeslog.NOTSET)
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        flags = to_json(
            self,
            return_dict=True,
            wts=StoreSpec.value_isfixed_isactive(only_fixed=True),
        )

        self.hot_side_outlet.unfix()
        self.cold_side_outlet.unfix()

        for t in self.flowsheet().time:
            for sb in (
                self.hot_side.properties_in[t],
                self.cold_side.properties_in[t],
            ):
                state_var_ids = set()
                for _name, var in sb.define_state_vars().items():
                    for idx in var:
                        state_var_ids.add(id(var[idx]))
                for v in sb.component_data_objects(Var, descend_into=False):
                    if id(v) not in state_var_ids and v.fixed:
                        v.unfix()

        super().initialize_build(*args, **kwargs)

        from_json(
            self,
            sd=flags,
            wts=StoreSpec.value_isfixed_isactive(only_fixed=True),
        )

        opt = get_solver(solver, optarg)
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info(
            "Initialization Complete: {}".format(idaeslog.condition(res))
        )

    @property
    def default_costing_method(self):
        return cost_heat_exchanger
