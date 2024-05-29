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

from copy import deepcopy

# Import Pyomo libraries
from pyomo.environ import (
    Set, 
    Var,
    Param,
    value,
    log, 
    units as pyunits,
)
from pyomo.common.config import ConfigValue, In
# Import IDAES cores
from idaes.core import declare_process_block_class

import idaes.core.util.scaling as iscale

from idaes.core.util.exceptions import InitializationError, ConfigurationError
from watertap.unit_models.ion_exchange.ion_exchange_base import (
    IonExchangeBaseData,
    IonExchangeType,
)

__author__ = "Kurban Sitterley"


@declare_process_block_class("IonExchangeMultiComp")
class IonExchangeMultiCompData(IonExchangeBaseData):
    """
    Ion exchange multi component Clark model
    """

    CONFIG = IonExchangeBaseData.CONFIG()

    CONFIG.declare(
        "reactive_ions",
        ConfigValue(
            default=list(),
            domain=list,
            description="Designates other reactive species",
        ),
    )

    CONFIG.declare(
        "number_traps",
        ConfigValue(
            default=5,
            domain=int,
            description="Number of trapezoids to use for steady-state effluent concentration estimation",
        ),
    )

    def build(self):
        super().build()

        prop_in = self.process_flow.properties_in[0]
        regen = self.regeneration_stream[0]
        comps = self.config.property_package.component_list
        target_component = self.config.target_component

        target_component = self.config.target_component
        reactive_ions = self.config.reactive_ions

        assert target_component in reactive_ions
        self.target_component_set = Set(
            initialize=[target_component]
        )  
        inerts = comps - self.target_component_set
        self.reactive_ion_set = Set(initialize=reactive_ions)
        self.inert_set = Set(initialize=comps - self.reactive_ion_set)

        if len(self.target_component_set) > 1:
            raise ConfigurationError(
                f"IonExchange0D can only accept a single target ion but {len(self.target_component_set)} were provided."
            )
        if self.config.property_package.charge_comp[target_component].value > 0:
            self.ion_exchange_type = IonExchangeType.cation
        elif self.config.property_package.charge_comp[target_component].value < 0:
            self.ion_exchange_type = IonExchangeType.anion
        else:
            raise ConfigurationError("Target ion must have non-zero charge.")
        
        for j in inerts:
            self.process_flow.mass_transfer_term[:, "Liq", j].fix(0)
            regen.get_material_flow_terms("Liq", j).fix(0)

        self.num_traps = self.config.number_traps
        self.trap_disc = range(self.num_traps + 1)
        self.trap_index = self.trap_disc[1:]

        self.eps = Param(initialize=1e-4, mutable=True)

        self.c_trap_min = Param(
            self.reactive_ion_set,
            initialize=0.01,
            mutable=True,
            doc="Minimum relative breakthrough concentration for estimating area under curve",
        )

        self.N_Sc = Var(
            self.reactive_ion_set,
            initialize=700,
            units=pyunits.dimensionless,
            doc="Schmidt number",
        )

        self.N_Sh = Var(
            self.reactive_ion_set,
            initialize=30,
            units=pyunits.dimensionless,
            doc="Sherwood number",
        )

        self.c_norm = Var(
            self.reactive_ion_set,
            initialize=0.5,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Dimensionless (relative) concentration [Ct/C0] of target ion",
        )

        self.c_traps = Var(
            self.reactive_ion_set,
            self.trap_disc,
            initialize=0.5,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Normalized breakthrough concentrations for estimating area under breakthrough curve",
        )

        self.tb_traps = Var(
            self.reactive_ion_set,
            self.trap_disc,
            initialize=1e6,
            bounds=(0, None),
            units=pyunits.second,
            doc="Breakthrough times for estimating area under breakthrough curve",
        )

        self.traps = Var(
            self.reactive_ion_set,
            self.trap_index,
            initialize=0.01,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Trapezoid areas for estimating area under breakthrough curve",
        )

        for c in self.reactive_ion_set:
            self.c_traps[(c, 0)].fix(0)
            self.tb_traps[(c, 0)].fix(0)

        self.c_norm_avg = Var(
            self.reactive_ion_set,
            initialize=0.25,
            bounds=(0, 2),
            units=pyunits.dimensionless,
            doc="Sum of trapezoid areas",
        )

        self.freundlich_n = Var(
            self.reactive_ion_set,
            initialize=1.5,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Freundlich isotherm exponent",
        )

        self.mass_transfer_coeff = Var(  # k_T
            self.reactive_ion_set,
            initialize=0.001,
            units=pyunits.s**-1,
            bounds=(0, None),
            doc="Mass transfer coefficient for Clark model (kT)",
        )

        self.bv_50 = Var(  # BV_50
            self.reactive_ion_set,
            initialize=2e5,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Bed volumes of feed at 50 percent breakthrough",
        )

        @self.Expression(self.reactive_ion_set, doc="Breakthrough concentration")
        def c_breakthru(b, j):
            return b.c_norm[j] * prop_in.conc_mass_phase_comp["Liq", j]
        
        @self.Constraint(self.reactive_ion_set, doc="Schmidt number")
        def eq_Sc(b, j):  # Eq. 3.359, Inglezakis + Poulopoulos
            return (
                b.N_Sc[j]
                == prop_in.visc_k_phase["Liq"] / prop_in.diffus_phase_comp["Liq", j]
            )

        @self.Constraint(self.reactive_ion_set, doc="Sherwood number")
        def eq_Sh(b, j):  # Eq. 3.346, Inglezakis + Poulopoulos
            return (
                b.N_Sh[j]
                == b.Sh_A
                * b.bed_porosity**b.Sh_exp_A
                * b.N_Re**b.Sh_exp_B
                * b.N_Sc[j] ** b.Sh_exp_C
            )

        @self.Constraint(
            self.reactive_ion_set, doc="Mass transfer for regeneration stream"
        )
        def eq_mass_transfer_regen(b, j):
            return (
                regen.get_material_flow_terms("Liq", j)
                == -b.process_flow.mass_transfer_term[0, "Liq", j]
            )

        @self.Constraint(doc="Bed volumes at breakthrough")
        def eq_bv(b):
            return b.breakthrough_time * b.loading_rate == b.bv * b.bed_depth

        @self.Constraint(
            self.reactive_ion_set, doc="Clark equation with fundamental constants"
        )  # Croll et al (2023), Eq.9
        def eq_clark(b, j):
            left_side = (
                (b.mass_transfer_coeff[j] * b.bed_depth * (b.freundlich_n[j] - 1))
                / (b.bv_50[j] * b.vel_bed)
            ) * (b.bv_50[j] - b.bv)

            right_side = log(
                ((1 / b.c_norm[j]) ** (b.freundlich_n[j] - 1) - 1)
                / (2 ** (b.freundlich_n[j] - 1) - 1)
            )
            return left_side - right_side == 0

        @self.Constraint(
            self.reactive_ion_set,
            self.trap_index,
            doc="Evenly spaced c_norm for trapezoids",
        )
        def eq_c_traps(b, j, k):
            if k == max(b.traps_index):
                return b.c_traps[j, k] == b.c_norm[j]
            else:
                return b.c_traps[j, k] == b.c_trap_min[j] + (b.trap_disc[k] - 1) * (
                    (b.c_norm[j] - b.c_trap_min[j]) / (b.num_traps - 1)
                )
            
        @self.Expression(self.reactive_ion_set, self.trap_disc, doc="BV for trapezoids")
        def bv_traps(b, j, k):
            if k == 0:
                return 0
            elif k == max(self.trap_index):
                return b.bv
            else:
                return (b.tb_traps[j, k] * b.vel_bed) / b.bed_depth
            
        @self.Constraint(
            self.reactive_ion_set,
            self.trap_index,
            doc="Breakthru time calc for trapezoids",
        )
        def eq_tb_traps(b, j, k):
            # if (j, k) == (target_ion, max(self.trap_index)):
            if k == max(self.trap_index):
                return b.tb_traps[j, k] == b.t_breakthru
            else:
                bv_traps = (b.tb_traps[j, k] * b.vel_bed) / b.bed_depth
                left_side = (
                    (b.mass_transfer_coeff[j] * b.bed_depth * (b.freundlich_n[j] - 1))
                    / (b.bv_50[j] * b.vel_bed)
                ) * (b.bv_50[j] - bv_traps)

                right_side = log(
                    ((1 / b.c_traps[j, k]) ** (b.freundlich_n[j] - 1) - 1)
                    / (2 ** (b.freundlich_n[j] - 1) - 1)
                )
                return left_side - right_side == 0

        @self.Constraint(
            self.reactive_ion_set, self.trap_index, doc="Area of trapezoids"
        )
        def eq_traps(b, j, k):
            return b.traps[j, k] == (
                b.tb_traps[j, k] - b.tb_traps[j, k - 1]
            ) / b.tb_traps[j, self.num_traps] * (
                (b.c_traps[j, k] + b.c_traps[j, k - 1]) / 2
            )

        @self.Constraint(
            self.reactive_ion_set, doc="Average relative effluent concentration"
        )
        def eq_c_norm_avg(b, j):
            return b.c_norm_avg[j] == sum(b.traps[j, k] for k in b.trap_index)

        @self.Constraint(
            self.target_component_set,
            doc="CV mass transfer term",
        )
        @self.Constraint(
            self.reactive_ion_set,
            doc="CV mass transfer term",
        )
        def eq_mass_transfer_reactive_ions(b, j):
            return (1 - b.c_norm_avg[j]) * prop_in.get_material_flow_terms(
                "Liq", j
            ) == -b.process_flow.mass_transfer_term[0, "Liq", j]


    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if iscale.get_scaling_factor(self.freundlich_n) is None:
            iscale.set_scaling_factor(self.freundlich_n, 0.1)

        if iscale.get_scaling_factor(self.mass_transfer_coeff) is None:
            iscale.set_scaling_factor(self.mass_transfer_coeff, 10)

        if iscale.get_scaling_factor(self.bv_50) is None:
            iscale.set_scaling_factor(self.bv_50, 1e-5)

        if iscale.get_scaling_factor(self.tb_traps) is None:
            sf = iscale.get_scaling_factor(self.breakthrough_time)
            iscale.set_scaling_factor(self.tb_traps, sf)

        if iscale.get_scaling_factor(self.c_traps) is None:
            iscale.set_scaling_factor(self.c_traps, 1)

        if iscale.get_scaling_factor(self.traps) is None:
            iscale.set_scaling_factor(self.traps, 1e3)

        if iscale.get_scaling_factor(self.c_norm_avg) is None:
            iscale.set_scaling_factor(self.c_norm_avg, 1e2)

        for ind, c in self.eq_clark.items():
            if iscale.get_scaling_factor(c) is None:
                iscale.constraint_scaling_transform(c, 1e-2)

        for ind, c in self.eq_traps.items():
            if iscale.get_scaling_factor(c) is None:
                iscale.constraint_scaling_transform(c, 1e2)

        # elif self.config.isotherm == IsothermType.freundlich:
        #     var_dict[f"BV at Breakthrough"] = self.bv
        #     var_dict[f"BV at 50% Breakthrough"] = self.bv_50
        #     var_dict[f"Freundlich n"] = self.freundlich_n
        #     var_dict[f"Clark Mass Transfer Coeff."] = self.mass_transfer_coeff