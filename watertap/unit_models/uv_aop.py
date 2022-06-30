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

# Import Pyomo libraries
from pyomo.environ import (
    Block,
    Set,
    Var,
    Param,
    Suffix,
    NonNegativeReals,
    Reference,
    units as pyunits,
)
from pyomo.environ import *
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util import get_solver
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.util.math import smooth_min, smooth_max

_log = idaeslog.getLogger(__name__)


@declare_process_block_class("Ultraviolet0D")
class Ultraviolet0DData(UnitModelBlockData):
    """
    Standard UV Unit Model Class:
    - zero dimensional model
    - steady state only
    - single liquid phase only
    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False. UV units do not support dynamic
    behavior.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False. UV units do not have defined volume, thus
    this must be False.""",
        ),
    )
    CONFIG.declare(
        "material_balance_type",
        ConfigValue(
            default=MaterialBalanceType.useDefault,
            domain=In(MaterialBalanceType),
            description="Material balance construction flag",
            doc="""Indicates what type of mass balance should be constructed,
    **default** - MaterialBalanceType.useDefault.
    **Valid values:** {
    **MaterialBalanceType.useDefault - refer to property package for default
    balance type
    **MaterialBalanceType.none** - exclude material balances,
    **MaterialBalanceType.componentPhase** - use phase component balances,
    **MaterialBalanceType.componentTotal** - use total component balances,
    **MaterialBalanceType.elementTotal** - use total element balances,
    **MaterialBalanceType.total** - use total material balance.}""",
        ),
    )
    CONFIG.declare(
        "energy_balance_type",
        ConfigValue(
            default=EnergyBalanceType.useDefault,
            domain=In(EnergyBalanceType),
            description="Energy balance construction flag",
            doc="""Indicates what type of energy balance should be constructed,
    **default** - EnergyBalanceType.useDefault.
    **Valid values:** {
    **EnergyBalanceType.useDefault - refer to property package for default
    balance type
    **EnergyBalanceType.none** - exclude energy balances,
    **EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
    **EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
    **EnergyBalanceType.energyTotal** - single energy balance for material,
    **EnergyBalanceType.energyPhase** - energy balances for each phase.}""",
        ),
    )
    CONFIG.declare(
        "momentum_balance_type",
        ConfigValue(
            default=MomentumBalanceType.pressureTotal,
            domain=In(MomentumBalanceType),
            description="Momentum balance construction flag",
            doc="""Indicates what type of momentum balance should be constructed,
    **default** - MomentumBalanceType.pressureTotal.
    **Valid values:** {
    **MomentumBalanceType.none** - exclude momentum balances,
    **MomentumBalanceType.pressureTotal** - single pressure balance for material,
    **MomentumBalanceType.pressurePhase** - pressure balances for each phase,
    **MomentumBalanceType.momentumTotal** - single momentum balance for material,
    **MomentumBalanceType.momentumPhase** - momentum balances for each phase.}""",
        ),
    )
    CONFIG.declare(
        "has_pressure_change",
        ConfigValue(
            default=False,
            domain=In([True, False]),
            description="Pressure change term construction flag",
            doc="""Indicates whether terms for pressure change should be
    constructed,
    **default** - False.
    **Valid values:** {
    **True** - include pressure change terms,
    **False** - exclude pressure change terms.}""",
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
    **default** - useDefault.
    **Valid values:** {
    **useDefault** - use default package from parent model or flowsheet,
    **PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
    and used when constructing these,
    **default** - None.
    **Valid values:** {
    see property package for documentation.}""",
        ),
    )

    def build(self):
        # Call UnitModel.build to setup dynamics
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        if list(self.config.property_package.phase_list) != ["Liq"]:
            raise ConfigurationError(
                "UV model only supports one liquid phase ['Liq'],"
                "the property package has specified the following phases {}".format(
                    [p for p in self.config.property_package.phase_list]
                )
            )

        units_meta = self.config.property_package.get_metadata().get_derived_units

        # TODO: update IDAES such that solvent and solute lists are automatically created on the parameter block
        self.solvent_list = Set()
        self.solute_list = Set()
        for c in self.config.property_package.component_list:
            comp = self.config.property_package.get_component(c)
            try:
                if comp.is_solvent():
                    self.solvent_list.add(c)
                if comp.is_solute():
                    self.solute_list.add(c)
            except TypeError:
                raise ConfigurationError(
                    "UV model only supports one solvent and one or more solutes,"
                    "the provided property package has specified a component '{}' "
                    "that is not a solvent or solute".format(c)
                )
        if len(self.solvent_list) > 1:
            raise ConfigurationError(
                "UV model only supports one solvent component,"
                "the provided property package has specified {} solvent components".format(
                    len(self.solvent_list)
                )
            )

        # Add unit parameters
        self.inactivation_rate = Var(
            self.config.property_package.phase_list,
            self.config.property_package.solute_set,
            initialize=2.5e-4,
            bounds=(1e-18, 100),
            domain=NonNegativeReals,
            units=units_meta("time") ** 2 * units_meta("mass") ** -1,
            doc="Inactivation rate coefficient with respect to uv dose.",
        )

        self.rate_constant = Var(
            self.config.property_package.phase_list,
            self.config.property_package.solute_set,
            initialize=2.5e-3,
            bounds=(1e-18, 100),
            domain=NonNegativeReals,
            units=units_meta("time") ** -1,
            doc="Overall pseudo-first order rate constant.",
        )

        self.photolysis_rate_constant = Var(
            self.config.property_package.phase_list,
            self.config.property_package.solute_set,
            initialize=2e-3,
            bounds=(1e-18, 100),
            domain=NonNegativeReals,
            units=units_meta("time") ** -1,
            doc="Pseudo-first order rate constant for direct photolysis of component.",
        )

        self.reaction_rate_constant = Var(
            self.config.property_package.phase_list,
            self.config.property_package.solute_set,
            initialize=5e-4,
            bounds=(1e-18, 100),
            domain=NonNegativeReals,
            units=units_meta("time") ** -1,
            doc="Pseudo-first order rate constant for indirect photolysis of component.",
        )

        self.dens_solvent = Param(
            initialize=1000,
            units=units_meta("mass") * units_meta("length") ** -3,
            doc="Pure water density",
        )

        # Add uv variables
        self.uv_dose = Var(
            initialize=5000,
            bounds=(1e-18, 10000),
            domain=NonNegativeReals,
            units=units_meta("mass") * units_meta("time") ** -2,
            doc="UV dose.",
        )
        self.uv_intensity = Var(
            initialize=10,
            bounds=(1e-18, 10000),
            domain=NonNegativeReals,
            units=units_meta("mass") * units_meta("time") ** -3,
            doc="Average intensity of UV light.",
        )
        self.exposure_time = Var(
            initialize=500,
            bounds=(1e-18, 10000),
            domain=NonNegativeReals,
            units=units_meta("time"),
            doc="Exposure time of UV light.",
        )
        self.UVT = Var(
            initialize=0.9,
            bounds=(0, 1),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="UV transmittance.",
        )
        self.UVA = Var(
            initialize=0.045,
            bounds=(0, 10),
            domain=NonNegativeReals,
            units=pyunits.cm**-1,
            doc="UV absorbance.",
        )

        # Add electricity parameters
        self.electricity_demand_phase_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.solute_set,
            initialize=1,
            bounds=(0, None),
            units=units_meta("mass")
            * units_meta("length") ** 2
            * units_meta("time") ** -3,
            doc="Electricity demand per component",
        )

        # self.electricity_demand_minimum = Var(
        #     self.flowsheet().config.time,
        #     initialize=1,
        #     bounds=(0, None),
        #     units=units_meta("mass")
        #     * units_meta("length") ** 2
        #     * units_meta("time") ** -3,
        #     doc="Minimum electricity demand of unit",
        # )

        self.electrical_efficiency_phase_comp = Var(
            self.flowsheet().time,
            self.config.property_package.phase_list,
            self.config.property_package.solute_set,
            initialize=1,
            bounds=(0, None),
            units=units_meta("mass")
            * units_meta("length") ** -1
            * units_meta("time") ** -2,
            doc="Electricity efficiency per log order reduction (EE/O)",
        )

        self.lamp_efficiency = Var(
            initialize=0.3,
            bounds=(0, 1),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Lamp efficiency",
        )

        # Build control volume for UV unit
        self.control_volume = ControlVolume0DBlock(
            default={
                "dynamic": False,
                "has_holdup": False,
                "property_package": self.config.property_package,
                "property_package_args": self.config.property_package_args,
            }
        )

        self.control_volume.add_state_blocks(has_phase_equilibrium=False)

        self.control_volume.add_material_balances(
            balance_type=self.config.material_balance_type, has_mass_transfer=True
        )

        self.control_volume.add_energy_balances(
            balance_type=self.config.energy_balance_type, has_enthalpy_transfer=True
        )

        self.control_volume.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change,
        )

        # Add Ports
        self.add_inlet_port(name="inlet", block=self.control_volume)
        self.add_outlet_port(name="outlet", block=self.control_volume)

        # References for control volume
        # pressure change
        if (
            self.config.has_pressure_change is True
            and self.config.momentum_balance_type != "none"
        ):
            self.deltaP = Reference(self.control_volume.deltaP)

        # UV dose
        @self.Constraint(
            doc="Constraint for UV dose",
        )
        def eq_uv_dose(b):
            return b.uv_dose == b.uv_intensity * b.exposure_time

        @self.Constraint(
            doc="Constraint for UV absorbance",
        )
        def eq_UVA(b):
            return b.UVA == -log10(b.UVT) / (1 * pyunits.cm)

        # rate constant
        @self.Constraint(
            self.config.property_package.phase_list,
            self.config.property_package.solute_set,
            doc="Constraint for pseudo-first order rate constant with respect to uv intensity",
        )
        def eq_rate_constant(b, p, j):
            return b.rate_constant[p, j] == b.uv_intensity * b.inactivation_rate[p, j]

        @self.Constraint(
            self.config.property_package.phase_list,
            self.config.property_package.solute_set,
            doc="Constraint for pseudo-first order rate constant with respect to direct and indirect photolysis",
        )
        def eq_overall_rate_constant(b, p, j):
            return (
                b.rate_constant[p, j]
                == b.photolysis_rate_constant[p, j] + b.reaction_rate_constant[p, j]
            )

        # mass transfer
        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Constraints for solvent and solute concentration in outlet stream.",
        )
        def eq_outlet_conc(b, t, p, j):
            prop_in = b.control_volume.properties_in[t]
            prop_out = b.control_volume.properties_out[t]
            comp = self.config.property_package.get_component(j)
            if comp.is_solvent():
                return prop_out.get_material_flow_terms(
                    p, j
                ) == prop_in.get_material_flow_terms(p, j)
            elif comp.is_solute():
                return prop_out.get_material_flow_terms(
                    p, j
                ) == prop_in.get_material_flow_terms(p, j) * exp(
                    pyunits.convert(
                        -b.uv_dose * b.inactivation_rate[p, j],
                        to_units=pyunits.dimensionless,
                    )
                )

        # electricity
        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.solute_set,
            doc="Constraints for electricity demand of the UV reactor.",
        )
        def eq_electricity_demand_phase_comp(b, t, p, j):
            prop_in = b.control_volume.properties_in[t]
            return b.electricity_demand_phase_comp[t, p, j] == (
                b.electrical_efficiency_phase_comp[t, p, j]
                * prop_in.flow_vol
                * log10(
                    1
                    / exp(
                        pyunits.convert(
                            -b.uv_dose * b.inactivation_rate[p, j],
                            to_units=pyunits.dimensionless,
                        )
                    )
                )
                / b.lamp_efficiency
            )

        # TODO: add minimum electricity demand for multiple solutes
        # @self.Constraint(
        #     self.flowsheet().config.time,
        #     doc="Constraints for minimum electricity demand of the UV reactor.",
        # )
        # def eq_minimum_electricity_demand(b, t):
        #     return b.electricity_demand_minimum[t] == smooth_max(b.electricity_demand_phase_comp[t, p, j])

    def initialize_build(
        blk, state_args=None, outlvl=idaeslog.NOTSET, solver=None, optarg=None
    ):
        """
        General wrapper for pressure changer initialization routines

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None)
            solver : str indicating which solver to use during
                     initialization (default = None)

        Returns: None
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")
        # Set solver options
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Initialize holdup block
        flags = blk.control_volume.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )
        init_log.info_high("Initialization Step 1 Complete.")
        # ---------------------------------------------------------------------
        # Initialize permeate
        # Set state_args from inlet state
        if state_args is None:
            state_args = {}
            state_dict = blk.control_volume.properties_in[
                blk.flowsheet().config.time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        init_log.info_high("Initialization Step 2 Complete.")

        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.control_volume.release_state(flags, outlvl + 1)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

    def _get_performance_contents(self, time_point=0):
        # TODO: add other performance constants
        var_dict = {}
        if hasattr(self, "deltaP"):
            var_dict["Pressure Change"] = self.deltaP[time_point]

        return {"vars": var_dict}

    # TODO: add costing
    # def get_costing(self, module=None, **kwargs):
    #     self.costing = Block()
    #     module.UV_costing(self.costing, **kwargs)

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # TODO: require users to set scaling factor for uv dose, uv_intensity and exposure time
        # setting scaling factors for variables
        # these variables should have user input, if not there will be a warning
        if iscale.get_scaling_factor(self.uv_intensity) is None:
            sf = iscale.get_scaling_factor(self.uv_intensity, default=0.1, warning=True)
        iscale.set_scaling_factor(self.uv_intensity, sf)

        if iscale.get_scaling_factor(self.exposure_time) is None:
            sf = iscale.get_scaling_factor(
                self.exposure_time, default=1e-2, warning=True
            )
            iscale.set_scaling_factor(self.exposure_time, sf)

        if iscale.get_scaling_factor(self.uv_dose) is None:
            sf = iscale.get_scaling_factor(self.uv_dose, default=1e-3, warning=True)
        iscale.set_scaling_factor(self.uv_dose, sf)

        if iscale.get_scaling_factor(self.inactivation_rate) is None:
            sf = iscale.get_scaling_factor(
                self.inactivation_rate, default=1e4, warning=True
            )
        iscale.set_scaling_factor(self.inactivation_rate, sf)

        if iscale.get_scaling_factor(self.rate_constant) is None:
            sf = iscale.get_scaling_factor(
                self.rate_constant, default=1e3, warning=True
            )
        iscale.set_scaling_factor(self.rate_constant, sf)

        if iscale.get_scaling_factor(self.photolysis_rate_constant) is None:
            sf = iscale.get_scaling_factor(
                self.photolysis_rate_constant, default=1e3, warning=True
            )
        iscale.set_scaling_factor(self.photolysis_rate_constant, sf)

        if iscale.get_scaling_factor(self.reaction_rate_constant) is None:
            sf = iscale.get_scaling_factor(
                self.reaction_rate_constant, default=1e3, warning=True
            )
        iscale.set_scaling_factor(self.reaction_rate_constant, sf)

        if iscale.get_scaling_factor(self.electrical_efficiency_phase_comp) is None:
            sf = iscale.get_scaling_factor(
                self.electrical_efficiency_phase_comp, default=1e-6, warning=True
            )
        iscale.set_scaling_factor(self.electrical_efficiency_phase_comp, sf)

        # these variables do not typically require user input,
        # will not override if the user does provide the scaling factor
        if iscale.get_scaling_factor(self.lamp_efficiency) is None:
            iscale.set_scaling_factor(self.lamp_efficiency, 10)

        if iscale.get_scaling_factor(self.UVT) is None:
            iscale.set_scaling_factor(self.UVT, 10)

        if iscale.get_scaling_factor(self.UVA) is None:
            iscale.set_scaling_factor(self.UVA, 100)

        if iscale.get_scaling_factor(self.dens_solvent) is None:
            sf = iscale.get_scaling_factor(
                self.control_volume.properties_in[0].dens_mass_phase["Liq"]
            )
            iscale.set_scaling_factor(self.dens_solvent, sf)

        for (t, p, j), v in self.electricity_demand_phase_comp.items():
            if iscale.get_scaling_factor(v) is None:
                removal = -iscale.get_scaling_factor(
                    self.uv_dose
                ) * iscale.get_scaling_factor(self.inactivation_rate[p, j])
                sf = (
                    iscale.get_scaling_factor(
                        self.electrical_efficiency_phase_comp[t, p, j]
                    )
                    * (1 / log10(1 / exp(removal)))
                    * iscale.get_scaling_factor(
                        self.control_volume.properties_in[t].flow_vol
                    )
                    / iscale.get_scaling_factor(self.lamp_efficiency)
                )
                iscale.set_scaling_factor(v, sf)

        for (t, p, j), v in self.control_volume.mass_transfer_term.items():
            if iscale.get_scaling_factor(v) is None:
                sf = iscale.get_scaling_factor(
                    self.control_volume.properties_in[t].get_material_flow_terms(p, j)
                )
                comp = self.config.property_package.get_component(j)
                if comp.is_solute:
                    sf *= 1e2  # solute typically has mass transfer 2 orders magnitude less than flow
                iscale.set_scaling_factor(v, sf)

        # TODO: update IDAES control volume to scale mass_transfer and enthalpy_transfer
        for ind, v in self.control_volume.mass_transfer_term.items():
            (t, p, j) = ind
            if iscale.get_scaling_factor(v) is None:
                sf = iscale.get_scaling_factor(
                    self.control_volume.mass_transfer_term[t, p, j]
                )
                iscale.constraint_scaling_transform(
                    self.control_volume.material_balances[t, j], sf
                )

        for t, v in self.control_volume.enthalpy_transfer.items():
            if iscale.get_scaling_factor(v) is None:
                sf = iscale.get_scaling_factor(
                    self.control_volume.properties_in[t].enth_flow
                )
                iscale.set_scaling_factor(v, sf)
                iscale.constraint_scaling_transform(
                    self.control_volume.enthalpy_balances[t], sf
                )

        # transforming constraints
        for c in self.eq_uv_dose.values():
            if iscale.get_scaling_factor(self.uv_dose) is None:
                sf = iscale.get_scaling_factor(
                    self.uv_intensity
                ) * iscale.get_scaling_factor(self.exposure_time)
            else:
                sf = iscale.get_scaling_factor(self.uv_dose)
            iscale.constraint_scaling_transform(c, sf)

        for c in self.eq_UVA.values():
            if iscale.get_scaling_factor(self.UVA) is None:
                sf = -log10(iscale.get_scaling_factor(self.UVT))
            else:
                sf = iscale.get_scaling_factor(self.UVA)
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_rate_constant.items():
            if iscale.get_scaling_factor(self.rate_constant) is None:
                sf = iscale.get_scaling_factor(
                    self.uv_intensity
                ) * iscale.get_scaling_factor(self.inactivation_rate[ind])
            else:
                sf = iscale.get_scaling_factor(self.rate_constant[ind])
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_overall_rate_constant.items():
            if iscale.get_scaling_factor(self.rate_constant) is None:
                sf = iscale.get_scaling_factor(
                    self.photolysis_rate_constant[ind]
                ) + iscale.get_scaling_factor(self.reaction_rate_constant[ind])
            else:
                sf = iscale.get_scaling_factor(self.rate_constant[ind])
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_outlet_conc.items():
            (t, p, j) = ind
            sf = iscale.get_scaling_factor(
                self.control_volume.properties_in[t].get_material_flow_terms(p, j)
            )
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_electricity_demand_phase_comp.items():
            (t, p, j) = ind
            sf = iscale.get_scaling_factor(self.electricity_demand_phase_comp[t, p, j])
            iscale.constraint_scaling_transform(c, sf)
