#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
Initial property package for seawater system
"""

# Import Pyomo libraries
from pyomo.environ import (
    check_optimal_termination,
    Constraint,
    Suffix,
    value,
    Var,
)
from pyomo.environ import units as pyunits
from pyomo.contrib.solver.common.util import NoSolutionError

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    StateBlock,
)
from idaes.core.base.components import Solute
from idaes.core.initialization.block_triangularization import (
    BlockTriangularizationInitializer,
)
from idaes.core.util.initialization import (
    fix_state_vars,
    revert_state_vars,
    solve_indexed_blocks,
)
import idaes.logger as idaeslog
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
)
from idaes.core.util.exceptions import (
    ConfigurationError,
    InitializationError,
    PropertyPackageError,
)
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import BurntToast

# Import WaterTAP libraries
from watertap.core.solvers import get_solver
from watertap.core.util.scaling import transform_property_constraints

from watertap.property_models.seawater_prop_pack import (
    SeawaterParameterData,
    SeawaterStateBlockData,
)
from watertap.property_models.bt_initialization_mixin import (
    BlockTriangularizationInitializationMixin,
)

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("SeawaterFTPwParameterBlock")
class SeawaterFTPwParameterData(SeawaterParameterData):
    """References
    Parameter block for a seawater property package. Two components: 'H2O' and 'TDS', and one phase: 'Liq'.

    This package was developed from the following references:
    - K.G.Nayar, M.H.Sharqawy, L.D.Banchik, and J.H.Lienhard V, "Thermophysical properties of seawater: A review and
    new correlations that include pressure dependence,"Desalination, Vol.390, pp.1 - 24, 2016.
    doi: 10.1016/j.desal.2016.02.024(preprint)
    - Mostafa H.Sharqawy, John H.Lienhard V, and Syed M.Zubair, "Thermophysical properties of seawater: A review of
    existing correlations and data,"Desalination and Water Treatment, Vol.16, pp.354 - 380, April 2010.
    (2017 corrections provided at http://web.mit.edu/seawater)
    Diffusivity for NaCl is being used temporarily based on
    Bartholomew & Mauter (2019) https://doi.org/10.1016/j.memsci.2018.11.067
    """

    CONFIG = SeawaterParameterData.CONFIG()

    def build(self):
        """
        Callable method for Block construction.
        """
        super().build()

    def _get_state_block_class(self):
        return SeawaterFTPwStateBlock

    def _set_default_scaling(self):
        self.set_default_scaling("temperature", 1e-2)
        self.set_default_scaling("pressure", 1e-6)
        self.set_default_scaling("mass_frac_phase_comp", 1, index=("Liq", "H2O"))
        self.set_default_scaling("mass_frac_phase_comp", 10, index=("Liq", "TDS"))
        self.set_default_scaling("flow_mass_phase", 1, index="Liq")
        self.set_default_scaling("dens_mass_phase", 1e-3, index="Liq")
        self.set_default_scaling("dens_mass_solvent", 1e-3)
        self.set_default_scaling("visc_d_phase", 1e3, index="Liq")
        self.set_default_scaling("osm_coeff", 1e0)
        self.set_default_scaling("enth_mass_phase", 1e-5, index="Liq")
        self.set_default_scaling("pressure_sat", 1e-5)
        self.set_default_scaling("cp_mass_phase", 1e-3, index="Liq")
        self.set_default_scaling("therm_cond_phase", 1e0, index="Liq")
        self.set_default_scaling("dh_vap_mass", 1e-6)
        self.set_default_scaling("diffus_phase_comp", 1e9)
        self.set_default_scaling("boiling_point_elevation_phase", 1e0, index="Liq")

    @classmethod
    def _trim_state_variables(cls, property_dict):
        property_dict["flow_mass_phase"]["method"] = None
        property_dict["mass_frac_phase_comp"]["method"] = None
        property_dict["temperature"]["method"] = None
        property_dict["pressure"]["method"] = None

    @classmethod
    def define_metadata(cls, obj):
        """Inherit property metadata from parent"""
        super().define_metadata(obj)


# We need to inherit from the mixin first so that its initialization methods
# overwrite the default StateBlock ones.
class _SeawaterFTPwStateBlock(
    BlockTriangularizationInitializationMixin,
    StateBlock,
):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """

    def fix_initialization_states(self, state_args=None):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """
        # Fix state variables
        flags = fix_state_vars(self, state_args)

        # Constraint on water concentration at outlet - unfix in these cases
        for idx, b in self.items():
            if b.config.defined_state is False:
                b.mass_frac_phase_comp["Liq", "H2O"].unfix()

        return flags


@declare_process_block_class(
    "SeawaterFTPwStateBlock", block_class=_SeawaterFTPwStateBlock
)
class SeawaterFTPwStateBlockData(SeawaterStateBlockData):
    """A seawater property package."""

    def build(self):
        """Callable method for Block construction."""
        # Do not call build from immediate parent class
        # but instead from its parent (StateBlockData)
        super(SeawaterStateBlockData, self).build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # Add state variables
        self.flow_mass_phase = Var(
            self.params.phase_list,
            initialize=1,
            bounds=(0.0, None),
            units=pyunits.kg / pyunits.s,
            doc="Mass flow rate",
        )

        self.mass_frac_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize={
                ("Liq", "H2O"): 0.965,
                ("Liq", "TDS"): 0.035,
            },
            bounds=(0.0, 1.001),
            units=pyunits.dimensionless,
            doc="Mass fraction",
        )

        self.temperature = Var(
            initialize=298.15,
            bounds=(273.15, 1000),
            units=pyunits.K,
            doc="Temperature",
        )

        self.pressure = Var(
            initialize=101325,
            bounds=(1e3, 5e7),
            units=pyunits.Pa,
            doc="Pressure",
        )

        if not self.config.defined_state:

            @self.Constraint(self.phase_list)
            def eq_sum_mass_frac_phase_comp(b, p):
                return sum(b.mass_frac_phase_comp[p, j] for j in b.component_list) == 1

    # -----------------------------------------------------------------------------
    # Property Methods

    def _flow_mass_phase(self):
        # Overwrite method from base class.
        raise BurntToast(
            "This function should never be called because flow_mass_phase should "
            "have been created as a state variable. If it does not exist as a "
            "state variable, please report this bug to the WaterTAP developers."
        )

    def _flow_mass_phase_comp(self):
        self.flow_mass_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize={("Liq", "H2O"): 0.965, ("Liq", "TDS"): 0.035},
            bounds=(0.0, None),
            units=pyunits.kg / pyunits.s,
            doc="Mass flow rate",
        )

        def rule_flow_mass_phase_comp(b, p, j):
            return b.flow_mass_phase_comp[p, j] == (
                b.flow_mass_phase[p] * b.mass_frac_phase_comp[p, j]
            )

        self.eq_flow_mass_phase_comp = Constraint(
            self.params.phase_list,
            self.params.component_list,
            rule=rule_flow_mass_phase_comp,
        )

    def _mass_frac_phase_comp(self):
        raise BurntToast(
            "This function should never be called because mass_frac_phase_comp should "
            "have been created as a state variable. If it does not exist as a "
            "state variable, please report this bug to the WaterTAP developers."
        )

    def define_state_vars(self):
        """Define state vars."""
        return {
            "flow_mass_phase": self.flow_mass_phase,
            "mass_frac_phase_comp": self.mass_frac_phase_comp,
            "temperature": self.temperature,
            "pressure": self.pressure,
        }

    # -----------------------------------------------------------------------------
    # Scaling methods
    def calculate_scaling_factors(self):
        # Don't call method from immediate parent class, but the
        # one from that class's parent class (StateBlockData)
        super(SeawaterStateBlockData, self).calculate_scaling_factors()

        # setting scaling factors for variables

        # default scaling factors have already been set with
        # idaes.core.property_base.calculate_scaling_factors()
        # for the following variables: flow_mass_phase, mass_frac_phase_comp, pressure,
        # temperature, dens_mass_phase, visc_d_phase, osm_coeff, and enth_mass_phase

        if self.is_property_constructed("flow_mass_phase_comp"):
            for (p, j), vardata in self.flow_mass_phase_comp.items():
                # These scaling factors should have been set either by the
                # user directly or from the default scaling factors
                sf_w = iscale.get_scaling_factor(
                    self.mass_frac_phase_comp[p, j], warning=True
                )
                sf_F = iscale.get_scaling_factor(self.flow_mass_phase[p], warning=True)
                iscale.set_scaling_factor(vardata, sf_w * sf_F, overwrite=False)

        # scaling factors for parameters
        # TODO this should live on the parameter block.
        for j, v in self.params.mw_comp.items():
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(self.params.mw_comp, 1e2)

        # these variables do not typically require user input,
        # will not override if the user does provide the scaling factor
        if self.is_property_constructed("pressure_osm_phase"):
            if iscale.get_scaling_factor(self.pressure_osm_phase["Liq"]) is None:
                iscale.set_scaling_factor(
                    self.pressure_osm_phase["Liq"],
                    iscale.get_scaling_factor(self.pressure),
                )

        if self.is_property_constructed("flow_vol_phase"):
            sf = iscale.get_scaling_factor(
                self.flow_mass_phase_comp["Liq", "H2O"]
            ) / iscale.get_scaling_factor(self.dens_mass_phase["Liq"])
            iscale.set_scaling_factor(self.flow_vol_phase, sf)

        if self.is_property_constructed("flow_vol"):
            sf = iscale.get_scaling_factor(self.flow_vol_phase["Liq"])
            iscale.set_scaling_factor(self.flow_vol, sf)

        if self.is_property_constructed("conc_mass_phase_comp"):
            for j in self.params.component_list:
                sf_dens = iscale.get_scaling_factor(self.dens_mass_phase["Liq"])
                if (
                    iscale.get_scaling_factor(self.conc_mass_phase_comp["Liq", j])
                    is None
                ):
                    if j == "H2O":
                        # solvents typically have a mass fraction between 0.5-1
                        iscale.set_scaling_factor(
                            self.conc_mass_phase_comp["Liq", j], sf_dens
                        )
                    elif j == "TDS":
                        iscale.set_scaling_factor(
                            self.conc_mass_phase_comp["Liq", j],
                            sf_dens
                            * iscale.get_scaling_factor(
                                self.mass_frac_phase_comp["Liq", j]
                            ),
                        )

        if self.is_property_constructed("flow_mol_phase_comp"):
            for j in self.params.component_list:
                if (
                    iscale.get_scaling_factor(self.flow_mol_phase_comp["Liq", j])
                    is None
                ):
                    sf = iscale.get_scaling_factor(self.flow_mass_phase_comp["Liq", j])
                    sf /= iscale.get_scaling_factor(self.params.mw_comp[j])
                    iscale.set_scaling_factor(self.flow_mol_phase_comp["Liq", j], sf)

        if self.is_property_constructed("mole_frac_phase_comp"):
            for j in self.params.component_list:
                if (
                    iscale.get_scaling_factor(self.mole_frac_phase_comp["Liq", j])
                    is None
                ):
                    if j == "TDS":
                        sf = iscale.get_scaling_factor(
                            self.flow_mol_phase_comp["Liq", j]
                        ) / iscale.get_scaling_factor(
                            self.flow_mol_phase_comp["Liq", "H2O"]
                        )
                        iscale.set_scaling_factor(
                            self.mole_frac_phase_comp["Liq", j], sf
                        )
                    elif j == "H2O":
                        iscale.set_scaling_factor(
                            self.mole_frac_phase_comp["Liq", j], 1
                        )

        if self.is_property_constructed("molality_phase_comp"):
            for j in self.params.component_list:
                if isinstance(getattr(self.params, j), Solute):
                    if (
                        iscale.get_scaling_factor(self.molality_phase_comp["Liq", j])
                        is None
                    ):
                        sf = iscale.get_scaling_factor(
                            self.mass_frac_phase_comp["Liq", j]
                        )
                        sf /= iscale.get_scaling_factor(self.params.mw_comp[j])
                        iscale.set_scaling_factor(
                            self.molality_phase_comp["Liq", j], sf
                        )

        if self.is_property_constructed("enth_flow"):
            iscale.set_scaling_factor(
                self.enth_flow,
                iscale.get_scaling_factor(self.flow_mass_phase_comp["Liq", "H2O"])
                * iscale.get_scaling_factor(self.enth_mass_phase["Liq"]),
            )

        if self.is_property_constructed("boiling_point_elevation_phase"):
            iscale.set_scaling_factor(self.boiling_point_elevation_phase["Liq"], 1)

        # transforming constraints
        if self.is_property_constructed("eq_sum_mass_frac_phase_comp"):
            for p, condata in self.eq_sum_mass_frac_phase_comp.items():
                # Register the fact that this constraint is well-scaled
                iscale.constraint_scaling_transform(condata, 1)
        # property relationships with no index, simple constraint
        v_str_lst_simple = [
            "dens_mass_solvent",
            "osm_coeff",
            "pressure_sat",
            "dh_vap_mass",
        ]
        for v_str in v_str_lst_simple:
            if self.is_property_constructed(v_str):
                v = getattr(self, v_str)
                sf = iscale.get_scaling_factor(v, default=1, warning=True)
                c = getattr(self, "eq_" + v_str)
                iscale.constraint_scaling_transform(c, sf)

        if self.is_property_constructed("pressure_osm_phase"):
            sf = iscale.get_scaling_factor(
                self.pressure_osm_phase["Liq"], default=1, warning=True
            )
            iscale.constraint_scaling_transform(self.eq_pressure_osm_phase["Liq"], sf)

        # transforming constraints
        transform_property_constraints(self)
