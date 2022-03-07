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

from copy import deepcopy
# Import Pyomo libraries
from pyomo.environ import (Block,
                           Set,
                           Var,
                           Param,
                           Constraint,
                           Expression,
                           Suffix,
                           NonNegativeReals,
                           Reference,
                           units as pyunits)
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (ControlVolume0DBlock,
                        declare_process_block_class,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        UnitModelBlockData,
                        useDefault,
                        MaterialFlowBasis)
from idaes.core.util import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.constants import Constants
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog


_log = idaeslog.getLogger(__name__)

# when using this file the name "Filtration" is what is imported
@declare_process_block_class("Crystallization")
class CrystallizationData(UnitModelBlockData):
    """
    Zero order crystallization model
    """
    # CONFIG are options for the unit model, this simple model only has the mandatory config options
    CONFIG = ConfigBlock()

    CONFIG.declare("dynamic", ConfigValue(
        domain=In([False]),
        default=False,
        description="Dynamic model flag - must be False",
        doc="""Indicates whether this model will be dynamic or not,
    **default** = False. The filtration unit does not support dynamic
    behavior, thus this must be False."""))

    CONFIG.declare("has_holdup", ConfigValue(
        default=False,
        domain=In([False]),
        description="Holdup construction flag - must be False",
        doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False. The filtration unit does not have defined volume, thus
    this must be False."""))

    CONFIG.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations,
    **default** - useDefault.
    **Valid values:** {
    **useDefault** - use default package from parent model or flowsheet,
    **PhysicalParameterObject** - a PhysicalParameterBlock object.}"""))

    CONFIG.declare("property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
    and used when constructing these,
    **default** - None.
    **Valid values:** {
    see property package for documentation.}"""))


    def build(self):
        super().build()

        solvent_set = self.config.property_package.solvent_set
        solute_set = self.config.property_package.solute_set

        # this creates blank scaling factors, which are populated later
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # Next, get the base units of measurement from the property definition
        units_meta = self.config.property_package.get_metadata().get_derived_units

        # Add unit variables
        ## TO-DO: Add slurry density?

        self.k_param = Param(
            initialize=0.04,
            units=pyunits.m/pyunits.s,
            doc="Sounders-Brown constant, set at 0.04 m/s based on Dutta et al. \
            Lewis et al suggests 0.024 m/s, while Tavare suggests about 0.06 m/s ")

        self.hex_max_temp_increase = Param(
            initialize=4,
            units=pyunits.K,
            doc="Maximum temperature difference between inlet and outlet of a crystallizer heat exchanger.\
            Lewis et al. suggests 1-2 degC but use 5degC in example; Tavare example used 4 degC.\
            Default is 4 degC")


        self.crystallization_yield = Var(
            solute_set,
            initialize=0.5,
            bounds=(1e-8, 1),
            units=pyunits.dimensionless,
            doc="Crystallizer solids yield"
            )

        self.product_volumetric_solids_fraction = Var(
            initialize=0.25,
            bounds=(1e-8, 1),
            units=pyunits.dimensionless,
            doc = "Volumetric fraction of solids in slurry product (i.e. solids-liquid mixture)."
            )

        self.temperature_operating = Var(
            initialize=298.15,
            bounds=(273, 1000),
            units=pyunits.K,
            doc="Crystallizer operating temperature: boiling point of the solution."
            )

        self.pressure_operating = Var(
            initialize=1e3,
            bounds=(0.001, 1e6),
            units=pyunits.Pa,
            doc="Operating pressure of the crystallizer."
            )

        self.magma_density = Var(
            initialize=250,
            bounds=(1, 5000),
            units=pyunits.kg/pyunits.m**3,
            doc="Magma density, i.e. mass of crystals per unit volume of suspension"
            )

        self.slurry_density = Var(
            initialize=1000,
            bounds=(1, 5000),
            units=pyunits.kg/pyunits.m**3,
            doc="Suspension density, i.e. density of solid-liquid mixture before separation"
            )

        self.work_mechanical = Var(
            self.flowsheet().config.time,
            initialize=1e5,
            bounds=(-5e6, 5e6),
            units=pyunits.kJ/pyunits.s,
            doc="Crystallizer thermal energy requirement"
            )

        self.d_crystallizer = Var(
            initialize=3,
            bounds=(0, 25),
            units=pyunits.m,
            doc="Minimum diameter of crystallizer"
            )

        self.hex_circulation_flow_vol = Var(
            initialize=1,
            bounds=(0, 100),
            units=pyunits.m**3/pyunits.s,
            doc="Minimum circulation flow rate through crystallizer heat exchanger"
            )

        self.relative_supersaturation = Var(
            solute_set,
            initialize=0.1,
            bounds=(0, 100),
            units=pyunits.dimensionless
        )

        # Add state blocks for inlet, outlet, and waste
        # These include the state variables and any other properties on demand
        # Add inlet block
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True  # inlet block is an inlet
        self.properties_in = self.config.property_package.state_block_class(
            self.flowsheet().config.time,  # time domain for the state block, just 0 in this case
            doc="Material properties of inlet",
            default=tmp_dict)

        # Add outlet and waste block
        tmp_dict["defined_state"] = False  # outlet and waste block is not an inlet
        self.properties_out = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of liquid outlet",
            default=tmp_dict)

        self.properties_solids = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of solid crystals at outlet",
            default=tmp_dict)

        self.properties_vapor = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of water vapour at outlet",
            default=tmp_dict)

        # Add ports - oftentimes users interact with these rather than the state blocks
        self.add_port(name='inlet', block=self.properties_in)
        self.add_port(name='outlet', block=self.properties_out)
        self.add_port(name='solids', block=self.properties_solids)
        self.add_port(name='vapor', block=self.properties_vapor)



        # Add constraints
        # 1. Material balances

        ## (i) Solute mass balances
        @self.Constraint(self.config.property_package.component_list,
                         doc="Mass balance for components")
        def eq_mass_balance_solutes(b, j):
            if j in solute_set: 
                return (b.properties_in[0].flow_mass_phase_comp['Liq', j] + b.properties_in[0].flow_mass_phase_comp['Sol', j] \
                    == b.properties_out[0].flow_mass_phase_comp['Liq', j] + b.properties_solids[0].flow_mass_phase_comp['Sol', j])
            else:
                return Constraint.Skip

        ## (ii) Solvent mass balances
        @self.Constraint(doc="Mass balance for solvent")
        def eq_mass_balance_solvent(b):
                return (b.properties_in[0].flow_mass_phase_comp['Liq', 'H2O'] 
                    == b.properties_out[0].flow_mass_phase_comp['Liq', 'H2O']
                    + b.properties_vapor[0].flow_mass_phase_comp['Vap', 'H2O'])               

        # ## (iii) Total solids composition balance
        # @self.Constraint(doc="Total crystal mass")
        # def eq_total_crystal_mass(b):
        #     return (sum(b.properties_solids[0].mass_frac_phase_comp['Sol', j] for j in solute_set) == 1)


        # 2. Constraint on outlet liquid composition based on solubility requirements
        @self.Constraint(self.config.property_package.component_list, doc="Solubility vs mass fraction constraint")
        def eq_solubility_massfrac_equality_constraint(b, j):
            if j in solute_set:
                return (b.properties_out[0].mass_frac_phase_comp['Liq', j] - b.properties_out[0].solubility_mass_frac_phase_comp['Liq', j] == 0)
            else:
                return Constraint.Skip


        # 3. Performance equations 
        # (a) based on yield
        @self.Constraint(self.config.property_package.component_list,
                         doc="Component salt yield equation")
        def eq_removal_balance(b, j):
            if j in solvent_set:
                return Constraint.Skip
            else:
                return (b.properties_in[0].flow_mass_phase_comp['Liq', j] * b.crystallization_yield[j]
                    == b.properties_in[0].flow_mass_phase_comp['Liq', j] - b.properties_out[0].flow_mass_phase_comp['Liq', j])

        # (b) Volumetric fraction constraint
        @self.Constraint(doc="Solid volumetric fraction in outlet: constraint, 1-E")
        def eq_vol_fraction_solids(b):
            return (self.product_volumetric_solids_fraction 
                == b.properties_solids[0].flow_vol / (b.properties_solids[0].flow_vol + b.properties_out[0].flow_vol)
                )

        # (c) Magma density constraint
        @self.Constraint(doc="Slurry magma density")
        def eq_dens_magma(b):
            return (self.magma_density == 
                b.properties_solids[0].dens_mass_solute['Sol'] * self.product_volumetric_solids_fraction
                )

        # (d) Operating pressure constraint
        @self.Constraint(doc="Operating pressure constraint")
        def eq_operating_pressure_constraint(b):
            return (self.pressure_operating - b.properties_out[0].pressure_sat == 0)

        # (d) Relative supersaturation 
        @self.Constraint(solute_set, doc="Relative supersaturation created via evaporation, g/g (solution)")
        def eq_relative_supersaturation(b, j):
            #  mass_frac_after_evap = SOLIDS IN + LIQUID IN - VAPOUR OUT
            mass_frac_after_evap = b.properties_in[0].flow_mass_phase_comp['Liq', j] \
            / (sum(b.properties_in[0].flow_mass_phase_comp['Liq', j] for j in solute_set)
                + b.properties_in[0].flow_mass_phase_comp['Liq', 'H2O']
                - b.properties_vapor[0].flow_mass_phase_comp['Vap', 'H2O']
                )
            # return (b.relative_supersaturation[j] * b.properties_out[0].solubility_mass_frac_phase_comp['Liq', j] == 
            # (mass_frac_after_evap - b.properties_out[0].solubility_mass_frac_phase_comp['Liq', j])
            # )
            return (b.relative_supersaturation[j] == 
                (mass_frac_after_evap - b.properties_out[0].solubility_mass_frac_phase_comp['Liq', j])
                / b.properties_out[0].solubility_mass_frac_phase_comp['Liq', j]
                )






        # 4. Fix flows of empty solid, liquid and vapour streams
        # (i) Fix solids: liquid and vapour flows must be zero
        @self.Constraint(self.properties_solids[0].phase_component_set,
            doc="Empty components in solids stream")
        def eq_solids_composition(b, p, j):
            if p != 'Sol':
                return (b.properties_solids[0].flow_mass_phase_comp[p, j] == 1e-8 * pyunits.kg/pyunits.s)
            else:
                return Constraint.Skip

        # (ii) Fix liquids: solid and vapour flows must be zero
        @self.Constraint(self.properties_out[0].phase_component_set,
            doc="Empty components in liquid stream")
        def eq_liquid_composition(b, p, j):
            if p != 'Liq':
                return (b.properties_out[0].flow_mass_phase_comp[p, j] == 1e-8 * pyunits.kg/pyunits.s)
            else:
                return Constraint.Skip

        # (iii) Fix vapor: solid and vapour flows must be zero. Vapour is also pure solvent.
        @self.Constraint(self.properties_vapor[0].phase_component_set,
            doc="Empty components in vapour stream")
        def eq_vapor_composition(b, p, j):
            if p != 'Vap':
                return (b.properties_vapor[0].flow_mass_phase_comp[p, j] == 1e-8 * pyunits.kg/pyunits.s)
            else:
                return Constraint.Skip


        # 5. Add an energy balance for the system
        ## (iii) Enthalpy balance: based on Lewis et al. Enthalpy is exothermic, and the hC in property package is -ve
        @self.Constraint(doc="Enthalpy balance over crystallization system")
        def eq_enthalpy_balance(b):
            return (b.properties_in[0].enth_flow 
                - b.properties_out[0].enth_flow - b.properties_vapor[0].enth_flow - b.properties_solids[0].enth_flow
                + self.work_mechanical[0]
                - (b.properties_solids[0].flow_vol_phase['Sol'] * b.properties_solids[0].dens_mass_solute['Sol'] * b.properties_solids[0].dh_crystallization['NaCl'])
                == 0
                )


        # 6. Pressure and temperature balances - what is the pressure of the outlet solid and vapour?
        # TO-DO: Figure out actual liquid and solid pressures.
        @self.Constraint()    
        def eq_p_con1(b):
            return b.properties_in[0].pressure - b.properties_out[0].pressure == 0

        @self.Constraint()    
        def eq_p_con2(b):
            return b.properties_in[0].pressure - b.properties_solids[0].pressure == 0

        @self.Constraint()    
        def eq_p_con3(b):
            return b.properties_vapor[0].pressure - self.pressure_operating == 0

        @self.Constraint()    
        def eq_T_con1(b):
            return self.temperature_operating - b.properties_solids[0].temperature == 0

        @self.Constraint()    
        def eq_T_con2(b):
            return self.temperature_operating - b.properties_vapor[0].temperature == 0

        @self.Constraint()
        def eq_T_con3(b):
            return self.temperature_operating - b.properties_out[0].temperature == 0


        # 7. Heat exchanger minimum circulation flow rate calculations - see Lewis et al. or Tavare et al.
        @self.Constraint(doc="Constraint on mimimum circulation rate through crystallizer heat exchanger")
        def eq_minimum_hex_circulation_rate_constraint(b):
            rho_cp_dt = self.hex_max_temp_increase * \
            (
                b.product_volumetric_solids_fraction * b.properties_solids[0].dens_mass_solute['Sol'] * b.properties_solids[0].cp_solute['Sol'] 
                + (1 - b.product_volumetric_solids_fraction) * b.properties_out[0].dens_mass_phase['Liq'] * b.properties_out[0].cp_phase['Liq']
                )
            return (
                b.hex_circulation_flow_vol * rho_cp_dt == pyunits.convert(b.work_mechanical[0], to_units=pyunits.J / pyunits.s)
                )

        # 8. Suspension density
        @self.Constraint(doc="Slurry density calculation")
        def eq_slurry_density(b):
            return (self.slurry_density == 
                b.product_volumetric_solids_fraction * b.properties_solids[0].dens_mass_solute['Sol'] 
                + (1 - b.product_volumetric_solids_fraction) * b.properties_out[0].dens_mass_phase['Liq']
                )



        # 9. Minimum diameter of evaporation zone
        def eq_max_allowable_velocity(b):
            return b.k_param * \
            (b.properties_out[0].dens_mass_phase['Liq']/ b.properties_vapor[0].dens_mass_solvent['Vap']) ** 0.5
        
        self.v_max = Expression(rule=eq_max_allowable_velocity, doc='maximum allowable vapour linear velocity in m/s') 

        @self.Constraint(doc="Mimimum diameter of evaporation zone")
        def eq_vapor_head_diameter_constraint(b):
            return self.d_crystallizer == \
            (4 * b.properties_vapor[0].flow_vol_phase['Vap'] / (Constants.pi * b.v_max)) ** 0.5






    def initialize(
                blk,
                state_args=None,
                outlvl=idaeslog.NOTSET,
                solver=None,
                optarg=None):
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

            opt = get_solver(solver, optarg)

            # ---------------------------------------------------------------------
            # Initialize holdup block
            flags = blk.properties_in.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                state_args=state_args,
                hold_state=True)
            init_log.info_high("Initialization Step 1 Complete.")
            # ---------------------------------------------------------------------
            # Initialize other state blocks
            # Set state_args from inlet state
            if state_args is None:
                state_args = {}
                state_dict = blk.properties_in[
                    blk.flowsheet().config.time.first()].define_port_members()

                for k in state_dict.keys():
                    if state_dict[k].is_indexed():
                        state_args[k] = {}
                        for m in state_dict[k].keys():
                            state_args[k][m] = state_dict[k][m].value
                    else:
                        state_args[k] = state_dict[k].value

            blk.properties_out.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                state_args=state_args,
            )

            state_args_solids = deepcopy(state_args)
            for p, j in blk.properties_solids.phase_component_set:
                if p == 'Sol':
                    state_args_solids['flow_mass_phase_comp'][p, j] = state_args['flow_mass_phase_comp']['Liq', j]
                elif p == 'Liq' or p == 'Vap':
                    state_args_solids['flow_mass_phase_comp'][p, j] = 1e-8
            blk.properties_solids.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                state_args=state_args_solids,
            )

            state_args_vapor = deepcopy(state_args)
            for p, j in blk.properties_vapor.phase_component_set:
                if p == 'Vap':
                    state_args_vapor['flow_mass_phase_comp'][p, j] = state_args['flow_mass_phase_comp']['Liq', j]
                elif p == 'Liq' or p == 'Sol':
                    state_args_vapor['flow_mass_phase_comp'][p, j] = 1e-8
            blk.properties_vapor.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                state_args=state_args_vapor,
            )
            init_log.info_high("Initialization Step 2 Complete.")
            # ---------------------------------------------------------------------
            # Solve unit
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = opt.solve(blk, tee=slc.tee)
            init_log.info_high(
                "Initialization Step 3 {}.".format(idaeslog.condition(res)))
            # ---------------------------------------------------------------------
            # Release Inlet state
            blk.properties_in.release_state(flags, outlvl=outlvl)
            init_log.info(
                "Initialization Complete: {}".format(idaeslog.condition(res))
            )

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        iscale.set_scaling_factor(self.crystallization_yield, 1)
        iscale.set_scaling_factor(self.product_volumetric_solids_fraction, 1)
        iscale.set_scaling_factor(self.temperature_operating, iscale.get_scaling_factor(self.properties_in[0].temperature))
        iscale.set_scaling_factor(self.pressure_operating, 1e-3)#iscale.get_scaling_factor(self.properties_in[0].pressure_sat))
        iscale.set_scaling_factor(self.magma_density, iscale.get_scaling_factor(self.properties_in[0].dens_mass_phase['Liq']))
        iscale.set_scaling_factor(self.slurry_density, iscale.get_scaling_factor(self.properties_in[0].dens_mass_phase['Liq']))
        # Temporary scale for work - scale by vapor mass*enthalpy because vapor typically has the largest enthalpy
        iscale.set_scaling_factor(self.work_mechanical[0], 
                iscale.get_scaling_factor(self.properties_in[0].flow_mass_phase_comp['Vap', 'H2O']) 
                * iscale.get_scaling_factor(self.properties_in[0].enth_mass_solvent['Vap'])) 


        # transforming constraints
        for ind, c in self.eq_T_con1.items():
            sf = iscale.get_scaling_factor(self.properties_in[0].temperature)
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_T_con2.items():
            sf = iscale.get_scaling_factor(self.properties_in[0].temperature)
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_T_con3.items():
            sf = iscale.get_scaling_factor(self.properties_in[0].temperature)
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_p_con1.items():
            sf = iscale.get_scaling_factor(self.properties_in[0].pressure)
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_p_con2.items():
            sf = iscale.get_scaling_factor(self.properties_in[0].pressure)
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_p_con3.items():
            sf = iscale.get_scaling_factor(self.properties_in[0].pressure)
            iscale.constraint_scaling_transform(c, sf)

        for j, c in self.eq_mass_balance_solutes.items():
            sf = iscale.get_scaling_factor(self.properties_in[0].flow_mass_phase_comp['Liq', j])
            iscale.constraint_scaling_transform(c, sf)

        for j, c in self.eq_mass_balance_solvent.items():
            sf = iscale.get_scaling_factor(self.properties_in[0].flow_mass_phase_comp['Liq', 'H2O'])
            iscale.constraint_scaling_transform(c, sf)

        for j, c in self.eq_solubility_massfrac_equality_constraint.items():
            iscale.constraint_scaling_transform(c, 1e0)

        for j, c in self.eq_dens_magma.items():
            sf = iscale.get_scaling_factor(self.properties_in[0].dens_mass_solvent['Liq'])
            iscale.constraint_scaling_transform(c, sf)


        for j, c in self.eq_removal_balance.items():
            sf = iscale.get_scaling_factor(self.properties_in[0].flow_mass_phase_comp['Liq', j])
            iscale.constraint_scaling_transform(c, sf)


    #     # Scaling for energy balance needed

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Feed Inlet": self.inlet,
                "Liquid Outlet": self.outlet,
                "Vapor Outlet": self.vapor,
                "Solid Outlet": self.solids,
            },
            time_point=time_point,
        )

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        var_dict["Operating Temperature (K)"] = self.temperature_operating
        var_dict["Operating Pressure (Pa)"] = self.pressure_operating
        var_dict["Magma density of solution (Kg/m**3)"] = self.magma_density
        var_dict["Slurry density (Kg/m3)"] = self.slurry_density
        var_dict["Heat requirement (kW)"] = self.work_mechanical[time_point]
        var_dict["Crystallizer diameter (m)"] = self.d_crystallizer
        var_dict["Heat exchanger solution-side flow rate (m**3/s)"] = self.hex_circulation_flow_vol
        var_dict["Vol. frac. of solids in suspension, 1-E"] = self.product_volumetric_solids_fraction
        for j in self.config.property_package.solute_set:
            yield_mem_name=f'{j} yield (fraction)'
            var_dict[yield_mem_name] = (self.crystallization_yield[j])
            supersat_mem_name=f'{j} relative supersaturation (g/g solution)'
            var_dict[supersat_mem_name] = (self.relative_supersaturation[j])

        return {"vars": var_dict}