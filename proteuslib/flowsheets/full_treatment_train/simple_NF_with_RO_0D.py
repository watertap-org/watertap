###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################
from pyomo.environ import ConcreteModel, TransformationFactory,\
                          Expression, SolverFactory, Constraint, \
                          Objective, Var, NonNegativeReals, maximize, value, ConstraintList
from pyomo.network import Arc, SequentialDecomposition
from pyomo.environ import units as pyunits
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import Mixer, Separator, Product, Feed
from idaes.generic_models.unit_models.separator import SplittingType, EnergySplittingType
from idaes.generic_models.unit_models.mixer import MixingType
from idaes.core.util.model_statistics import fixed_variables_set, fixed_variables_in_activated_equalities_set, degrees_of_freedom, unfixed_variables_in_activated_equalities_set, activated_equalities_set
import idaes.core.util.scaling as iscale
import proteuslib.property_models.NaCl_CaSO4_prop_pack as nacl_caso4
from proteuslib.util.initialization import check_solve, check_dof
import copy

# from pyomo.environ import (TerminationCondition,
#                            value,
#                            Constraint,
#                            Expression,
#                            Objective,
#                            Param,
#                            TransformationFactory,
#                            units as pyunits)
# import pyomo.util.infeasible as infeas
from idaes.core.util import get_solver
# from idaes.core.util.model_statistics import degrees_of_freedom
# from idaes.core.util.initialization import (solve_indexed_blocks,
#                                             propagate_state,
#                                             fix_state_vars,
#                                             revert_state_vars)
# from idaes.generic_models.unit_models.mixer import MomentumMixingType
import idaes.logger as idaeslog
#
from idaes.generic_models.unit_models.translator import Translator
import proteuslib.property_models.seawater_prop_pack as tds_props
from proteuslib.unit_models.reverse_osmosis_0D import (ReverseOsmosis0D,
                                                       ConcentrationPolarizationType,
                                                       MassTransferCoefficient,
                                                       PressureChangeType)
# from proteuslib.unit_models.pressure_exchanger import PressureExchanger
from proteuslib.unit_models.pump_isothermal import Pump
# import proteuslib.flowsheets.RO_with_energy_recovery.financials as financials

# Set up logger
_log = idaeslog.getLogger(__name__)

# set up solver
solver = get_solver(options={'nlp_scaling_method': 'user-scaling'})

def build_flowsheet():
    ### Build the flowsheet model object ###
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties_NF = nacl_caso4.PropParameterBlock()

    # ==========================================================================
    # Pretreatment train: nanofiltration
    ### Add the mixer block ###
    m.fs.pre_mixer = Mixer(default={
            "property_package": m.fs.properties_NF,
            "inlet_list": ['recycle', 'bypass'],
            "energy_mixing_type": MixingType.none
            })
    ### Add the mixer block ###
    m.fs.post_mixer = Mixer(default={
            "property_package": m.fs.properties_NF,
            "inlet_list": ['pretreatment', 'bypass'],
            "energy_mixing_type": MixingType.none
            })

    ### Add the splitter block ###
    m.fs.bypass = Separator(default={
            "property_package": m.fs.properties_NF,
            "outlet_list": ['pretreatment', 'bypass'],
            "ideal_separation": False,
            "split_basis": SplittingType.totalFlow,
            "energy_split_basis": EnergySplittingType.equal_temperature
            })

    ### Add the separator block ###
    m.fs.NF = Separator(default={
            "property_package": m.fs.properties_NF,
            "outlet_list": ['retentate', 'permeate'],
            "ideal_separation": False,
            "split_basis": SplittingType.componentFlow,
            "energy_split_basis": EnergySplittingType.equal_temperature
            })

    ### Add the splitter block ###
    m.fs.split = Separator(default={
            "property_package": m.fs.properties_NF,
            "outlet_list": ['disposal', 'recycle'],
            "ideal_separation": False,
            "split_basis": SplittingType.totalFlow,
            "energy_split_basis": EnergySplittingType.equal_temperature
            })

    # ==========================================================================
    # Translate NaCl-CaSO4 to TDS property model
    m.fs.properties_RO = tds_props.SeawaterParameterBlock()
    m.fs.translator = Translator(
        default={"inlet_property_package": m.fs.properties_NF,
                 "outlet_property_package": m.fs.properties_RO})
    #TODO: put these constraints on the translator block
    m.fs.translator.constraints = ConstraintList()
    m.fs.translator.constraints.add(m.fs.translator.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'] ==
                                        m.fs.translator.outlet.flow_mass_phase_comp[0, 'Liq', 'H2O'])

    m.fs.translator.constraints.add(m.fs.translator.inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl']
                                   + m.fs.translator.inlet.flow_mass_phase_comp[0, 'Liq', 'CaSO4']
                                   == m.fs.translator.outlet.flow_mass_phase_comp[0, 'Liq', 'TDS'])
    m.fs.translator.constraints.add(m.fs.translator.inlet.temperature[0]
                                   == m.fs.translator.outlet.temperature[0])
    m.fs.translator.constraints.add(m.fs.translator.inlet.pressure[0]
                                       == m.fs.translator.outlet.pressure[0])

    # ==========================================================================
    # Primary treatment train: Reverse osmosis
    m.fs.HP_pump = Pump(default={'property_package': m.fs.properties_RO})
    m.fs.RO = ReverseOsmosis0D(default={
        "property_package": m.fs.properties_RO,
        "has_pressure_change": True,
        "pressure_change_type": PressureChangeType.calculated,
        "mass_transfer_coefficient": MassTransferCoefficient.calculated,
        "concentration_polarization_type": ConcentrationPolarizationType.calculated
    })
    m.fs.permeate = Product(default={'property_package': m.fs.properties_RO})
    m.fs.RO_disposal = Product(default={'property_package': m.fs.properties_RO})
    # ==========================================================================
    # Nanofiltration connections
    #
    #            ----------------(S2)--------------------------------|
    #            |                                                   v
    #  feed --> bypass--(S1)--> pre_mixer --(S3)--> NF --(S5)--> post_mixer -(S7)-> to translator
    #                               ^               |
    #                               |              (S4)
    #                             (S6)              |
    #                               |               V
    #                            recycle <-------split
    #                                               |
    #                                               V
    #                                            disposal
    #
    m.fs.S1 = Arc(source=m.fs.bypass.pretreatment, destination=m.fs.pre_mixer.bypass)
    m.fs.S2 = Arc(source=m.fs.bypass.bypass,       destination=m.fs.post_mixer.bypass)
    m.fs.S3 = Arc(source=m.fs.pre_mixer.outlet,    destination=m.fs.NF.inlet)
    m.fs.S4 = Arc(source=m.fs.NF.retentate,        destination=m.fs.split.inlet)
    m.fs.S5 = Arc(source=m.fs.NF.permeate,         destination=m.fs.post_mixer.pretreatment)
    m.fs.S6 = Arc(source=m.fs.split.recycle,       destination=m.fs.pre_mixer.recycle)
    m.fs.S7 = Arc(source=m.fs.post_mixer.outlet, destination=m.fs.translator.inlet)
    # ==========================================================================
    # Reverse osmosis connections
    #
    #
    # translator--(S8)--> HP_pump --(S9)--> RO_stage --(S11)--> Permeate
    #                                              |
    #                                            (S10)
    #                                              |
    #                                              V
    #                                           disposal
    #
    m.fs.S8 = Arc(source=m.fs.translator.outlet, destination=m.fs.HP_pump.inlet)
    m.fs.S9 = Arc(source=m.fs.HP_pump.outlet, destination=m.fs.RO.inlet)
    m.fs.S10 = Arc(source=m.fs.RO.retentate, destination=m.fs.RO_disposal.inlet)
    m.fs.S11 = Arc(source=m.fs.RO.permeate, destination=m.fs.permeate.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.system_water_recovery = Var(
            initialize=0.5,
            bounds=(0, 1),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc='System Water Recovery')

    ### Add water recovery constraints ###
    m.fs.eq_system_recovery = Constraint(
        expr=m.fs.system_water_recovery ==
             sum(m.fs.RO.permeate.flow_mass_phase_comp[0, 'Liq', k]
                  for k in m.fs.properties_RO.component_list)
              / sum(
                         m.fs.bypass.inlet.flow_mass_phase_comp[0, 'Liq', j] for j in
                         m.fs.properties_NF.component_list))
    #TODO: add SI var and constraints for CaSO4
    scaling_mass_frac = 0.002

    m.fs.eq_no_scaling = Constraint(expr=m.fs.post_mixer.outlet.flow_mass_phase_comp[0, 'Liq', 'CaSO4']
                                         <= scaling_mass_frac * (1 - m.fs.RO.recovery_mass_phase_comp[0, 'Liq', 'H2O']))

    return m

def set_dof(m):
    mass_flow = 1
    x_NaCl = 0.005
    x_CaSO4 = 0.0015
    feed_temp = 298.15
    m.fs.bypass.inlet.pressure.fix(101325)
    m.fs.bypass.inlet.temperature.fix(feed_temp)
    m.fs.bypass.inlet.flow_mass_phase_comp[0.0, 'Liq', 'NaCl'].fix(x_NaCl * mass_flow)
    m.fs.bypass.inlet.flow_mass_phase_comp[0.0, 'Liq', 'CaSO4'].fix(x_CaSO4 * mass_flow)

    m.fs.feed_H2O = Constraint(
        expr=m.fs.bypass.inlet.flow_mass_phase_comp[0.0, 'Liq', 'H2O'] ==
             (
             1.0
             - m.fs.bypass.inlet.flow_mass_phase_comp[0.0, 'Liq', 'NaCl']
             - m.fs.bypass.inlet.flow_mass_phase_comp[0.0, 'Liq', 'CaSO4']
             )
             * mass_flow)


    ### Fix mixer temp out ###
    m.fs.pre_mixer.outlet.temperature.fix(feed_temp)
    m.fs.post_mixer.outlet.temperature.fix(feed_temp)

    ### set NF recover ratio ###
    m.fs.NF.split_fraction[0, "permeate", "H2O"].fix(0.85)
    m.fs.NF.split_fraction[0, "permeate", "NaCl"].fix(0.8)
    m.fs.NF.split_fraction[0, "permeate", "CaSO4"].fix(0.2) # Parameter of Interest

    ### set waste ratio
    m.fs.split.split_fraction[0, "recycle"].fix(0.8)

    ### set bypass ratio
    m.fs.bypass.split_fraction[0, "bypass"].fix(0.5)

    ### Fix Pump and RO variables
    m.fs.HP_pump.efficiency_pump.fix(0.8)  # pump efficiency [-]
    m.fs.HP_pump.control_volume.properties_out[0].pressure.fix(20e5)
    m.fs.RO.A_comp.fix(4.2e-12)  # membrane water permeability coefficient [m/s-Pa]
    m.fs.RO.B_comp.fix(3.5e-8)  # membrane salt permeability coefficient [m/s]
    m.fs.RO.channel_height.fix(1e-3)  # channel height in membrane stage [m]
    m.fs.RO.spacer_porosity.fix(0.85)  # spacer porosity in membrane stage [-]
    m.fs.RO.permeate.pressure[0].fix(101325)  # atmospheric pressure [Pa]
    m.fs.RO.recovery_vol_phase[0, 'Liq'].fix(0.4)  # solvent recovery rate [-]
    m.fs.RO.velocity_io[0, 'out'].fix(0.1)

    m.fs.properties_NF.set_default_scaling('flow_mass_phase_comp', 1/mass_flow, index=('Liq', 'H2O'))
    m.fs.properties_NF.set_default_scaling('flow_mass_phase_comp', 1/mass_flow * 1e2, index=('Liq', 'NaCl'))
    m.fs.properties_NF.set_default_scaling('flow_mass_phase_comp', 1/mass_flow * 1e4, index=('Liq', 'CaSO4'))
    m.fs.properties_RO.set_default_scaling('flow_mass_phase_comp', 1e1, index=('Liq', 'H2O'))
    m.fs.properties_RO.set_default_scaling('flow_mass_phase_comp', 1e3, index=('Liq', 'TDS'))
    iscale.set_scaling_factor(m.fs.HP_pump.control_volume.work, 1e-3)
    iscale.calculate_scaling_factors(m)

    badly_scaled_var_list = list(iscale.badly_scaled_var_generator(m))
    assert badly_scaled_var_list == []

    check_dof(m, fail_flag=True)

def initialize_flowsheet(m):
    m.fs.HP_pump.initialize()

    m.fs.RO.feed_side.properties_in[0].flow_mass_phase_comp['Liq', 'H2O'] = \
        0.97657
        # value(m.fs.translator.outlet.flow_mass_phase_comp[0, 'Liq', 'H2O'])
    m.fs.RO.feed_side.properties_in[0].flow_mass_phase_comp['Liq', 'TDS'] = \
        0.0060476
        # value(m.fs.translator.outlet.flow_mass_phase_comp[0, 'Liq', 'TDS'])
    m.fs.RO.feed_side.properties_in[0].temperature = \
        value(m.fs.translator.outlet.temperature[0])
    m.fs.RO.feed_side.properties_in[0].pressure = \
        value(m.fs.HP_pump.control_volume.properties_out[0].pressure)

    m.fs.RO.area.setub(100)
    m.fs.RO.initialize(fail_on_warning=False, ignore_dof=False)


def solve_flowsheet(m):
    results = solver.solve(m, tee=False)
    check_solve(results, checkpoint="Simulation Solve", fail_flag=True)

def optimize_flowsheet(m):
    #TODO: add objective function

    # [print(i) for i in fixed_variables_set(m)]
    # print(len(fixed_variables_set(m)))
    # print('===========================================================')
    # [print(i) for i in fixed_variables_in_activated_equalities_set(m)]
    # print(len(fixed_variables_in_activated_equalities_set(m)))

    m.fs.RO.velocity_io[0.0, 'out'].unfix()
    m.fs.HP_pump.control_volume.properties_out[0.0].pressure.unfix()
    m.fs.RO.recovery_vol_phase[0, 'Liq'].unfix()

    m.fs.bypass.split_fraction[0.0, 'bypass'].unfix()
    m.fs.split.split_fraction[0, "recycle"].unfix()
    # m.fs.RO.design_permeate_conc = Constraint(
    #     expr=m.fs.RO.permeate_side.properties_mixed[0].mass_frac_phase_comp['Liq', 'TDS'] <= 5e-6)
    m.fs.objective = Objective(expr=-m.fs.system_water_recovery)


    results = solver.solve(m, tee=False)
    check_solve(results, checkpoint="Final optimization", fail_flag=True)

def print_results(m):
    m.fs.RO.report()
    min_recycle_fraction, RO_water_recovery, System_water_recovery, NF_inflow, bypass_fraction = [
        m.fs.split.split_fraction[0, "recycle"].value,
        m.fs.RO.recovery_vol_phase[0, 'Liq'].value,
        m.fs.system_water_recovery.value,
        m.fs.NF.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'].value,
        m.fs.bypass.split_fraction[0, "bypass"].value,
    ]
    print("recycle fraction:      {:>5.2f} %".format(min_recycle_fraction * 100))
    print("bypassed fraction:     {:>5.2f} %".format(bypass_fraction * 100))
    print("RO water recovery:     {:>5.2f} %".format(RO_water_recovery * 100))
    print("System water recovery: {:>5.2f} %".format(System_water_recovery * 100))
    print("NF inflow:             {:>5.2f} kg/s".format(NF_inflow))
    print("-------------------------------------")

def main():
    m = build_flowsheet()

    set_dof(m)

    initialize_flowsheet(m)

    solve_flowsheet(m)

    optimize_flowsheet(m)

    print_results(m)


if __name__ == "__main__":
    main()