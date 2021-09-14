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
from pyomo.environ import (
    Block, Constraint, Expression, Var, Param, Reals, NonNegativeReals, units as pyunits)
from idaes.core.util.exceptions import ConfigurationError


# TODO: choose year --> 2018 probably (use CEPCI)
# TODO: in example flowsheets --> build_costing and use **kwargs to build flowsheet
# TODO: make kwargs dict
# Todo: have options for PX types/pump types (for example) or use more generic approach with conditionals
# mixers, splitters, pumps, erds, RO, NF, stoich reactor (lime softening), equilibrium reactor (chlorination)

def add_costing_param_block(self):
    self.costing_param = Block()
    b = self.costing_param

    b.load_factor = Var(
        initialize=0.9,
        doc='Load factor [fraction of uptime]')
    b.factor_total_investment = Var(
        initialize=2,
        doc='Total investment factor [investment cost/equipment cost]')
    b.factor_MLC = Var(
        initialize=0.03,
        doc='Maintenance-labor-chemical factor [fraction of investment cost/year]')
    b.factor_capital_annualization = Var(
        initialize=0.1,
        doc='Capital annualization factor [fraction of investment cost/year]')
    b.factor_membrane_replacement = Var(
        initialize=0.2,
        doc='Membrane replacement factor [fraction of membrane replaced/year]')
    b.electricity_cost = Var(
        initialize=0.07,
        doc='Electricity cost [$/kWh]')
    b.RO_mem_cost = Var(
        initialize=30,
        doc='Membrane cost [$/m2]')
    b.NF_mem_cost = Var(
        initialize=15,
        doc='Membrane cost [$/m2]')
    b.hp_pump_cost = Var(
        initialize=53 / 1e5 * 3600,
        doc='High pressure pump cost [$/W]')
    b.pxr_cost = Var(
        initialize=535,
        doc='Pressure exchanger cost [$/(m3/h)]')
    b.chemical_lime_cost = Var(
        #TODO: add "real" value instead of dummy val for lime cost per kg
        initialize=1,
        doc='Lime cost [$/kg]')

    # traditional parameters are the only Vars on the block and should be fixed
    for v in b.component_objects(Var, descend_into=True):
        for i in v:
            if v[i].value is None:
                raise ConfigurationError(
                    "{} parameter {} was not assigned"
                    " a value. Please check your configuration "
                    "arguments.".format(b.name, v.local_name))
            v[i].fix()


def get_system_costing(self):
    if not hasattr(self, 'costing'):
        self.costing = Block()
    b = self.costing

    b.capital_cost_total = Var(
        initialize=1e3,
        domain=NonNegativeReals,
        doc='Total capital cost [$]')
    b.investment_cost_total = Var(
        initialize=1e3,
        domain=NonNegativeReals,
        doc='Total investment cost [$]')
    b.operating_cost_MLC = Var(
        initialize=1e3,
        domain=NonNegativeReals,
        doc='Maintenance-labor-chemical operating cost [$/year]')
    b.operating_cost_total = Var(
        initialize=1e3,
        domain=NonNegativeReals,
        doc='Total operating cost [$/year]')
    b.LCOW = Var(
        initialize=1,
        domain=NonNegativeReals,
        doc='Levelized cost of water [$/m3]')

    capital_cost_var_lst = []
    operating_cost_var_lst = []
    for b_unit in self.component_objects(Block, descend_into=True):
        if hasattr(b_unit, 'costing'):
            capital_cost_var_lst.append(b_unit.costing.capital_cost)
            operating_cost_var_lst.append(b_unit.costing.operating_cost)
    operating_cost_var_lst.append(b.operating_cost_MLC)

    b.eq_capital_cost_total = Constraint(
        expr=b.capital_cost_total == sum(capital_cost_var_lst))
    b.eq_investment_cost_total = Constraint(
        expr=(b.investment_cost_total ==
              b.capital_cost_total * self.costing_param.factor_total_investment))
    b.eq_operating_cost_MLC = Constraint(
        expr=(b.operating_cost_MLC ==
              b.investment_cost_total * self.costing_param.factor_MLC))
    b.eq_operating_cost_total = Constraint(
        expr=b.operating_cost_total == sum(operating_cost_var_lst))
    b.eq_LCOW = Constraint(
        expr=b.LCOW == (b.investment_cost_total * self.costing_param.factor_capital_annualization
                        + b.operating_cost_total) / (self.annual_water_production / (pyunits.m ** 3 / pyunits.year)))


def _make_vars(self):
    # build generic costing variables (all costing models need these vars)
    self.capital_cost = Var(initialize=1e5,
                            domain=NonNegativeReals,
                            doc='Unit capital cost [$]')
    self.operating_cost = Var(initialize=1e5,
                              domain=Reals,
                              bounds=(0, 1e6),
                              doc='Unit operating cost [$/year]')


def ReverseOsmosis_costing(self):
    _make_vars(self)

    b_RO = self.parent_block()
    b_fs = b_RO.parent_block()

    # capital cost
    self.eq_capital_cost = Constraint(
        expr=self.capital_cost == b_fs.costing_param.RO_mem_cost * b_RO.area / pyunits.m ** 2)

    # operating cost
    self.eq_operating_cost = Constraint(
        expr=self.operating_cost == b_fs.costing_param.factor_membrane_replacement
             * b_fs.costing_param.RO_mem_cost * b_RO.area / pyunits.m ** 2)


def Nanofiltration_costing(self):
    ''' This method is being added for the nanofiltration step in the pre-treatment section of the full treatment train'''

    _make_vars(self)

    b_NF = self.parent_block()
    b_fs = b_NF.parent_block()

    # capital cost
    self.eq_capital_cost = Constraint(
        expr=self.capital_cost == b_fs.costing_param.NF_mem_cost * b_NF.area / pyunits.m ** 2)

    # operating cost
    self.eq_operating_cost = Constraint(
        expr=self.operating_cost == b_fs.costing_param.factor_membrane_replacement
             * b_fs.costing_param.NF_mem_cost * b_NF.area / pyunits.m ** 2)

def Separator_costing(self):
    # TODO: Get separator cost from Tim's single stage MD paper

    pass

def Mixer_costing(self, mixer_type='default'):
    _make_vars(self)

    b_m = self.parent_block()
    b_fs = b_m.parent_block()

    if mixer_type == 'default':
        #TODO: Get mixer cost from Tim's single stage MD paper
        # self.eq_capital_cost = Constraint()
        self.capital_cost.fix(0)
        self.operating_cost.fix(0)

    elif mixer_type == 'naocl_mixer':
        '''Cost estimation of chlorination step for disinfection in post-treatment
        Digitized Fig. 4.19 in Voutchkov, 2018 using WebPlotDigitizer, https://apps.automeris.io/wpd/,
         September 2021. Fig. 4.19 provides construction cost as a function of daily desalination plant capacity.
         Curves for sodium hypochlorite and chlorine dioxide are provided, but only NaOCl data were extracted.
         Data were converted to specific construction costs as a function of capacity to get the provided cost curve
         for capex (of the form a*X**b).
         Since cost figures are reported for the year 2018, the capex cost constraint is assumed to be in 2018 USD;the 
         cost escalation factor, cost_esc, can be modified to account for changes over time.'''

        # NaOCl specific capex ($/m3/day) = 479.87 * x ** (-0.396) ; x is plant capacity (m3/day)
        # TODO: may need to touch flow_vol while building naocl_mixer_unit. Double-check. Alternative: use flow_vol of RO final permeate
        self.cost_esc = Param(initialize=1, mutable=True)
        self.eq_capital_cost = Constraint(expr=self.capital_cost ==
                                               479.87
                                               * (b_m.inlet_stream_state[0].flow_vol*3600*24) ** 0.604
                                               * self.cost_esc)

        # Sodium hypochlorite cost taken from WaterTAP (2020 USD) which assumes 15% purity
        self.naocl_cost = Param(initialize=0.23, mutable=True, units=pyunits.kg**-1)
        self.naocl_purity = Param(initialize=0.15, mutable=True)
        #TODO: no electricity cost included -- would be based on pump work, expected to be negligible for the time being
        self.eq_operating_cost = Constraint(expr=self.operating_cost ==
                                                 b_m.naocl_stream.flow_mol[0]
                                                 * b_m.naocl_stream.mole_frac_comp[0, "OCl_-"]
                                                 * 74.44e-3 * pyunits.kg / pyunits.mol
                                                 * self.naocl_cost
                                                 / self.naocl_purity
                                                 * 3600 * 8760 * pyunits.s)

    elif mixer_type == 'lime_softening':
        '''Cost estimation of lime addition for precipitation step in pretreatment
        Digitized Fig. 5.59 in McGivney & Kawamura, 2008 using WebPlotDigitizer, https://apps.automeris.io/wpd/,
         September 2021. WaterTAP provides a similar equation for cost capacity curve based on the same reference.
         This is suspected to be due to digitization of the reference data although WaterTAP's documentation indicates that
         cost indices were accounted for. The original reference cites an ENR CCI = 8889, representative of the construction cost index 
         for Los Angeles in April 2007. Since recent year data for ENR's CCI are currently unknown, WaterTAP's equation will
         be used; note that WT's documentation may need to correct reporting of units/methodology.
         ------------------------------------------------------------------------------------------
         Cost capacity equations to consider:
         McGivney & Kawamura: 
         12985*x**0.5901
         =====================================================
         WaterTAP:
         16972*x**0.5435
         =====================================================
         Manual digitization and fitting relationships to original reference:
         1) Power law fit (coefficient of determination = 0.9721)
         13310*x**0.5855 
         2) Polynomial fit (coefficient of determination = 1)
         -1.4071*x**2 + 1661.9*x + 40782
         Although the polynomial fit matches the data more closely than other relationships, going above ~ 700 lb/day would 
         result in erroneous decrease in cost.
         '''
        # x is converts mol/s to lb/day
        self.lime_lbs_per_day = Expression(expr=2.205 * 3600 * 24 * 74.09e-3
            * b_m.lime_stream.flow_mol[0].value
            * b_m.lime_stream.mole_frac_comp[0, "Ca(OH)2"].value / pyunits.mol / pyunits.s)

        self.cost_esc = Param(initialize=1, mutable=True)
        self.eq_capital_cost = Constraint(expr=self.capital_cost == 16972 * self.lime_lbs_per_day ** 0.5435 * self.cost_esc)

        # Calcium hydroxide (lime) cost taken from WaterTAP (2020 USD) which assumes 100% purity
        self.caoh2_cost = Param(initialize=0.15, mutable=True, units=pyunits.kg**-1)
        self.caoh2_purity = Param(initialize=1, mutable=True)
        #TODO: no electricity cost included -- would be based on pump work, expected to be negligible for the time being
        self.eq_operating_cost = Constraint(expr=self.operating_cost == b_m.lime_stream.flow_mol[0]
                                                 * b_m.lime_stream.mole_frac_comp[0, "Ca(OH)2"]
                                                 * 74.093e-3 * pyunits.kg / pyunits.mol
                                                 * self.caoh2_cost
                                                 / self.caoh2_purity
                                                 * 3600 * 8760 * pyunits.s)

def PressureExchanger_costing(self):
    _make_vars(self)

    b_PXR = self.parent_block()
    b_fs = b_PXR.parent_block()

    # capital cost
    self.eq_capital_cost = Constraint(
        expr=self.capital_cost == b_fs.costing_param.pxr_cost
             * b_PXR.low_pressure_side.properties_in[0].flow_vol * 3600 / (pyunits.m ** 3 / pyunits.s))

    # operating cost
    self.operating_cost.fix(0)


def pressure_changer_costing(self, pump_type="centrifugal"):
    _make_vars(self)

    b_PC = self.parent_block()
    b_fs = b_PC.parent_block()

    self.purchase_cost = Var()
    self.cp_cost_eq = Constraint(expr=self.purchase_cost == 0)

    if pump_type == 'High pressure':
        # capital cost
        self.eq_capital_cost = Constraint(
            expr=self.capital_cost == b_fs.costing_param.hp_pump_cost * b_PC.work_mechanical[0] / pyunits.W)

        # operating cost
        self.eq_operating_cost = Constraint(
            expr=self.operating_cost == (b_PC.work_mechanical[0] / pyunits.W
                                         * 3600 * 24 * 365 * b_fs.costing_param.load_factor)
                 * b_fs.costing_param.electricity_cost / 3600 / 1000)

    elif pump_type == 'Pressure exchanger':
        # capital cost
        b_cv_in = b_PC.control_volume.properties_in[0]
        self.eq_capital_cost = Constraint(
            expr=(self.capital_cost == b_fs.costing_param.erd_cost['A']
                  * (sum(b_cv_in.flow_mass_comp[j] / (pyunits.kg / pyunits.s)
                         for j in b_PC.config.property_package.component_list)
                     / (b_cv_in.dens_mass / (pyunits.kg / pyunits.m ** 3)) * 3600) ** 0.58))

        # operating cost
        self.operating_cost.setlb(-1e6)
        self.eq_operating_cost = Constraint(
            expr=self.operating_cost == (b_PC.work_mechanical[0] / pyunits.W
                                         * 3600 * 24 * 365 * b_fs.costing_param.load_factor)
                 * b_fs.costing_param.electricity_cost / 3600 / 1000)
    #TODO: add cost relationship for low pressure pump, intended for NF pump