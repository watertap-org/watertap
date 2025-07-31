#################################################################################
# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

# Model example , which needs to run in QSD-San evironment

from exposan.metab import create_system

__author__ = "Maojian Wang"

inf_fr = 9.18
temp = 24.58
hrt = 11.39
n_stages = 2
reactor_type = "FB"
gas_extraction = "M"
t_span = 200
# set model
sys = create_system(
    n_stages=n_stages,  # number of stages
    reactor_type=reactor_type,  # PB for packed bed, FB for fluidized bed, or UASB
    gas_extraction=gas_extraction,
    Q=inf_fr,  # influent flowrate in m3/d
    T=temp,  # reactor temperature in degree C
    tot_HRT=hrt,  # total HRT in d
)
# run model
sys.simulate(state_reset_hook="reset_cache", t_span=(0, t_span), method="BDF")
