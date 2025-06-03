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

# input space generation

import csv
import pandas as pd
from exposan.metab import create_system


__author__ = "Maojian Wang"


def get_input_data(filename=None):
    if filename == None:
        # It is a test mode
        # Define data
        input_data = {
            "inf_fr": [5, 5, 5],
            "temp": [20, 25, 30],  # First column data
            "hrt": [12, 13, 14],  # Second column data
        }
        # Create the DataFrame
        df = pd.DataFrame(input_data)
    else:
        df = pd.read_csv(filename, header=0)

    return df


# define output
import pandas as pd


def get_eff_fr(case=None, df=None):
    if case is None:
        print(" The system is off")
    else:
        eff_dg = case.outs[3]
        r = eff_dg._info(
            layout=None,
            T=None,
            P=None,
            flow=None,
            composition=None,
            N=100,
            IDs=None,
        )
        # print(r)
        comp = [i.strip() for i in r.split("flow (kmol/hr):")[1].split("\n")]
        # print(comp)
        # print(len(comp))

        fr_dict = {}
        for element in comp:
            element_list = element.split()
            # print(element_list)
            fr_dict[element_list[0]] = [float(element_list[-1])]
            # print(fr_dict)
        if df == None:
            df = pd.DataFrame(fr_dict)
            # print(df)
        else:
            df.loc[len(df)] = fr_dict
            # NEED to TEST
    return df


def get_ch4_fr(case=None, df=None):
    if case is None:
        print(" The system is off")
    else:
        # bge : biogas extracted in reactor 2
        eff_dg = case.outs[2]
        r = eff_dg._info(
            layout=None,
            T=None,
            P=None,
            flow=None,
            composition=None,
            N=100,
            IDs=None,
        )
        # print(r)
        comp = [i.strip() for i in r.split("flow (kmol/hr):")[1].split("\n")]
        # print(comp)
        # print(len(comp))

        fr_dict = {}
        for element in comp:
            element_list = element.split()
            # print(element_list)
            fr_dict[element_list[0]] = [float(element_list[-1])]
            # print(fr_dict)
        if df == None:
            df = pd.DataFrame(fr_dict)
            # print(df)
        else:
            df.loc[len(df)] = fr_dict
            # NEED to TEST
    return df


def get_h2_fr(case=None, df=None):
    if case is None:
        print(" The system is off")
    else:
        # bg2: biogas from reactor 2
        eff_dg = case.outs[1]
        r = eff_dg._info(
            layout=None,
            T=None,
            P=None,
            flow=None,
            composition=None,
            N=100,
            IDs=None,
        )
        # print(r)
        comp = [i.strip() for i in r.split("flow (kmol/hr):")[1].split("\n")]
        # print(comp)
        # print(len(comp))

        fr_dict = {}
        for element in comp:
            element_list = element.split()
            # print(element_list)
            fr_dict[element_list[0]] = [float(element_list[-1])]
            # print(fr_dict)
        if df == None:
            df = pd.DataFrame(fr_dict)
            # print(df)
        else:
            df.loc[len(df)] = fr_dict
            # NEED to TEST
    return df


def get_r1_ex_biogas_fr(case=None, df=None):
    if case is None:
        print(" The system is off")
    else:
        # bgs: biogas extracted from reactor 1
        eff_dg = case.outs[4]
        r = eff_dg._info(
            layout=None,
            T=None,
            P=None,
            flow=None,
            composition=None,
            N=100,
            IDs=None,
        )
        # print(r)
        comp = [i.strip() for i in r.split("flow (kmol/hr):")[1].split("\n")]
        # print(comp)
        # print(len(comp))

        fr_dict = {}
        for element in comp:
            element_list = element.split()
            # print(element_list)
            fr_dict[element_list[0]] = [float(element_list[-1])]
            # print(fr_dict)
        if df == None:
            df = pd.DataFrame(fr_dict)
            # print(df)
        else:
            df.loc[len(df)] = fr_dict
            # NEED to TEST
    return df


def get_mass_flowrate(case=None, df=None, stream=None):
    # df = None
    # stream = eff_dg
    keys = str(stream.components)[20:-2].split(",")
    values = [[float(x)] for x in list(stream.conc.to_array())]
    fr_dict = dict(zip(keys, values))
    if len(stream.state) != len(stream.conc.to_array()):
        fr_dict["Volumetric Flowrate"] = [stream.state[-1]]
    print(fr_dict)
    if df == None:
        df = pd.DataFrame(fr_dict)
        # print(df)
    else:
        df.loc[len(df)] = fr_dict
    return df


def collect_results(case=None, results=None, mass=True):
    if mass:
        eff_fr = get_mass_flowrate(case, None, case.outs[3])
        ch4_fr = get_mass_flowrate(case, None, case.outs[2])
        h2_fr = get_mass_flowrate(case, None, case.outs[1])
        bgs_fr = get_mass_flowrate(case, None, case.outs[4])
    else:
        eff_fr = get_eff_fr(case)
        ch4_fr = get_ch4_fr(case)
        h2_fr = get_h2_fr(case)
        bgs_fr = get_r1_ex_biogas_fr(case)

    ch4_fr.columns = ch4_fr.columns.to_series().apply(lambda x: "bge2_" + x)
    h2_fr.columns = h2_fr.columns.to_series().apply(lambda x: "bgr2_" + x)
    bgs_fr.columns = bgs_fr.columns.to_series().apply(lambda x: "bge1_" + x)
    result = pd.concat([eff_fr, ch4_fr, h2_fr, bgs_fr], axis=1)
    if not isinstance(results, pd.DataFrame):
        results = result
    else:
        results = pd.concat([results, result], axis=0)
    print(results)
    return results


def run_model(df):
    for idx in df.index:
        # Changine input variables
        inf_fr = df.loc[idx, "inf_fr"]
        temp = df.loc[idx, "temp"]
        hrt = df.loc[idx, "hrt"]

        # Fixed input variables
        n_stages = 2
        reactor_type = "FB"
        gas_extraction = "M"
        t_span = 200

        # set model
        sys = create_system(
            n_stages=n_stages,  # number of stages
            reactor_type=reactor_type,  # PB for packed bed, FB for fluidized bed, or UASB
            gas_extraction=gas_extraction,  # M for membrane gas extraction, V for vacuum extraction, P for passive venting
            Q=inf_fr,  # influent flowrate in m3/d
            T=temp,  # reactor temperature in degree C
            tot_HRT=hrt,  # total HRT in d
        )

        # run model
        sys.simulate(state_reset_hook="reset_cache", t_span=(0, t_span), method="BDF")

        # collect output data
        if idx == 0:
            output_data = None
        output_data = collect_results(case=sys, results=output_data)

    return output_data


def export_output_data(df, path=None, filename="output_data.csv"):
    df.to_csv(filename)
    print("The output data is readay")


if __name__ == "__main__":

    input_data = get_input_data(filename="input_data.csv")
    output_data = run_model(input_data)
    export_output_data(output_data)
