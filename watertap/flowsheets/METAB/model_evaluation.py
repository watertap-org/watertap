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

import pandas as pd

# third-party
try:
    import exposan
    from exposan.metab import create_system
except ImportError:
    exposan = None


__author__ = "Marcus Holly"


def get_input_data(filename=None):
    """
    Load input variable data from a CSV file or return a hardcoded default dataset.
    If no filename is provided, a small default dataset with three sample points
    is returned, covering influent flowrate, temperature, and hydraulic retention
    time.

    Args:
        filename (str, optional): Path to a CSV file containing input variable
            columns

    Returns:
        pd.DataFrame: DataFrame with one row per sample point
    """
    if filename == None:
        input_data = {
            "inf_fr": [5, 5, 5],
            "temp": [20, 25, 30],  # First column data
            "hrt": [12, 13, 14],  # Second column data
        }
        # Create the DataFrame
        df = pd.DataFrame(input_data)
    else:
        df = pd.read_csv(filename, header=0, skipinitialspace=True)

    return df


def get_eff_fr(case=None, df=None):
    """
    Extract effluent molar flowrates from a METAB system and append them to a DataFrame

    Args:
        case: An EXPOsan METAB system object with an ``outs`` attribute
        df: Existing DataFrame to append the new row
            to. If ``None``, a new DataFrame is created from the extracted
            flowrates.

    Returns:
        pd.DataFrame: DataFrame containing one row of effluent molar flowrates,
        with component names as columns
    """
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
        comp = [i.strip() for i in r.split("flow (kmol/hr):")[1].split("\n")]

        fr_dict = {}
        for element in comp:
            element_list = element.split()
            if not element_list:
                continue
            fr_dict[element_list[0]] = [float(element_list[-1])]
        if df is None:
            df = pd.DataFrame(fr_dict)
        else:
            df.loc[len(df)] = fr_dict
    return df


def get_ch4_fr(case=None, df=None):
    """
    Extract methane stream molar flowrates from a METAB system and append them to a DataFrame

    Args:
        case: An EXPOsan METAB system object with an ``outs`` attribute
        df: Existing DataFrame to append the new row
            to. If ``None``, a new DataFrame is created from the extracted
            flowrates.

    Returns:
        pd.DataFrame: DataFrame containing one row of CH4 stream molar
        flowrates, with component names as columns
    """
    if case is None:
        print(" The system is off")
    else:
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
        comp = [i.strip() for i in r.split("flow (kmol/hr):")[1].split("\n")]

        fr_dict = {}
        for element in comp:
            element_list = element.split()
            if not element_list:
                continue
            fr_dict[element_list[0]] = [float(element_list[-1])]
        if df is None:
            df = pd.DataFrame(fr_dict)
        else:
            df.loc[len(df)] = fr_dict
    return df


def get_h2_fr(case=None, df=None):
    """
    Extract hydrogen stream molar flowrates from a METAB system and append them to a DataFrame.

    Args:
        case: An EXPOsan METAB system object with an ``outs`` attribute
        df: Existing DataFrame to append the new row
            to.

    Returns:
        pd.DataFrame: DataFrame containing one row of H2 stream molar
        flowrates, with component names as columns
    """
    if case is None:
        print(" The system is off")
    else:
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
        comp = [i.strip() for i in r.split("flow (kmol/hr):")[1].split("\n")]

        fr_dict = {}
        for element in comp:
            element_list = element.split()
            if not element_list:
                continue
            fr_dict[element_list[0]] = [float(element_list[-1])]
        if df is None:
            df = pd.DataFrame(fr_dict)
        else:
            df.loc[len(df)] = fr_dict
    return df


def get_r1_ex_biogas_fr(case=None, df=None):
    """
    Extract reactor 1 extracted biogas molar flowrates from a METAB system and append them to a DataFrame.

    Args:
        case: An EXPOsan METAB system object with an ``outs`` attribute
        df: Existing DataFrame to append the new row
            to

    Returns:
        pd.DataFrame: DataFrame containing one row of reactor 1 biogas molar
        flowrates, with component names as columns
    """
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
        comp = [i.strip() for i in r.split("flow (kmol/hr):")[1].split("\n")]

        fr_dict = {}
        for element in comp:
            element_list = element.split()
            if not element_list:
                continue
            fr_dict[element_list[0]] = [float(element_list[-1])]
        if df is None:
            df = pd.DataFrame(fr_dict)
        else:
            df.loc[len(df)] = fr_dict
    return df


def get_mass_flowrate(case=None, df=None, stream=None):
    """
    Extract mass flowrates and volumetric flowrate from a stream and append them to a DataFrame.

    Reads component mass flowrates directly from a stream's state vector,
    using the stream's component list as column names. A ``Volumetric Flowrate``
    column is appended from the last element of the state vector. Columns where
    all values are zero are dropped from the result.

    Args:
        case: An EXPOsan METAB system object
        df: Existing DataFrame to append the new row
            to
        stream: An EXPOsan stream object with ``components`` and ``state``
            attributes

    Returns:
        pd.DataFrame: DataFrame with one row of mass flowrates and volumetric
        flowrate, with zero-only columns removed.
    """
    keys = str(stream.components).split("(")[1].rstrip(")").split(",")
    keys = [k.strip() for k in keys]
    values = [[float(x)] for x in list(stream.state[:-1])]
    fr_dict = dict(zip(keys, values))
    fr_dict["Volumetric Flowrate"] = [stream.state[-1]]
    print(fr_dict)
    if df is None:
        df = pd.DataFrame(fr_dict)
    else:
        df.loc[len(df)] = fr_dict

    df = df.loc[:, (df != 0).any(axis=0)]

    return df


def collect_results(case=None, results=None, mass=True):
    """
    Gathers flowrates from four output streams — effluent, methane biogas
    (reactor 2 extraction), hydrogen biogas (reactor 1 extraction), and reactor
    1 extracted biogas — and concatenates them into a single wide-format row.
    Column names are prefixed to distinguish streams:

    - ``eff_``  : effluent digestate (no prefix applied)
    - ``bge2_`` : biogas extraction from reactor 2 (CH4-rich)
    - ``bgr2_`` : biogas from reactor 2 hydrogen stream
    - ``bge1_`` : biogas extraction from reactor 1

    Args:
        case: An EXPOsan METAB system object whose ``outs`` streams will be
            read.
        results (pd.DataFrame, optional): Accumulator DataFrame from previous
            iterations.
        mass (bool, optional): If ``True``, use ``get_mass_flowrate`` for all
            streams. If ``False``, use the molar flowrate helpers
            (``get_eff_fr``, ``get_ch4_fr``, ``get_h2_fr``,
            ``get_r1_ex_biogas_fr``).

    Returns:
        pd.DataFrame: Updated accumulator with the current run's output
        appended as a new row
    """
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
    return results


def run_model(df):
    """
    Iterates over all rows of the input DataFrame, configures a two-stage
    fluidized-bed METAB system with membrane gas extraction for each set of
    conditions, simulates it using the BDF solver, and accumulates the output
    stream flowrates into a single results DataFrame.

    Fixed model parameters applied to every run:

    - ``n_stages``       : 2
    - ``reactor_type``   : ``"FB"`` (fluidized bed)
    - ``gas_extraction`` : ``"M"`` (membrane)
    - ``t_span``         : 200 (simulation time in days)

    Args:
        df (pd.DataFrame): Input DataFrame with one row per simulation run.
            Must contain columns ``inf_fr`` (influent flowrate, m³/d),
            ``temp`` (reactor temperature, °C), and ``hrt`` (total hydraulic
            retention time, d).

    Returns:
        pd.DataFrame: Accumulated output flowrates from all simulation runs,
        as returned by ``collect_results``. Returns an empty dict if
        ``exposan`` is not installed or no rows are processed.
    """
    output_data = {}
    for idx in df.index:
        # Changing input variables
        inf_fr = df.loc[idx, "inf_fr"]
        temp = df.loc[idx, "temp"]
        hrt = df.loc[idx, "hrt"]

        # Fixed input variables
        n_stages = 2
        reactor_type = "FB"
        gas_extraction = "M"
        t_span = 200

        if exposan is not None:
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
            sys.simulate(
                state_reset_hook="reset_cache", t_span=(0, t_span), method="BDF"
            )

            # collect output data
            if idx == 0:
                output_data = None
            output_data = collect_results(case=sys, results=output_data)

    return output_data


def export_output_data(df, filename=None):
    """
    Export a DataFrame of output results to a CSV file.

    Args:
        df: Output data to export
        filename (str, optional): Full path (including filename and extension)
            for the output CSV file. Passed directly to ``df.to_csv()``.
            Defaults to ``None``.

    Returns:
        None
    """
    df.to_csv(filename)
    print("The output data is ready")


if __name__ == "__main__":

    input_data = get_input_data(filename="input_data.csv")
    output_data = run_model(input_data)
    export_output_data(output_data)
