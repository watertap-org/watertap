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

import numpy as np
import pyomo.environ as pyo
import sys

from scipy.interpolate import griddata
from enum import Enum, auto
from abc import abstractmethod, ABC
from idaes.core.util import get_solver

from idaes.surrogate.pysmo import sampling
from pyomo.common.collections import ComponentSet
from pyomo.common.tee import capture_output

np.set_printoptions(linewidth=200)

# class _WriterBase(ABC):
#     def __init__(self):
#         pass
#
#     @abstractmethod
#     def strip_extension(self, file_name, extension):
#         if file_name.lower().endswith(extension):
#             return file_name[: -len(extension)], extension
#         else:
#             return file_name, None
#
#     @abstractmethod
#     def process_results_filename(self, results_file_name):
#         # Get the directory path
#         dirname = os.path.dirname(results_file_name)
#         # Get the file name without the extension
#         known_extensions = [".h5", ".csv"]
#         for ext in known_extensions:
#             fname_no_ext, extension = _strip_extension(results_file_name, ext)
#             if extension is not None:
#                 break
#
#         return dirname, fname_no_ext, extension
#
# class CSVWriter(_WriterBase):
#     def __init__(self):
#         pass

class ParameterSweepWriter:

    def __init__(self,
        results_file_name = None,
        debugging_data_dir = None,
        comm = None,
        write_h5 = False,
        write_csv = False,
        interpolate_nan_outputs = False):

    # ================================================================

    def _strip_extension(self, file_name, extension):
        if file_name.lower().endswith(extension):
            return file_name[: -len(extension)], extension
        else:
            return file_name, None

    # ================================================================

    def _process_results_filename(self, results_file_name):
        # Get the directory path
        dirname = os.path.dirname(results_file_name)
        # Get the file name without the extension
        known_extensions = [".h5", ".csv"]
        for ext in known_extensions:
            fname_no_ext, extension = _strip_extension(results_file_name, ext)
            if extension is not None:
                break

        return dirname, fname_no_ext, extension

    # ================================================================

    def save_results(self, sweep_params, local_values, global_values, local_results_dict,
        global_results_dict, global_results_arr, results_file_name, write_csv,
        write_h5, debugging_data_dir, comm, rank, num_procs, interpolate_nan_outputs,
    ):

        if results_file_name is not None:
            dirname, fname_no_ext, extension = _process_results_filename(results_file_name)
            if rank == 0:
                # Create the directories for the results and/or for debugging
                if extension == ".h5" and write_h5 is False:
                    warnings.warn(
                        "An H5 results filename was provided. Outputs will be created",
                        UserWarning,
                    )
                    write_h5 = True
                elif extension == ".csv" and write_csv is False:
                    warnings.warn(
                        "A CSV results filename was provided. Outputs will be created",
                        UserWarning,
                    )
                    write_csv = True

                if (write_h5 or write_csv) and dirname != "":
                    os.makedirs(dirname, exist_ok=True)
                elif not write_h5 and not write_csv:
                    warnings.warn(
                        "A results filename was provided but neither options to write H5 or csv was selected. No file will be written.",
                        UserWarning,
                    )

                if debugging_data_dir is not None:
                    os.makedirs(debugging_data_dir, exist_ok=True)
        else:
            dirname = None
            fname_no_ext = None

        if num_procs > 1:  # pragma: no cover
            comm.Barrier()

        # Handle values in the debugging data_directory
        if debugging_data_dir is not None:
            _write_debug_data(
                sweep_params,
                local_values,
                local_results_dict,
                debugging_data_dir,
                rank,
                write_h5,
                write_csv,
            )

        global_save_data = _write_to_csv(
            sweep_params,
            global_values,
            global_results_dict,
            global_results_arr,
            rank,
            write_csv,
            dirname,
            fname_no_ext,
            interpolate_nan_outputs,
        )

        if rank == 0 and write_h5:
            # Save the data of output dictionary
            _write_outputs(global_results_dict, dirname, fname_no_ext, txt_options="keys")

        return global_save_data

    def _write_debug_data(self, sweep_params, local_values, local_results_dict,
        debugging_data_dir, rank, write_h5, write_csv):

        if write_h5:
            fname_h5 = f"local_results_{rank:03}.h5"
            _write_output_to_h5(local_results_dict, debugging_data_dir, fname_h5)
        if write_csv:
            fname_csv = f"local_results_{rank:03}.csv"

            data_header = ",".join(itertools.chain(sweep_params))
            local_results = np.zeros(
                (np.shape(local_values)[0], len(local_results_dict["outputs"])),
                dtype=np.float64,
            )
            for i, (key, item) in enumerate(local_results_dict["outputs"].items()):
                data_header = ",".join([data_header, key])
                local_results[:, i] = item["value"][:]

            local_save_data = np.hstack((local_values, local_results))

            # Save the local data
            np.savetxt(
                os.path.join(debugging_data_dir, fname_csv),
                local_save_data,
                header=data_header,
                delimiter=", ",
                fmt="%.6e",
            )


    def _write_outputs(self, output_dict, output_directory, fname_no_ext, txt_options="metadata"):

        h5_results_file = fname_no_ext + ".h5"

        _write_output_to_h5(output_dict, output_directory, h5_results_file)

        # We will also create a companion txt file by default which contains
        # the metadata of the h5 file in a user readable format.
        txt_fname = fname_no_ext + ".txt"
        txt_fpath = os.path.join(output_directory, txt_fname)
        if "solve_successful" in output_dict.keys():
            output_dict.pop("solve_successful")
        if txt_options == "metadata":
            my_dict = copy.deepcopy(output_dict)
            for key, value in my_dict.items():
                for subkey, subvalue in value.items():
                    subvalue.pop("value")
        elif txt_options == "keys":
            my_dict = {}
            for key, value in output_dict.items():
                my_dict[key] = list(value.keys())
        else:
            my_dict = output_dict

        with open(txt_fpath, "w") as log_file:
            pprint.pprint(my_dict, log_file)

    def _write_output_to_h5(self, output_dict, output_directory, fname):

        fpath = os.path.join(output_directory, fname)
        f = h5py.File(fpath, "w")
        for key, item in output_dict.items():
            grp = f.create_group(key)
            if key != "solve_successful":
                for subkey, subitem in item.items():
                    subgrp = grp.create_group(subkey)
                    for subsubkey, subsubitem in subitem.items():
                        if subsubkey[0] != "_":
                            if subsubkey == "lower bound" and subsubitem is None:
                                subgrp.create_dataset(subsubkey, data=np.finfo("d").min)
                            elif subsubkey == "upper bound" and subsubitem is None:
                                subgrp.create_dataset(subsubkey, data=np.finfo("d").max)
                            else:
                                subgrp.create_dataset(subsubkey, data=subsubitem)
            elif key == "solve_successful":
                grp.create_dataset(key, data=output_dict[key])

        f.close()
