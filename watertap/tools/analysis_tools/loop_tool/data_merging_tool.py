#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################


import h5py
import numpy as np
import time
import os
import datetime
from datetime import datetime

__author__ = "Alexander V. Dudchenko (SLAC)"


def merge_data_into_file(
    file_name,
    backup_file_name,
    directory,
    expected_solved_values=None,
    min_solve_values=None,
    force_rerun=None,
):
    """
    This function checks if there is a sim file, and a back up file.
    if there is a sim file, it backs it up, and checks if full solution set already exsts and copies it over
    other wise it creates a new group in actual file, and adds new solved data set to it
    """
    run_sweep = True
    if os.path.isfile(file_name) == False:
        create_h5_file(file_name)
    h5file = h5py.File(file_name, "a")
    try:
        # check if there is a back up
        if isinstance(backup_file_name, str):
            f_old_solutions = h5py.File(backup_file_name, "r")
            solved_values = sum(
                np.array(
                    f_old_solutions[directory]["solve_successful"]["solve_successful"][
                        ()
                    ]
                )
            )
        else:
            solved_values = 0
        if force_rerun == False:
            print("Forced to not run")
            run_sweep = False
        elif force_rerun == None:
            if min_solve_values != None:
                if min_solve_values > solved_values:
                    run_sweep = False
            elif expected_solved_values == solved_values:
                run_sweep = False

        if run_sweep:
            h5file.create_group(directory)
        elif os.path.isfile(backup_file_name):
            h5file.copy(f_old_solutions[directory], directory)
            f_old_solutions.close()
    except (KeyError, ValueError, FileNotFoundError):
        try:
            h5file.create_group(directory)
        except ValueError:
            pass
    h5file.close()
    return run_sweep


def create_backup_file(file_name, backup_name, h5_dir):
    """used to created file and back up file"""

    if backup_name is None and os.path.isfile(file_name):
        h5file = h5py.File(file_name, "r")
        print(h5_dir, h5file.get(h5_dir))
        if h5_dir in h5file:
            h5file.close()
            date = datetime.now().strftime("%d_%m-%H_%M_%S")
            backup_name = file_name + "_{}_{}".format(date, ".bak")
            if os.path.isfile(backup_name) == False:
                os.rename(file_name, backup_name)
        else:
            h5file.close()
    return backup_name


def create_h5_file(file_name):
    """used to created h5file"""
    if os.path.isfile(file_name) == False:
        for i in range(60):
            try:
                f = h5py.File(file_name, "w")
                f.close()
                break
            except:
                print("could not creat h5 file")
                time.sleep(0.01)  # Waiting to see if file is free to create again
