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


class LocalResults:
    """
    Class representing the results of one process's run of an optimization routine.
    Parameters:
    - process_number: a unique number identifying the process in the group
    - parameters: a list containing an entry for each parameter group run by the process
    - results: the results of the optimization routine run for all parameter groups in the list
    """

    def __init__(self, process_number, parameters, results):
        self.process_number = process_number
        self.parameters = parameters
        self.results = results
