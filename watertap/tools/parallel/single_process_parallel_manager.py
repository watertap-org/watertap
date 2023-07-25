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

from watertap.tools.parallel.parallel_manager import build_and_execute, ParallelManager
from watertap.tools.parallel.results import LocalResults


class SingleProcessParallelManager(ParallelManager):
    def __init__(self, **kwargs):
        self.results = None

    def is_root_process(self):
        return True

    def get_rank(self):
        return self.ROOT_PROCESS_RANK

    def number_of_worker_processes(self):
        return 1

    def sync_with_peers(self):
        pass

    def sync_data_with_peers(self, data):
        pass

    def scatter(
        self,
        do_build,
        do_build_kwargs,
        do_execute,
        all_parameters,
    ):

        self.results = LocalResults(
            self.ROOT_PROCESS_RANK,
            all_parameters,
            build_and_execute(do_build, do_build_kwargs, do_execute, all_parameters),
        )

    def gather(self):
        return [self.results]

    def results_from_local_tree(self, results):
        return results
