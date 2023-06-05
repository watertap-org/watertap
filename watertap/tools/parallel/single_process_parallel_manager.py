from watertap.tools.parallel.parallel_manager import ParallelManager, run_sweep
from watertap.tools.parallel.results import LocalResults


class SingleProcessParallelManager(ParallelManager):
    def __init__(self):
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

    def scatter(self, param_sweep_obj, common_params, all_parameter_combos):
        self.results = LocalResults(
            self.ROOT_PROCESS_RANK,
            all_parameter_combos,
            run_sweep(param_sweep_obj, common_params, all_parameter_combos),
        )

    def gather(self):
        return [self.results]

    def results_from_local_tree(self, results):
        return results
