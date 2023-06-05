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
