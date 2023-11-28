Parallel Manager
====================

.. index::
    pair: watertap.tools.parallel.parallel_manager;parallel_manager

.. currentmodule:: watertap.tools.parallel.parallel_manager

Overview
--------

The parallel manager is a high level wrapper designed for the parameter sweep
tool that allows the use of different parallel paradigms for running
parameter sweep in parallel, for example, on a supercomputer, using a unified
API. Currently, supported parallel backends include:

    #. MPI
    #. Multiprocessing
    #. Concurrent Futures
    #. Ray Core

Parallel manager was implemented to support different workflows of WaterTAP
users and developers while reducing the burden to develop software. In 
particular, it supports collective operations using the various paradigms.
These operations are listed below:

    #. Barrier : Waits until all parallel workers arrive at the same point in code
    #. Gather : Gathers an array or list on one or all parallel workers.
    #. Scatter : Divides up the array or list across all the parallel workers
    #. Broadcast : Sends the array or list to all parallel workers
    #. Reduce : Collects the data from all parallel worker and performs a mathematical operation. The current behavior sums the values from all workers.

Usage
--------

End User
~~~~~~~~

If the user is simply invoking the parameter sweep tool, they only need to
specify parallel backend by using the keyword argument ``parallel_back_end``
when constructing a parameter sweep object. Valid keyword values are

    #. "MPI"
    #. "ConcurrentFutures"
    #. "MultiProcessing"
    #. "RayIo"

If any other keyword is provided, the parallel manager reverts to serial. Note
that, when using MPI, the end user must still use the ``mpiexec`` or ``mpirun``
command when executing the parameter sweep from the command line.

Adding Parallel Manager to Your Code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to use the parallel manager in your code, you will need to adapt the
following statement. 

.. code:: python

    from watertap.tools.parallel.parallel_manager_factory import create_parallel_manager
    parallel_manager = create_parallel_manager(
            parallel_manager_class=None,
            number_of_subprocesses=2,
            parallel_back_end="ConcurrentFutures",
        )

The above code return 1 of 5 parallel manager classes. There are currently 5 parallel manager classes:

    #. ``SingleProcessParallelManager`` : To enable serial implementation
    #. ``MPIParallelManager`` : For MPI
    #. ``ConcurrentFuturesParallelManager`` : For python concurrent futures
    #. ``MultiprocessingParallelManager`` : For more advanced python multiprocessing
    #. ``RayIoParallelManager`` : For using Ray

Alternatively, the user can specify the parameter sweep class directly, for example

.. code:: python

    from watertap.tools.parallel.single_process_parallel_manager import SingleProcessParallelManager
    s_parallel_manager = create_parallel_manager(
            parallel_manager_class=SingleProcessParallelManager,
            number_of_subprocesses=1,
        )

Each of the parallel manager classes inherits from the base ``ParallelManager``
class and defines the abstract methods. The methods include:

* ``is_root_process`` : Return whether the current process should be considered as the root of the parallel group it is a part of.
* ``get_rank`` : Get the process's rank number in its parallel peer group.
* ``number_of_worker_processes`` : Return how many total processes are running local simulations.
* ``sync_with_peers`` : Implement a synchronization point with any peer processes.
* ``sync_array_with_peers``: Synchronize an array with all peers processes. The data parameter is either

    - (if root) the source of truth that will be broadcast to all processes
    - (if not root) an empty buffer that will hold the data sent from root once the function returns
* ``sync_pyobject_with_peers`` : Synchronize a python object with all peer processes. The ``obj`` parameter is either

    - (if root) the source of truth that will be broadcast to all processes
    - (if not root) ignored

    It is different from ``sync_array_with_peers`` in that it returns the synced object rather than
    using an existing buffer.
* ``combine_data_with_peers`` : Combine the data from each peer into a list. The return value will be a list of all the data elements provided by each process that calls this function. With multiple processes, this must:

    - act as a synchronization point
    - return the list in the same order on each process
* ``gather_arrays_to_root`` : Gather the data in the send buffer to the root process. Parameters are:

    - ``sendbuf``: the data to be sent from this process
    - ``recvbuf_spec``: a tuple of (receive buffer, list of sizes) where the list of sizes is how much data will come from each process. Ignored if not the root process.
* ``sum_values_and_sync`` : Sum a list of values from each process and sync the result to all processes. Parameters are:

    - ``sendbuf``: an array of values the local process is contributing to the sum
    - ``recvbuf``: an array of a single value that will hold the summed result (same on each process) when the function returns.

    .. note::
        This is a specific case of a global reduce (with ``sum()`` as the reducing function).
        If more than ``sum()`` is needed in the future, this function should be extended to receive
        an Op argument.

* ``scatter`` : Scatter the specified execution out, as defined by the implementation's parallelism, for a list of parameters. Arguments are

    - ``do_build``: a function that builds the arguments necessary for the execution function. It is expected to return a list that will be exploded and passed into the do_execute function as arguments.
    - ``do_build_kwargs``: a dictionary of keyword arguments for the do_build function
    - ``do_execute``: the execution function. It is expected to take in the list of local parameters as its first argument. Any arguments after that should match the list returned by the ``do_build`` function.
    - ``all_parameters``: a list of all parameters to be run. It is included so that different implementations of the parallel manager can make decisions about splitting and parallelization.

* ``gather`` : Gather the results of the computation that was kicked off via a previous scatter. Returns a list of ``LocalResults``, representing the results for each process. Each result will be the return value of the do_execute function from ``scatter()`` for one process.

* ``results_from_local_tree`` : Given a list of ``LocalResults`` objects, return a sublist of the ones the current process is responsible for.

The developer may use these methods to enable parallelism in their code using the parallel manager.


Adding Features to the Parallel Manager
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An advanced user may add features as necessary to the parallel manager. We ask
that they follow appropriate code development practices and include testing and
documentation. A lot of features are currently missing from the parallel 
manager that limits its wider use. We encourage developers to make PRs to the
parallel manager.