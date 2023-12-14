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
    #. Reduce : Collects the data from all parallel workers and performs a mathematical operation. The current behavior sums the values from all workers.

Usage
--------

End User
~~~~~~~~

If the user is simply invoking the parameter sweep tool, they only need to
specify parallel backend by using the keyword argument ``parallel_back_end``
when constructing a parameter sweep object. Valid keyword values are

    #. "MPI" : For use on distributed memory systems such as a cluster or supercomputer. We have run the parameter sweep on 100 compute nodes with 10400 total workers.
    #. "ConcurrentFutures" : For use on shared memory systems (e.g. your workstation) with flow sheets that initialize quickly.
    #. "MultiProcessing" : For use on shared memory systems with flow sheets that take a long time to initialize.
    #. "RayIo" : For distributed memory systems where MPI is not stable, or hybrid computing.

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
class and defines the abstract methods. The methods can be found in the
class documentation below. 

* :mod:`watertap.tools.parallel.parallel_manager`
* :mod:`watertap.tools.parallel.parallel_manager_factory`
* :mod:`watertap.tools.parallel.concurrent_futures_parallel_manager`
* :mod:`watertap.tools.parallel.mpi_parallel_manager`
* :mod:`watertap.tools.parallel.multiprocessing_parallel_manager`
* :mod:`watertap.tools.parallel.ray_io_parallel_manager`
* :mod:`watertap.tools.parallel.single_process_parallel_manager`


Adding Features to the Parallel Manager
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An advanced user may add features as necessary to the parallel manager. We ask
that they follow appropriate code development practices and include testing and
documentation. A lot of features are currently missing from the parallel 
manager that limits its wider use. We encourage developers to make PRs to the
parallel manager.