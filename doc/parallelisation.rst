.. _parallelisation:

Parallelisation
===============

Minim can be parallelised at two levels:

- **MPI** distributes a state across multiple processors, which may be on
  different nodes of a cluster. This is the primary form of
  parallelisation and is used to study large systems.
- **OpenMP** shares the work of a single processor between multiple
  threads on the same node. This is used to speed up a few key
  operations without any communication cost.

Both are enabled by default when the library is built with ``make``.


Running with MPI
----------------

Parallelisation is enabled at compile time with the ``-DPARALLEL`` flag,
which the ``Makefile`` includes by default, so no extra setup is needed.
Since MPI requires initialising, user programs must call
``mpi.init`` before creating any states:

.. code-block:: cpp

    int main(int argc, char **argv) {
      mpi.init(&argc, &argv);
      ...
    }

The number of mpi ranks and the current rank can be reference by
``mpi.size`` and ``mpi.rank``, respectively.

A program is then run in the usual way with ``mpirun``::

    mpirun -np 4 ./run.exe

Everything proceeds identically whether the program is run on one
processor or many: the energy and gradient functions always return the
total values over the whole system, and the coordinates can always be
read and written in full. The parallelisation is abstracted away by the
state's :class:`~minim.Communicator`, which handles all of the
communication between processors in the background.

The degrees of freedom of each :class:`~minim.State` are distributed
over the processors automatically when the state is created; there is no
need to manually split the data. All of the built-in minimisers and
potentials work in parallel without modification. By default a state is
spread over all of the MPI ranks, but a subset can also be used by
passing the ``ranks`` argument to its constructor, which allows
different states to be minimised simultaneously on different processors.

For grid potentials, the number of processors along each grid dimension
can be controlled with ``setCommArray``, e.g.::

    PhaseField pot;
    pot.setGridSize({nx, ny, 1});
    pot.setCommArray({2, 2, 1});

which splits the grid over a 2×2 array of processors (and so requires at
least four ranks). By default it will attempt to automatically set the
number along each dimension to make the sub-grids as cubic as possible.


Using OpenMP
------------

OpenMP is enabled at compile time with the ``-fopenmp`` flag, which the
``Makefile`` includes by default. It is used to parallelise some of the
more intensive computations that are local to each processor, such as
element loops in the phase field potentials, the vector operations used
by L-BFGS, and the dot products of grid data.

The number of threads is controlled in the usual way with the
``OMP_NUM_THREADS`` environment variable::

    OMP_NUM_THREADS=4 mpirun -np 4 ./run.exe

which runs four MPI ranks with four OpenMP threads each.

Running under Slurm
-------------------

On a cluster managed by Slurm, the three modes of running are:

**MPI only** — one MPI rank per core, no OpenMP threads. Use
``--ntasks`` to request the number of ranks and launch with ``srun``:

.. code-block:: bash

    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --ntasks=32
    #SBATCH --time=01:00:00

    srun ./run.exe

**OpenMP only** — a single MPI rank with several threads. Use
``--ntasks=1`` with ``--cpus-per-task`` for the thread count, and set
``OMP_NUM_THREADS`` to match:

.. code-block:: bash

    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --ntasks=1
    #SBATCH --cpus-per-task=16
    #SBATCH --time=01:00:00

    export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
    srun ./run.exe

**Hybrid MPI + OpenMP** — several ranks, each with several threads. The
total core count is ``--ntasks`` times ``--cpus-per-task``. Set
``OMP_NUM_THREADS`` to the number of threads *per rank*, and bind each
rank (and its threads) to its own cores so that threads belonging to
different ranks do not interfere:

.. code-block:: bash

    #!/bin/bash
    #SBATCH --nodes=2
    #SBATCH --ntasks=8
    #SBATCH --cpus-per-task=8
    #SBATCH --time=01:00:00

    export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
    export OMP_PROC_BIND=true
    export OMP_PLACES=cores
    srun --cpubind=cores ./run.exe

The same layouts apply when launching with ``mpirun`` instead of
``srun`` (e.g. ``mpirun -np 32 ./run.exe``); inside a Slurm job,
``mpirun`` will usually detect the allocated ranks automatically.

Beware of oversubscription: requesting more ranks or threads than the
allocated cores will slow the program down. In particular, leaving
``OMP_NUM_THREADS`` unset may cause each rank to spawn as many threads
as there are cores on the node, i.e. ``ntasks``-fold more than intended.


How the parallelisation is implemented
--------------------------------------

Data distribution
~~~~~~~~~~~~~~~~~

Each state holds a :class:`~minim.Communicator`, created automatically
from its potential, which defines how the degrees of freedom are divided
between the processors. Each processor is assigned one *block*: a set of
degrees of freedom specific to that processor. To compute the energy and
gradient of its block, a processor may also require degrees of freedom
owned by its neighbours, so it additionally stores a *halo* of degrees
of freedom from the neighbouring blocks. Together, the block plus halo
form the *proc* data held by that processor.

There are two types of communicator.

``CommGrid`` is used by potentials defined on a structured grid, where
the grid is split into a ``commArray`` of sub-grids, one per processor,
each surrounded by a halo of the given ``haloWidth``.

``CommUnstructured`` is used by element-based potentials such as
:class:`~minim.BarAndHinge`, where the degrees of freedom are split into
equal-sized blocks in the order they are given. For these, the energy
and gradient are broken into *elements*, each depending on only a few
degrees of freedom, and each element is assigned to the processor that
holds the most of its degrees of freedom. The remaining degrees of
freedom of that element then lie in the halo. Defining a potential in
terms of elements therefore allows the halo regions to be determined
automatically. For quantities that cannot be split into elements, such
as a system-wide volume constraint, a potential may instead define a
``blockEnergyGradient`` function that receives the communicator and
performs its own communication.

When creating a state, the communicator sets up the block and halo
layout, remaps the elements to local indices, and creates MPI derived
datatypes describing the halo and edge regions of the local data.

Communication
~~~~~~~~~~~~~

All communication is handled by the communicator. The key operations
are:

- ``communicate`` sends the edges of each block into the halo regions of
  the neighbouring processors, updating stale halo data (e.g. after
  randomly perturbing the coordinates).
- ``communicateAccumulate`` instead adds the values of a processor's
  halo onto the corresponding entries of the neighbours' blocks. This is
  how gradient contributions from elements that span two blocks are
  combined.
- ``gather`` and ``scatter`` convert between the global data and the
  local proc data.
- Reductions such as ``sum``, ``norm`` and ``dotProduct`` compute a
  result from the local blocks and sum it over all processors, so every
  processor receives the total.

Since halo entries are duplicated between processors, operations that
combine data across processors (like ``dotProduct``) only count each
block's entries once, excluding the halo.

The state's energy and gradient functions build on these operations.
The block functions compute a single processor's contribution without
any communication, so summing the block energies over the processors
gives the total. The plain ``energy`` / ``gradient`` functions do this
summing (and the equivalent gather for the gradient) internally, which
is why they can be used as if the program were serial. See the
:class:`~minim.State` documentation for the full set of functions and
their parallel behaviour.

Where MPI communication can be avoided, it is. For example, the
gradient returned by ``procGradient`` is correct over the whole proc
region (block plus halo), so the minimisers can take their steps using
only local data and only need to communicate when the next gradient is
required.

Hybrid MPI + OpenMP
~~~~~~~~~~~~~~~~~~~

OpenMP is used within a single processor to parallelise operations that
require no communication, such as loops over the local energy elements
and the vector arithmetic in L-BFGS. It is applied with ``#pragma omp
parallel for`` directives, so it requires no changes to user code. The
dot product over grid data is a special case: rather than excluding the
halo entries one by one, each processor sums over its interior and the
result is reduced over the processors with MPI.
