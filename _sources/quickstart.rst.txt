.. _quickstart:

Quick Start
===========

Installation
------------

Download the repository and call ``make`` in the root directory to compile
the library::

    make

To use in a program, include the ``minim.h`` header file and compile with the
``-lminim`` flag::

    mpic++ -I$(MINIM)/include -L$(MINIM)/bin -lminim -DPARALLEL script.cpp -o run.exe

Refer to the ``examples`` folder for simple demonstrations of how to use the
library.


A first program
---------------

A minimisation consists of a potential, a state, and a minimiser.
The potential defines the energy, the state holds the coordinates, and the
minimiser performs the iterations.

.. code-block:: cpp

    #include "minim.h"

    using namespace minim;

    int main(int argc, char **argv) {
      mpi.init(&argc, &argv);

      // Define the system
      Lj3d potential = Lj3d();
      State state(potential, {0,0,0, 2,0,0, 1,1,0, 5,1,1});
      state.convergence = 1e-4;

      // Minimise
      Lbfgs min = Lbfgs();
      auto result = min.minimise(state);

      print("Complete after", min.iter, "iterations.");
      return 0;
    }

The available minimisers are :doc:`L-BFGS <algorithms/lbfgs>`,
:doc:`FIRE <algorithms/fire>`, :doc:`gradient descent
<algorithms/grad_descent>` and :doc:`simulated annealing
<algorithms/anneal>`.

See also the ``examples/`` directory for other simple example scripts.


Configuring a potential
-----------------------

Potentials are configured via a fluent interface of setter methods.
For example, a droplet on a surface using the
:class:`~minim.PhaseField` potential:

.. code-block:: cpp

    PhaseField pot;
    pot.setGridSize({nx, ny, 1});
    pot.setSolid([](int x, int y, int z){ return (y==0); });
    pot.setContactAngle(60);
    pot.setVolumeFixed(true);

    State state(pot, initCoords());

See the :doc:`potentials <potentials/index>` pages for the available
potentials and their parameters.


Logging during minimisation
---------------------------

Pass a log format to ``minimise`` to print information during the
minimisation.
The format is ``[fields]-[iter]``, where fields can include:

- ``e`` — the energy
- ``g`` — the gradient norm
- ``eg`` — both

For example, ``min.minimise(state, "eg-100")`` prints the energy and
gradient norm every 100 iterations.


Parallelisation
---------------

Minim can distribute a state over multiple MPI processes and share the
work of each processor over multiple OpenMP threads.
See :doc:`parallelisation` for how to run in parallel, including under
Slurm, and how the parallelisation is implemented.
