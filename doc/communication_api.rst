.. _communication_api:

Communication API
=================

The classes handling distribution over MPI processes.

The global :class:`~minim.mpi` object provides basic MPI operations over
all processors, such as broadcasts and sums, and is used to initialise
MPI at the start of a program.

Each :class:`~minim.State` holds a
:class:`~minim.Communicator`, created automatically from its potential,
which handles the distribution of the state's data over the processors
it uses. It is rarely necessary to interact with the communicator
directly, but it can be accessed via ``state.comm`` to make use of its
data assignment and MPI reduction functions.

See :doc:`parallelisation` for how the parallelisation works and how to
use it.


.. doxygenvariable:: minim::mpi

.. doxygenclass:: minim::Mpi
   :members:

.. doxygenclass:: minim::Communicator
   :members:


Communicator subclasses
-----------------------

The communicator type is chosen by the potential: grid potentials use
``CommGrid``, while element-based potentials use ``CommUnstructured``.

.. doxygenclass:: minim::CommGrid
   :members:

.. doxygenclass:: minim::CommGrid2
   :members:

.. doxygenclass:: minim::CommGrid3
   :members:

.. doxygenclass:: minim::CommUnstructured
   :members:
