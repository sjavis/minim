.. _utils_api:

Utils API
=========

Helper utilities for common operations.
The MPI helper is documented separately in :ref:`communication_api`.


Printing
--------

The print functions can print any combination of arguments of different
types, including entire vectors.
``print`` outputs from the root processor only, so it is the usual
choice for program output, while ``printAll`` and ``printAllPlain``
output from every processor, prefixed with the processor rank or not.

.. doxygenfile:: print.h


Vector operations
-----------------

Element-wise arithmetic between vectors and scalars is provided by
operator overloads, e.g. ``a + b``, ``a * 2`` and ``a += b`` all work
as expected for vectors.

The ``vec`` namespace contains further vector operations in the style
of NumPy, such as ``vec::dotProduct``, ``vec::sum``, ``vec::norm`` and
``vec::slice``.

.. doxygenfile:: vec.h


Ranges
------

The range classes provide iterators over a multi-dimensional grid,
optionally with a halo.
``RangeX`` yields the grid indices at each position, while ``RangeI``
yields the flat index, e.g.::

    for (auto x : RangeX({nx, ny}, halo)) {
      // x is the vector of grid indices, including the halo
    }

.. doxygenclass:: minim::RangeX
   :members:

.. doxygenclass:: minim::RangeI
   :members:
