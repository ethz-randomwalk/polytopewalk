.. polytopewalk documentation master file, created by
   sphinx-quickstart on Tue May  9 22:50:40 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to polytopewalk's documentation!
=======================================

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   python_api
   cpp_api
   support


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Python API
===========

The Python bindings are created using `pybind11` and are documented here. This section contains the Python API documentation for the PolytopeWalk library.

.. toctree::
   :maxdepth: 4
   :caption: Python Modules

   py_init
   py_dense_walk
   py_sparse_walk
   py_utils


C++ API
========

The C++ API is documented using Doxygen. 
This section provides an overview of the C++ API.

Here is a list of important classes in the C++ API:

* `RandomWalk`
   * `BallWalk`
   * `HitAndRun`
   * `BarrierWalk`
      * `DikinWalk`
      * `VaidyaWalk`
      * `JohnWalk`
      * `DikinLSWalk`
* `SparseRandomWalk`
   * `SparseBallWalk`
   * `SparseHitAndRun`
   * `SparseBarrierWalk`
      * `SparseDikinWalk`
      * `SparseVaidyaWalk`
      * `SparseJohnWalk`
      * `SparseDikinLSWalk`


For more detailed documentation on these classes, see the following pages.

.. toctree::
   :maxdepth: 4
   :caption: C++ Modules

   cpp_init
   cpp_dense_walk
   cpp_sparse_walk
   cpp_utils