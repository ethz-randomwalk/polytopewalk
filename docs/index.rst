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

The Python bindings are created using `pybind11` and are documented here.

.. toctree::
   :maxdepth: 2
   :caption: Python Modules

   modules


C++ API
========

The C++ API is documented using Doxygen.

.. doxygenindex::


Specific C++ Classes (Examples)
====

.. doxygenclass:: RandomWalk
   :members:
   :protected-members:

.. doxygenclass:: BallWalk
   :members:
   :protected-members: