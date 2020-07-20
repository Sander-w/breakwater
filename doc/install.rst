************
Installation
************

This chapter explains the dependencies and the different ways to install the
package.

Dependencies
============

For the package to work some prerequisite packages must be installed on you
system. It might be that you have already installed these packages, since it
are quite common packages.

- `Python`_ : Version 3.6 or higher

- `NumPy`_ : The fundamental scientific programming package, it
  provides a multidimensional array type and many useful functions for
  numerical analysis.

- `SciPy`_ : Library of algorithms for mathematics, science and engineering.

- `Matplotlib`_ : Matplotlib is a comprehensive library for creating static,
  animated, and interactive visualizations in Python.

- `Pandas`_ : A fast, powerful, flexible and easy to use open source data
  analysis and manipulation tool

- `tabulate`_ : Pretty-print tabular data

The following are optional dependencies and are only required for some features:

- `XlsxWriter`_ : a module that can be used to write text, numbers, formulas
  and hyperlinks to multiple worksheets in an Excel 2007+ XLSX file. Required
  dependency for
  :obj:`bw.generate_excel <breakwater.utils.input_generator.generate_excel>`

- `Basemap`_ : The matplotlib basemap toolkit is a library for plotting 2D data
  on maps in Python. Required dependency for
  :obj:`breakwater.database.BreakwaterDatabase <breakwater.database.database.BreakwaterDatabase>`

.. _`Python`: https://www.python.org/

.. _`NumPy`: https://numpy.org/doc/1.18/

.. _`SciPy`: http://www.scipy.org/

.. _`Matplotlib`: https://matplotlib.org/

.. _`Pandas`: https://pandas.pydata.org/docs/

.. _`tabulate`: https://pypi.org/project/tabulate/

.. _`XlsxWriter`: https://xlsxwriter.readthedocs.io/

.. _`Basemap`: https://matplotlib.org/basemap/users/installing.html

Installation
============

:obj:`breakwater` can be installed from PyPI or with the source code from the
GitHub repository.

From PyPI
---------

The latest release is available at the `Python package index`_

.. code-block:: python

   pip install breakwater

From GitHub
-----------

Download or clone the source code from the `GitHub`_ repository. Then move into
the directory with the source code and run:

.. code-block:: python

  python setup.py install


.. _`GitHub`: https://github.com/Sander-w/breakwater/
.. _`Python package index`: https://pypi.org/


Bugs and Feature Requests
=========================

Problems with the installation, bugs in the code or feature request can be
reported at the `issue tracker`_ of the GitHub repository. Comments and
questions are always welcome and can be send to `pybreakwater@gmail.com`_.


.. _`issue tracker`: https://github.com/Sander-w/breakwater/issues
.. _`pybreakwater@gmail.com`: pybreakwater@gmail.com
