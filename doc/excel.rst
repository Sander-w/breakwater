*******************
Reading Excel files
*******************

Besides specifying the input in Python it is also possible to load a parametric
design or material from an Excel file. The Excel input file can be made by the
user in Excel, but it advised to use
:obj:`bw.generate_excel <breakwater.utils.input_generator.generate_excel>` to
generate the required Excel input file.

Creating Excel input file
=========================

.. autofunction:: breakwater.utils.input_generator.generate_excel

Parametric Design
=================

.. autofunction:: breakwater.design.read_configurations

Material
========

.. autofunction:: breakwater.material.read_grading

.. autofunction:: breakwater.material.read_units
