********
Material
********

The most important construction materials for breakwaters are quarry stone and
concrete, as this is the material used in the armour layer. Both types can be
used for designing breakwaters, the rock grading is defined in the first
section, whereafter the concrete armour units are defined.

Rock
====

.. autoclass:: breakwater.material.RockGrading
   :members:

Concrete Armour Units
=====================

The availability of rock is limited. Furthermore, natural blocks weighing more
than 10 tons are very rare. However, sometimes a block larger than 10 tons is
needed or a breakwater is constructed in an area with limited availability of
rock. In these cases artificial blocks, made of concrete, are used.

Several types of concrete armour units have been developed over the years.
Currently only Xbloc and XblocPlus from Delta Marine Consultants have been
incorporated in the package. However, it is possible to define another type
of concrete armour unit by using :py:class:`ConcreteArmour`.

General
-------

.. autoclass:: breakwater.material.ConcreteArmour
   :members:

Predefined types
----------------

.. autoclass:: breakwater.material.Xbloc
   :inherited-members:
   :members:


.. autoclass:: breakwater.material.XblocPlus
   :inherited-members:
   :members:

User defined Armour Unit
------------------------

As mentioned at the start of this section it is possible to define another type
of Armour Unit by using the general class
:obj:`bw.ConcreteArmour <breakwater.material.ConcreteArmour>`. In this
subsection an example is given how to define a custom armour unit, which can
be used in :obj:`bw.ConcreteRubbleMound <breakwater.rubble.ConcreteRubbleMound>`
to design a breakwater. The following code must be used to define the custom
armour unit:

.. code-block:: python

   # define class of the armour unit which inherits from bw.ConcreteArmour
   class CustomArmourUnit(bw.ConcreteArmour):
       """ Custom Armour Unit class """

       # define the init method, note that it is not required to give a
       # keyword argument.
       def __init__(self, rho=2400):
           """ See help(CustomArmourUnit) for more info """
           # define the armour units in the init method
           units = {
               1: {'D': 1, 'h': 1, 'Vc': 1},
               2: {'D': 2, 'h': 2, 'Vc': 2},
               3: {'D': 3, 'h': 3, 'Vc': 3},
               4: {'D': 4, 'h': 4, 'Vc': 4},
           }

           # call the init method of bw.ConcreteArmour with super and pass
           # the following arguments: kd factor for Hudson formula, name of
           # the armour unit for overtopping and rho for the density
           super().__init__(kd=10, units=units, name='CustomArmourUnit', rho=rho)

       # this is an optional step. It is possible to define a method
       # to compute a correction factor to be applied on the computed
       # required nominal diameter. Method must allow for kwargs input
       # as several parameters are passed to the method, it is not
       # necessary to use all of the passed parameters
       def correction_factor(self, h, Hs, Rc, **kwargs):
           # it is then possible to build logic in order to determine the
           # correction factor, for example:
           # define list to store more correction factors
           correction = []

           # check for water depth wave height relation
           if h > 2.5*Hs:
               correction.append(1.3)
           else:
               pass

           # check for low crested structure
           if Rc/Hs < 1:
               correction.append(1.5)
           else:
               pass

           # check if any correction factors have been added
           if any(correction):
               # return maximum value
               return max(correction)
           else:
               # return 1, in other words, no correction factor
               # make sure a correction_factor is always returned
               return 1

.. warning::
   When using one of the design classes, the supported armour units are limited
   to the armour units for which an influence factor for the permeability and
   roughness of the slope has been determined.

.. note::
   It is not required to add a method to compute a correction factor, the design
   classes do raise an exception if there is no method to compute a correction
   factor. However, in case you want to define a method to compute a correction
   factor the following kwargs can be used to determine the correction factor:

   +-----------------+-------+-----------------------------------------------------------------+
   | Argument        | Type  | Description                                                     |
   +=================+=======+=================================================================+
   | h               | float | water depth                                                     |
   +-----------------+-------+-----------------------------------------------------------------+
   | Hs              | float | significant wave height                                         |
   +-----------------+-------+-----------------------------------------------------------------+
   | Rc              | float | crest freeboard                                                 |
   +-----------------+-------+-----------------------------------------------------------------+
   | occurrence_hs   | bool  | | frequent occurrence of the near-design wave height during the |
   |                 |       | | during the lifetime of the structure                          |
   +-----------------+-------+-----------------------------------------------------------------+
   | slope           | float | slope of the armour layer in radians                            |
   +-----------------+-------+-----------------------------------------------------------------+
   | slope_foreshore | float | slope of the foreshore in radians                               |
   +-----------------+-------+-----------------------------------------------------------------+
   | permeability    | str   | permeability of the core {permeable, low, impermeable}          |
   +-----------------+-------+-----------------------------------------------------------------+
   | Dn              | float | nominal diameter of the armour units in the armour layer        |
   +-----------------+-------+-----------------------------------------------------------------+
   | layers          | int   | number of layers in the armour layer                            |
   +-----------------+-------+-----------------------------------------------------------------+
   | B               | float | crest width                                                     |
   +-----------------+-------+-----------------------------------------------------------------+
   | beta            | float | | angle between direction of wave approach and a line normal to |
   |                 |       | | the breakwater in degrees                                     |
   +-----------------+-------+-----------------------------------------------------------------+
