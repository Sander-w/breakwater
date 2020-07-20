**********************
Geotechnical stability
**********************

This chapter includes the implemented geotechnical equations. These computations
are automatically executed by the design classes when a :py:class:`Soil` is
specified. For the :py:class:`Caisson` class the :py:meth:`brinch_hansen` method
is used to compute the bearing capacity of the soil. In :py:class:`RockRubbleMound`
and :py:class:`ConcreteRubbleMound` a slip circle analysis is performed with
:py:class:`Bishop`. The required length of the scour protection is computed
for all design classes.

Define Soil
===========

.. autoclass:: breakwater.core.soil.Soil
   :members:

Slip Circle
===========

.. autoclass:: breakwater.core.bishop.Bishop
   :members:

.. autoclass:: breakwater.core.bishop.SlipCircle
   :members:

Scour
=====

.. autofunction:: breakwater.core.scour.scour_protection
