**********************************
Stability of Rock and Armour Units
**********************************

The stability of the breakwater is divided into three sections, the stability
of the armour layer, the toe and the substructure of the breakwater. With the
functions in substructure the required nominal diameter of the first underlayer
and filter layer can be computed.

Armour Layer
============
The armour layer is the outermost layer of the breakwater, the stability of
this layer is thus important for the overall stability of the breakwater. The
stability of this layer can be computed with either the Hudson formula or
the Van der Meer formula


Hudson
------

.. autofunction:: breakwater.core.stability.hudson

Van der Meer
------------

Three functions are implemented for the Van der Meer equations, one for deep
water conditions, one for shallow water conditions and a full functions. The
latter function takes into account the range in which the functions are
applicable. The applicability range is determined by the ratio of the water
depth over the significant wave height, see table 8.1 for these ranges.

.. figure:: _figures/vandermeer.png
   :align: center
   :alt: range of applicability of the Van der Meer equations for deep and
         shallow water

   Table 8.1: Range of applicability of the Van der Meer equations for deep
   and shallow water conditions (CIRIA, CUR, CETMEF, 2007, p. 579)

.. autofunction:: breakwater.core.stability.vandermeer

.. autofunction:: breakwater.core.stability.vandermeer_deep

.. autofunction:: breakwater.core.stability.vandermeer_shallow

Toe
===

.. autofunction:: breakwater.core.toe.toe_stability

.. autofunction:: breakwater.core.toe.toe_berm_stability

Substructure
============

.. autofunction:: breakwater.core.substructure.underlayer

.. autofunction:: breakwater.core.substructure.filter_layers

.. autofunction:: breakwater.core.substructure.layer_coefficient
