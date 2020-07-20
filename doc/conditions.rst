***********
Limit State
***********

The design conditions, such as the design wave height and the damage number, are
described in a limit state. A limit state is defined by standard
NEN-EN 1990:2002 (2002), as the state beyond which the structure no longer
fulfils the relevant design criteria. Two of the most common limit states are:

- The Ultimate Limit State (ULS), which is the state associated with collapse
  of the breakwater or with other similar forms of failure
- The Serviceability Limit State (SLS), which is the state that corresponds to
  the conditions beyond which specified service requirements for the breakwater
  or part of the breakwater are no longer met.

Define Limit State
==================

.. autoclass:: breakwater.conditions.LimitState
   :members:

Compute wave heights
====================

Design criteria for coastal structures require, among other parameters, a
certain type of characteristic wave height, for instance the significant wave
height. In deep water the behaviour of the waves can be determined with the
characteristic of the Rayleigh distribution. However, in shallow water the waves
are no longer Rayleigh distributed, but a designer still wants to determine the
significant wave height.

Two methods have been implemented to determine the nearshore wave
characteristics. The first is Battjes and Groenendijk, which determines the wave
characteristics from the output of a wave energy model
(Battjes and Groenendijk, 2000). Secondly, the formulation by Goda (2000) which
is based on experimental data.

Battjes and Groenendijk
-----------------------

.. autoclass:: breakwater.core.battjes.BattjesGroenendijk
   :members:

Goda wave heights
-----------------

.. autofunction:: breakwater.core.goda.goda_wave_heights

.. autofunction:: breakwater.utils.wave.shoaling_coefficient
