********
Tutorial
********

This tutorial guides you through some typical use cases. See table 1 for the
discussed tutorials.

Table 1: List of tutorials and description

+-------------------------------+--------------------+------------------------------------------------------+
| Tutorial                      | File               | Description                                          |
+===============================+====================+======================================================+
| Rubble Mound with Rock        | :obj:`example1.py` | | Using BattjesGroenendijk to compute the wave       |
|                               |                    | | heights and design with single :obj:`LimitState`   |
+-------------------------------+--------------------+------------------------------------------------------+
| Rubble Mound with Xbloc       | :obj:`example1.py` | | Using BattjesGroenendijk to compute the wave       |
|                               |                    | | heights and design with single :obj:`LimitState`   |
+-------------------------------+--------------------+------------------------------------------------------+
| Vertical with Rock            | :obj:`example2.py` | | Using goda_wave_heights to compute wave heights    |
|                               |                    | | and design with two :obj:`LimitState` objects      |
+-------------------------------+--------------------+------------------------------------------------------+
| Parametric design with Python | :obj:`example3.py` | | Parametric design for Rubble Mound breakwater      |
|                               |                    | | with Rock (RRM) and Xbloc (CRM)                    |
+-------------------------------+--------------------+------------------------------------------------------+
| Parametric design from Excel  | :obj:`example4.py` | | Generate Excel input file and make a parametric    |
|                               |                    | | design from Excel an input file                    |
+-------------------------------+--------------------+------------------------------------------------------+

Rubble Mound Design
===================

Two designs are elaborated on. First a conventional rubble mound breakwater with
rock as armour layer and secondly, one with armour units as armour layer. This
example can also be found in :obj:`example1.py`

Rock
----

First the package must be imported.

.. code-block:: python

   import breakwater as bw

Importing :obj:`breakwater` like this will automatically import the most commonly
used classes and functions. These are, all design classes, all materials, the
limit state and the functions to compute wave heights. Which is also the next
step. Compute a missing wave height, H2% in this case, and define a limit state.

.. code-block:: python

   battjes = bw.BattjesGroenendijk(Hm0=4.4, h=15, slope_foreshore=(1,100))
   H2_per = battjes.get_Hp(0.02)

   ULS = bw.LimitState(
       h=15, Hs=4.5, Hm0=4.4, H2_per=H2_per, Tp=9.4, Tm=8.8, T_m_min_1=9.7,
       Sd=5, Nod=2, q=20, label='ULS')

The next step is to define a material for the armour layer, in this case we will
first design a breakwater with rock as armour layer. We use the standard rock
grading defined in the standard NEN-EN 13383 (2002) with a density of 2650
:math:`kg/m^3`, which is the default value.

.. code-block:: python

   NEN = bw.RockGrading(rho=2650)

We are now ready to design a breakwater.

.. code-block:: python

   RRM = bw.RockRubbleMound(
       slope=(2,3), slope_foreshore=(1,100), rho_w=1025, B=5.5, N=2100,
       LimitState=ULS, Grading=NEN, Dn50_core=0.4)

The breakwater is designed, let’s  first inspect if any warnings were encountered
during the design process.

.. code-block:: python

   RRM.print_logger(level='warnings')

No warnings were encountered, so we are ready to explore the design. Let’s first
see at how the breakwater looks, and the details of the generated variants.

.. code-block:: python

   RRM.plot('all')
   RRM.print_variant('all')

We can see that two variants have been designed, with a different filter layer
but the same underlayer. Finally, let’s inspect the validity ranges as the
equations of van der Meer for the armour layer and toe are experimental formulae

.. code-block:: python

   RRM.check_validity()

From the table the validity ranges, used value and if the parameter is in range
can be read. It is up to the user to assess if the used values are reasonable.

Armour Units
------------

The first steps are the same as when rock was the armour layer, as can be seen
from the code.

.. code-block:: python

   import breakwater as bw

   battjes = bw.BattjesGroenendijk(Hm0=4.4, h=15, slope_foreshore=(1,100))
   H2_per = battjes.get_Hp(0.02)

   ULS = bw.LimitState(
       h=15, Hs=4.5, Hm0=4.4, H2_per=H2_per, Tp=9.4, Tm=8.8, T_m_min_1=9.7,
       Sd=5, Nod=2, q=20, label='ULS')

   NEN = bw.RockGrading(rho=2650)

It is still required to define a rock grading as the underlayer of the
breakwater is made out of rock. We will now define the armour units. It is
possible to define a custom armour layer with :py:class:`bw.ConcreteArmour`,
but in the example we will use a predefined one.

.. code-block:: python

   xbloc = bw.Xbloc()

We are now ready to design a breakwater with Xbloc as armour layer.

.. code-block:: python

   CRM = bw.ConcreteRubbleMound(
       slope=(2,3), slope_foreshore=(1,100), B=5.5, rho_w=1025, LimitState=ULS,
       ArmourUnit=xbloc, Grading=NEN, Dn50_core=0.4)

The breakwater is designed, again let’s first inspect if any warnings were
encountered during the design process.

.. code-block:: python

   CRM.print_logger(level='warnings')

Again, no warnings were encountered during the design, so we are ready to
explore the generated design. So let’s again plot and print all variants.

.. code-block:: python

   CRM.plot('all')
   CRM.print_variant('all')

These values can now be compared to the concept with rock as armour layer, or
further designed as the geotechnical stability is not designed by the tool.

Monolithic Breakwater Design
============================

In this section it is explained how a vertical breakwater with rock as
armour layer of the foundation can be designed. This example can also be
found in :obj:`example2.py`

First the package must be imported.

.. code-block:: python

   import breakwater as bw

Importing :obj:`breakwater` like this will automatically import the most commonly
used classes and functions. These are, all design classes, all materials, the
limit state and the functions to compute wave heights. Which is also the next
step. We will transform the deep water wave height to a design wave height,
Hmax, with the empirical formulae derived by Goda (2000).

.. code-block:: python

   H13_ULS, Hmax_ULS = bw.goda_wave_heights(
       h=15.1, d=12, Ho=5.3, T=9.4, slope_foreshore=(1,100))

   H13_SLS, Hmax_SLS = bw.goda_wave_heights(
       h=12.1, d=9, Ho=3.3, T=7.9, slope_foreshore=(1,100))

The next step is to define the limit state functions. In this example we will
define two limit states, the ultimate limit state (ULS) and a serviceability
limit state (SLS).

.. code-block:: python

   ULS = bw.LimitState(
       h=15.1, H13=H13_ULS, Hmax=Hmax_ULS, T13=9.4, q=30, label='ULS')
   SLS = bw.LimitState(
       h=12.1, H13=H13_SLS, Hmax=Hmax_SLS, T13=7.9, q=15, label='SLS')

   ULS.transform_periods(0.5)
   SLS.transform_periods(0.5)

As you can see we also used the method :py:meth:`transform_periods` to
transform the missing wave periods, especially :math:`T_{m-1.0}` as it is used
to compute the required crest freeboard. Since the assumption that we are in
deep water is not valid the method will display a LimitStateWarning, the better
practice is to derive the wave periods from a model like SWAN.

The next step is to define a material for the armour layer off the foundation,
in this example we will use the default rock grading of the NEN-EN 13383 (2002).

.. code-block:: python

   NEN = bw.RockGrading()

We are now ready to design a vertical breakwater!

.. code-block:: python

   RC = bw.Caisson(
       Pc=0.2, rho_c=2400, rho_fill=1600, rho_w=1000, Bm=8, hb=2, layers=2,
       BermMaterial=NEN, LimitState=[ULS, SLS], slope_foreshore=(1,100), mu=0.5,
       beta=15)

Let us first inspect the logger, we will print the info level as well in this
example because this allows us to see which overtopping formula is used.
Furthermore, we can also see if overturning or sliding was normative for the
computation of the required width.

.. code-block:: python

   RC.print_logger(level='info')

We now know some general info, so let’s explore the design. Let’s
see at how the breakwater looks, and the details of the generated variants.

.. code-block:: python

  RC.plot('all')
  RC.print_variant('all')

Design Automation
=================

In this section a design automation tool is used. First the design is made
fully in Python and secondly the design is made from and Excel input file. The
values used in this section are based on the ones used in `Rubble Mound Design`_.

With Python
-----------

*This example can also be found in* :obj:`example3.py`

Similar to the example for a rubble mound breakwater we again begin with by
defining the wave heights and materials to used.

.. code-block:: python

   battjes = bw.BattjesGroenendijk(Hm0=4.4, h=15, slope_foreshore=(1,100))
   H2_per = battjes.get_Hp(0.02)

   ULS = bw.LimitState(
       h=15, Hs=4.5, Hm0=4.4, H2_per=H2_per, Tp=9.4, Tm=8.8, T_m_min_1=9.7,
       Sd=5, Nod=2, q=20, label='ULS')

   NEN = bw.RockGrading(rho=2650)
   xbloc = bw.Xbloc()

Now that the hydraulic conditions and materials are defined, we can use the
:py:class:`Configurations` class. This means that multiple configurations of
the specified breakwater type(s) are designed, in the example the parameters
allowed to vary are: the slope, the width of the breakwater (B) and den nominal
diameter of the core (Dn50_core).

.. code-block:: python

   configs = bw.Configurations(
       structure=['RRM', 'CRM'], LimitState=ULS, rho_w=1025,
       slope_foreshore=(1,100), Grading=NEN, slope=((1,3), (3,4), 4), B=(5, 8, 4),
       Dn50_core=(0.2, 0.4, 3), N=2100, ArmourUnit=xbloc)

The first step is to see if, and which, warnings have been encountered during
the design process.

.. code-block:: python

   configs.show_warnings()

In the table we can see that a lot of warnings have been encountered, the most
are related to the hydraulic input. Instead of H13 we specified Hs in the
LimitState, but the design looks for H13. This value is not found and thus
Hs is used instead, which is fine as H13 is a definition of Hs. However, the
tool does inform the user about this behaviour.

As 96 concepts are a lot of concepts to filter by yourself, the concepts can be
exported to Design Explorer 2. This online tool allows the user to visually
filter the concepts, see
:obj:`to_design_explorer <breakwater.design.Configurations.to_design_explorer>`.
We export the generated concepts with the following parameters: the slope, the
class of the armour layer, the width of the breakwater and the crest freeboard.

.. code-block:: python

   configs.to_design_explorer(params=['slope', 'class armour', 'B', 'Rc'])

The last step is to save the generated concepts to a .breakwaters file. This
allows the user to access all the concepts at a latter moment without designing
them again.

.. code-block:: python

   configs.to_breakwaters('example3')

The concepts can be reloaded with
:obj:`bw.read_breakwaters <breakwater.design.read_breakwaters>`.

With Excel input file
---------------------

As mentioned in the introduction of this section a design can also be made
from an Excel input file. Note that this script is available in
:obj:`example4.py`.

.. warning::
   While it is possible to make your own Excel input file, you are advised to
   use :obj:`bw.generate_excel <breakwater.utils.input_generator.generate_excel>`
   to generate the required excel file.

The first step is to create the required Excel input file.

.. code-block:: python

   import breakwater as bw

   bw.generate_excel('config input.xlsx')

Now you can give the required input in an Excel file, so give your input in
the Excel file. We can now import the Excel file and make a design for multiple
configurations.

.. code-block:: python

   configs = bw.read_configurations(
       'config input.xlsx', structure=['RRM', 'CRM'], kd=16, name='Xbloc')

This function returns a :obj:`bw.Configurations <breakwater.design.Configurations>`
object. Therefore, it is now possible to show the warnings, export the
concepts to the design explorer, just as we have done in the previous example.

.. note::
   In this example we specified all classes of Xbloc in the Excel file. However,
   as these armour units are already defined in :obj:`breakwater` the better
   practice is to specify them in
   :obj:`bw.read_configurations <breakwater.design.read_configurations>` with
   the keyword argument ArmourUnits. This allows :obj:`breakwater` to compute
   the correction factor, which results in a better design.

.. warning::
   It is currently not possible to use the Excel input file with multiple
   limit states. Want to design with multiple limit states? See the previous
   subsection, `With Python`_
