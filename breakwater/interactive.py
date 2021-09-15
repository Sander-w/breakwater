from matplotlib.widgets import Slider

from .utils._kwarg_validator import _RM_vkwargs, _C_vkwargs
from .app.main import BreakwaterDesign
from .utils.exceptions import InputError, NotSupportedError
from .utils.cost import _process_cost


def interactive_design(
        structure, type, LimitState, Grading, ArmourUnit=None, BermMaterial=None,
        cost=None, Soil=None, display_warnings=True):
    """ Interactive Design Application

    Application to interactively design breakwaters. All parameters can
    be specified in the application, together with varying parameters.
    On the Interactive Design Page it is then possible to change the
    value of a varying parameter with a slider. The design can then be
    updated and a cross-section of the new design is shown.

    .. note
       The maximum number of structures that can be specified is 2.

    Parameters
    ----------
    structure : {'RRM', 'CRM', 'RC', 'CC'}
        structure for which conceptual designs must be generated. RRM
        for a rubble mound with rock as armour layer, CRM for a rubble
        mound with concrete armour units as armour layer, RC for a
        vertical (composite) breakwater with rock as armour layer for
        the foundation and CC for a vertical (composite) breakwater with
        concrete armour units as armour layer for the foundation.
    type : {'Material', 'CO2'}
        Which type of cost analysis to do
    LimitState : :py:class:`LimitState` or list of :py:class:`LimitState`
        ULS, SLS or another limit state defined with
        :py:class:`LimitState`
    Grading : :py:class:`RockGrading`
        standard rock grading defined in the NEN-EN 13383-1 or a user
        defined rock grading. Required for all parameters
    ArmourUnit : obj, optional, default: None
        armour unit class which inherits from :py:class:`ConcreteArmour`,
        for instance :py:class:`Xbloc` or :py:class:`XblocPlus`. This
        argument is used for RRM.
    BermMaterial : obj, optional, default: None
        armour unit class which inherits from :py:class:`ConcreteArmour`,
        for instance :py:class:`Xbloc` or :py:class:`XblocPlus`. This
        argument is used for CC.
    Soil : :py:class:`Soil`, optional, default: None
        by default Soil is None, which means that the geotechnical checks
        are not performed. By specifying a Soil object, the geotechnical
        checks are automatically performed. Optional argument for RC and
        CC.
    cost : dict, optional, default: None
        by default no cost are computed, when the cost must be computed
        all relevant cost must be included. Optional keys of the dict
        are: core_price for the price of the core material, unit_price
        for the price of armour units,  concrete_price for the price of
        concrete, fill_price for the price of the fill material of the
        caisson. If not included in the prices it is also possible to
        specify transport_price for the price of transporting rocks,
        dry_dock for the rent of a dry dock, which is divided through
        the length of the breakwater.
    display_warnings : bool, optional, default: True
        if warnings must be displayed when designing
    """
    # convert the input of structure to a list
    if isinstance(structure, list):
        # must be a list so no change
        structure = structure
    elif isinstance(structure, str):
        # convert single input to list
        structure = [structure]

    # convert single LimitState to list if needed
    if isinstance(LimitState, list):
        LimitStates = LimitState
    else:
        LimitStates = [LimitState]

    # create dict of the python input for the application
    python_input = {
        'LimitState': LimitStates, 'Grading': Grading, 'Soil': Soil}

    # process python input
    if 'CRM' in structure:
        if ArmourUnit is None:
            raise InputError(f'ArmourUnit is a required argument for CRM')
        else:
            python_input['ArmourUnit'] = ArmourUnit

    if 'CC' in structure:
        if BermMaterial is None:
            raise InputError(f'BermMaterial is a required argument for CC')
        else:
            python_input['BermMaterial'] = BermMaterial

    # dict with the structures to design
    to_design = {
        'RM': {'num': 0, 'RRM': False, 'CRM': False},
        'C': {'num': 0, 'RC': False, 'CC': False}
        }

    # dict to store the vkwargs in
    vkwargs = {}

    # check structure type
    if 'RRM' in structure and 'CRM' in structure:
        vkwargs.update(_RM_vkwargs(type='both'))
        to_design['RM']['RRM'] = True
        to_design['RM']['CRM'] = True
        to_design['RM']['num'] = 2

    elif 'RRM' in structure:
        vkwargs.update(_RM_vkwargs(type='Rock'))
        to_design['RM']['RRM'] = True
        to_design['RM']['num'] = 1

    elif 'CRM' in structure:
        vkwargs.update(_RM_vkwargs(type='ArmourUnit'))
        to_design['RM']['CRM'] = True
        to_design['RM']['num'] = 1

    if 'RC' in structure and 'CC' in structure:
        vkwargs.update(_C_vkwargs(type='both'))
        to_design['C']['RC'] = True
        to_design['C']['CC'] = True
        to_design['C']['num'] = 2

    elif 'RC' in structure:
        vkwargs.update(_C_vkwargs(type='Rock'))
        to_design['C']['RC'] = True
        to_design['C']['num'] = 1

    elif 'CC' in structure:
        vkwargs.update(_C_vkwargs(type='ArmourUnit'))
        to_design['C']['CC'] = True
        to_design['C']['num'] = 1

    # check number of structures
    # due to fixed lay out the maximum number is 2
    total_num = to_design['RM']['num'] + to_design['C']['num']
    if total_num > 2:
        # to many structures
        raise NotSupportedError(
            ('Too many structures have been specified, currently only 2 '
             'structures can be designed at the same time'))

    # delete python input from parameters as these are given in Python
    for parameter in python_input.keys():
        if parameter in vkwargs.keys():
            del vkwargs[parameter]

    # process the cost
    cost = _process_cost(structure, type, cost, Grading, validate=False)

    # start the app
    app = BreakwaterDesign(
        vkwargs, python_input, to_design, cost, display_warnings)
    app.mainloop()
