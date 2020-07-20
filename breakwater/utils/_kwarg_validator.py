import numpy as np

from ..material import RockGrading, ConcreteArmour
from ..conditions import LimitState
from ..core.soil import Soil
from .exceptions import InputError

def _RM_vkwargs(type):
    """ Valid kwargs table for a rubble mound breakwater

    Construct a dictionary with all the valid kwargs for a rubble mound
    breakwater.

    Parameters
    ----------
    type : {'Rock', 'ArmourUnit', 'both'}
        type of rubble mound breakwater, R for rock and C for concrete
        armour units.

    Returns
    -------
    vkwargs : dict
        dictionary with all the valid kwargs for a rubble mound
        breakwater
    """
    # vkwargs for both types of rubble mound breakwaters
    vkwargs = {
        'slope': {'Default': None,
                  'Validator': lambda value: isinstance(value, tuple)
                                             and len(value) == 2,
                  'Correct': 'tuple with length 2',
                  'Constant': False,
                  'Required': True,
                  'unit': '[-]'},

        'slope_foreshore': {'Default': None,
                            'Validator': lambda value: isinstance(value, tuple)
                                                       and len(value) == 2,
                            'Correct': 'tuple with length 2',
                            'Constant': True,
                            'Required': True,
                            'unit': '[-]'},

        'rho_w' : {'Default': None,
                   'Validator': lambda value: isinstance(value, float)
                                              or isinstance(value, int),
                   'Correct': 'int or float',
                   'Constant': True,
                   'Required': True,
                   'unit': '[kg/m³]'},

        'B' : {'Default': None,
               'Validator': lambda value: isinstance(value, float)
                                          or isinstance(value, int),
               'Correct': 'int or float',
               'Constant': False,
               'Required': True,
               'unit': '[m]'},

        'LimitState' : {'Default': None,
                        'Validator': lambda value: isinstance(value, LimitState),
                        'Correct': 'LimitState',
                        'Constant': True,
                        'Required': True,
                        'unit': '[-]'},

        'Grading' : {'Default': None,
                     'Validator': lambda value: isinstance(value, RockGrading),
                     'Correct': 'RockGrading',
                     'Constant': True,
                     'Required': True,
                     'unit': '[-]'},

        'Dn50_core' : {'Default': None,
                       'Validator': lambda value: isinstance(value, int)
                                                  or isinstance(value, float),
                       'Correct': 'int or float',
                       'Constant': False,
                       'Required': True,
                       'unit': '[m]'},

        'safety' : {'Default': 1,
                    'Validator': lambda value: isinstance(value, int)
                                               or isinstance(value, float),
                    'Correct': 'int or float',
                    'Constant': True,
                    'Required': False,
                    'unit': '[-]'},

        'slope_toe' : {'Default': (2,3),
                       'Validator': lambda value: isinstance(value, tuple)
                                                  and len(value) == 2,
                       'Correct': 'tuple with length 2',
                       'Constant': False,
                       'Required': False,
                       'unit': '[-]'},

        'B_toe' : {'Default': None,
                   'Validator': lambda value: isinstance(value, int)
                                              or isinstance(value, float),
                   'Correct': 'int or float',
                   'Constant': False,
                   'Required': False,
                   'unit': '[m]'},

        'layers_underlayer' : {'Default': 2,
                               'Validator': lambda value: isinstance(value, int),
                               'Correct': 'int',
                               'Constant': True,
                               'Required': False,
                               'unit': '[-]'},

        'beta' : {'Default': 0,
                  'Validator': lambda value: isinstance(value, float)
                                             or isinstance(value, int),
                  'Correct': 'int or float',
                  'Constant': True,
                  'Required': False,
                  'unit': '[degrees]'},

        'Soil': {'Default': None,
                 'Validator': lambda value: isinstance(value, Soil),
                 'Correct': 'Soil',
                 'Constant': True,
                 'Required': False,
                 'unit': '[-]'},

        'phi': {'Default': 40,
                'Validator': lambda value: isinstance(value, float)
                                           or isinstance(value, int),
                'Correct': 'int or float',
                'Constant': True,
                'Required': False,
                'unit': '[degrees]'},

    }

    # only vkwargs if type of structure is rock rubble mound
    vkwargs_RRM = {
        'N' : {'Default': None,
               'Validator': lambda value: isinstance(value, int),
               'Correct': 'int',
               'Constant': True,
               'Required': True,
               'unit': '[-]'},

        'layers_rock' : {'Default': 2,
                         'Validator': lambda value: isinstance(value, int),
                         'Correct': 'int',
                         'Constant': True,
                         'Required': False,
                         'unit': '[-]'},

        'vdm' : {'Default': 'max',
                 'Validator': lambda value: isinstance(value, str)
                                            and value in ['min', 'max', 'avg'],
                 'Correct': 'str and {min, max, avg}',
                 'Constant': True,
                 'Required': False,
                 'unit': '[-]'},
    }

    # only vkwargs if type of structure is concrete rubble mound
    vkwargs_CRM = {
        'ArmourUnit' : {'Default': None,
                        'Validator': lambda value: isinstance(value, ConcreteArmour),
                        'Correct': 'ConcreteArmour or inherit from ConcreteArmour',
                        'Constant': True,
                        'Required': True,
                        'unit': '[-]'},

        'layers_units' : {'Default': 1,
                          'Validator': lambda value: isinstance(value, int),
                          'Correct': 'int',
                          'Constant': True,
                          'Required': False,
                          'unit': '[-]'},

        'filter_rule' : {'Default': None,
                         'Validator': lambda value: value in ['Rock', 'Xbloc', 'XblocPlus'],
                         'Correct': 'str and {Rock, Xbloc or XblocPlus}',
                         'Constant': True,
                         'Required': False,
                         'unit': '[-]'},

    }

    if type == 'Rock' or type == 'both':
        vkwargs.update(vkwargs_RRM)

    if type == 'ArmourUnit' or type == 'both':
        vkwargs.update(vkwargs_CRM)

    return vkwargs

def _C_vkwargs(type):
    """ Valid kwargs table for a vertical (composite) breakwater

    Construct a dictionary with all the valid kwargs for a vertical
    (composite) breakwater

    Parameters
    ----------
    type : {'Rock', 'ArmourUnit', 'both'}
        type of rubble mound breakwater, R for rock and C for concrete
        armour units.

    Returns
    -------
    vkwargs : dict
        dictionary with all the valid kwargs for a vertical (composite)
        breakwater
    """
    vkwargs = {
        'Pc': {'Default': None,
               'Validator': lambda value: isinstance(value, float)
                                          or isinstance(value, int),
               'Correct': 'int or float',
               'Constant': False,
               'Required': True,
               'unit': '[-]'},

        'rho_c': {'Default': None,
                   'Validator': lambda value: isinstance(value, float)
                                              or isinstance(value, int),
                   'Correct': 'int or float',
                   'Constant': True,
                   'Required': True,
                   'unit': '[kg/m³]'},

        'rho_fill': {'Default': None,
                     'Validator': lambda value: isinstance(value, float)
                                                or isinstance(value, int),
                     'Correct': 'int or float',
                     'Constant': True,
                     'Required': True,
                     'unit': '[kg/m³]'},

        'rho_w' : {'Default': None,
                   'Validator': lambda value: isinstance(value, float)
                                              or isinstance(value, int),
                   'Correct': 'int or float',
                   'Constant': True,
                   'Required': True,
                   'unit': '[kg/m³]'},

        'Bm': {'Default': None,
               'Validator': lambda value: isinstance(value, float)
                                       or isinstance(value, int),
               'Correct': 'int or float',
               'Constant': False,
               'Required': True,
               'unit': '[m]'},

        'hb': {'Default': None,
               'Validator': lambda value: isinstance(value, float)
                                       or isinstance(value, int),
               'Correct': 'int or float',
               'Constant': False,
               'Required': True,
               'unit': '[m]'},

        'LimitState' : {'Default': None,
                        'Validator': lambda value: isinstance(value, LimitState),
                        'Correct': 'LimitState',
                        'Constant': True,
                        'Required': True,
                        'unit': '[-]'},

        'slope_foreshore': {'Default': None,
                            'Validator': lambda value: isinstance(value, tuple)
                                                       and len(value) == 2,
                            'Correct': 'tuple with length 2',
                            'Constant': True,
                            'Required': True,
                            'unit': '[-]'},

        'mu': {'Default': None,
               'Validator': lambda value: isinstance(value, float),
               'Correct': 'float',
               'Constant': True,
               'Required': True,
               'unit': '[-]'},

        'safety' : {'Default': 1,
                    'Validator': lambda value: isinstance(value, int)
                                               or isinstance(value, float),
                    'Correct': 'int or float',
                    'Constant': True,
                    'Required': False,
                    'unit': '[-]'},

        'SF_sliding': {'Default': 1.2,
                       'Validator': lambda value: isinstance(value, int)
                                                  or isinstance(value, float),
                       'Correct': 'int or float',
                       'Constant': True,
                       'Required': False,
                       'unit': '[-]'},

        'SF_turning': {'Default': 1.2,
                       'Validator': lambda value: isinstance(value, int)
                                                  or isinstance(value, float),
                       'Correct': 'int or float',
                       'Constant': True,
                       'Required': False,
                       'unit': '[-]'},

        'beta': {'Default': 0,
                 'Validator': lambda value: isinstance(value, int)
                                            or isinstance(value, float),
                 'Correct': 'int or float',
                 'Constant': True,
                 'Required': False,
                 'unit': '[degrees]'},

        'slope_foundation' : {'Default': (2,3),
                              'Validator': lambda value: isinstance(value, tuple)
                                                         and len(value) == 2,
                              'Correct': 'tuple with length 2',
                              'Constant': False,
                              'Required': False,
                              'unit': '[-]'},

        'lambda_': {'Default': [1,1,1],
                    'Validator': lambda value: isinstance(value, list)
                                               and len(value) == 3,
                    'Correct': 'list with length 3',
                    'Constant': True,
                    'Required': False,
                    'unit': '[-]'},

        'filter_rule' : {'Default': None,
                         'Validator': lambda value: value in ['Rock', 'Xbloc', 'XblocPlus'],
                         'Correct': 'str and {Rock, Xbloc or XblocPlus}',
                         'Constant': True,
                         'Required': False,
                         'unit': '[-]'},

        'Grading' : {'Default': None,
                     'Validator': lambda value: isinstance(value, RockGrading),
                     'Correct': 'RockGrading',
                     'Constant': True,
                     'Required': True,
                     'unit': '[-]'},

        'pe_max' : {'Default': 500,
                    'Validator': lambda value: isinstance(value, int)
                                               or isinstance(value, float),
                    'Correct': 'int or float',
                    'Constant': True,
                    'Required': False,
                    'unit': '[kPa]'},

        'Soil': {'Default': None,
                 'Validator': lambda value: isinstance(value, Soil),
                 'Correct': 'Soil',
                 'Constant': True,
                 'Required': False,
                 'unit': '[-]'},
    }

    # only vkwargs if type of structure is rock foundation for the caisson
    vkwargs_RC = {

        'layers_rock' : {'Default': 2,
                         'Validator': lambda value: isinstance(value, int),
                         'Correct': 'int',
                         'Constant': True,
                         'Required': False,
                         'unit': '[-]'},
    }

    # only vkwargs if type of structure is armour unit foundation for caisson
    vkwargs_CC = {

        'layers_units' : {'Default': 1,
                    'Validator': lambda value: isinstance(value, int),
                    'Correct': 'int',
                    'Constant': True,
                    'Required': False,
                    'unit': '[-]'},

        'BermMaterial': {'Default': None,
                         'Validator': lambda value: isinstance(value, RockGrading)
                                                    or isinstance(value, ConcreteArmour),
                         'Correct': 'RockGrading or ConcreteArmour',
                         'Constant': True,
                         'Required': True,
                         'unit': '[-]'},

    }

    if type == 'Rock' or type == 'both':
        vkwargs.update(vkwargs_RC)

    if type == 'ArmourUnit' or type == 'both':
        vkwargs.update(vkwargs_CC)

    return vkwargs

def _validate_kwargs(params, vkwargs):
    """ Validate if all required input is given

    Parameters
    ---------
    params : dict
        dictionary with the processed kwargs
    vkwargs : dict
        dictionary with the valid kwargs, default values and a validator
    """
    for key, value in params.items():
        required = vkwargs[key]['Required']
        if required and value is None:
            raise InputError(
                f'{key} is a required parameter for the specified structure')
        else:
            continue

def _process_varying_kwarg(key, value, correct):
    """" Process varying kwarg

    Parameters
    ----------
    key : str
        name of the kwarg
    value : any
        value of the kwarg
    correct : str
        correct format of the parameter if not varying

    Returns
    -------
    value : any
        processed value
    """
    # set bool for converting slope
    convert_slope_to_tuple = False

    if isinstance(value, tuple) and len(value) == 3:
        # parameter is a tuple and has length of 3
        # check if value[0] and value[1] are tuples
        if isinstance(value[0], tuple) and isinstance(value[1], tuple):
            # parameter is a slope thus preprocessing is needed
            # tuples must be converted to floats (angles)
            lower = value[0][0]/value[0][1]
            upper = value[1][0]/value[1][1]

            # set new format of value for varying slopes
            value = (lower, upper, value[2])

            # set bool to convert values back to tuples
            convert_slope_to_tuple = True

        # check if format is correct (min, max, num)
        if value[1] > value[0] and isinstance(value[2], int):
            # make linspace of value and return the linspace
            varying_values = np.linspace(value[0], value[1], value[2])

            # check if varying value must be converted
            if convert_slope_to_tuple:
                converted_slopes = []

                # convert values to tuple
                for varying_val in varying_values:
                    converted_slopes.append((1, 1/varying_val))

                # convert the list of converted slopes to a ndarray
                dt = np.dtype('int, float')
                varying_values = np.array(converted_slopes, dtype=dt)

            # return the varying values
            return varying_values
        else:
            # incorrect format of min, max, num
            raise InputError(
                (f'{key} is allowed to vary but tuple must be formatted as '
                '(min, max, num)'))
    else:
        # incorrect type of length
        given_type = type(value).__name__
        given_len = len(value)
        raise InputError(
            (f'{key} is allowed to vary but must be tuple with length 3, is '
             f'{given_type} with length {given_len}. Don\'t want a varying '
             f'parameter? Specify argument as {correct}'))

def _process_kwargs(kwargs, vkwargs):
    """ Process the kwargs

    Check if the input kwarg is in the valid kwargs table, verify that
    a keyword is valid and that the value has the correct type.
    Furthermore, validate if all kwargs for the type of structure are
    given.

    Parameters
    ----------
    kwargs : dict
        dictionary with kwargs and their values
    vkwargs : dict
        dictionary with the valid kwargs, default values and a validator

    Returns
    -------
    dict
        dictionary with the kwargs and their values, and default values
    """
    # initialize params from vkwargs
    params = {}
    for key, value in vkwargs.items():
        params[key] = value['Default']

    # validate the kwargs, and replace the value in params
    # if the kwarg is a valid kwarg
    for key in kwargs.keys():
        if key not in vkwargs:
            # invalid kwarg of current structure, do not raise error
            # kwarg might be for another type of structure
            continue
        else:
            value = kwargs[key]
            valid = vkwargs[key]['Validator'](value)
            if valid:
                # argument has valid format, thus add to params
                params[key] = value
            else:
                # argument is invalid, thus raise error
                incorrect = type(value).__name__
                correct = vkwargs[key]['Correct']

                # check if parameter is allowed to vary
                if not vkwargs[key]['Required'] and value is None:
                    params[key] = value
                elif not vkwargs[key]['Constant']:
                    # parameter is allowed to vary
                    processed_value = _process_varying_kwarg(key, value, correct)

                    params[key] = processed_value

                else:
                    # parameter not allowed to vary thus incorrect input
                    raise InputError(
                        (f'{incorrect} is invalid for kwarg {key}, must be '
                         f'{correct}'))

    # validate if all required parameters are given
    _validate_kwargs(params, vkwargs)

    return params
