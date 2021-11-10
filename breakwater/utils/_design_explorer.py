import numpy as np

from .exceptions import InputError

def _RM_vparams():
    """ Parameters of RM that can exported to the design explorer """
    vparams = {
        'B': {'param': 'B',
              'key': 'df'},

        'Dn50_core': {'param': 'Dn50_core',
                      'key': 'df'},

        'B_toe': {'param': '_B_toe',
                  'key': 'attr',
                  'nested': False},

        'Rc': {'param': 'Rc',
               'key': 'attr',
               'nested': False},

        'slope': {'param': '_input_arguments',
                  'key': 'attr',
                  'nested': True},

        'slope_toe': {'param': '_input_arguments',
                      'key': 'attr',
                      'nested': True},

        'computed Dn50 armour': {'param': 'computed Dn50',
                                 'key': 'armour'},

        'class armour': {'param': 'class',
                         'key': 'armour'},

        'class Dn50 armour': {'param': 'class Dn50',
                              'key': 'armour'},

        'computed Dn50 underlayer': {'param': 'computed Dn50',
                                     'key': 'underlayer'},

        'class underlayer': {'param': 'rock class',
                             'key': 'underlayer'},

        'class Dn50 underlayer': {'param': 'class Dn50',
                                  'key': 'underlayer'},

        'computed Dn50 filter': {'param': 'computed Dn50',
                                     'key': 'filter layer'},

        'class filter': {'param': 'rock class',
                         'key': 'filter layer'},

        'class Dn50 filter': {'param': 'class Dn50',
                              'key': 'filter layer'},

    }

    return vparams

def _C_vparams():
    """ Parameters of C that can be exported to the design explorer """
    vparams = {
        'Pc': {'param': 'Pc',
               'key': 'caisson'},

        'hb': {'param': 'hb',
               'key': 'caisson'},

        'h_acc': {'param': 'h_acc',
                  'key': 'caisson'},

        'Rc': {'param': 'Rc',
               'key': 'caisson'},

        'B': {'param': 'B',
               'key': 'caisson'},

        'Bm': {'param': 'Bm',
               'key': 'caisson'},

        'slope_foundation': {'param': '_input_arguments',
                             'key': 'attr',
                             'nested': True},

        'computed Dn50 armour': {'param': 'computed Dn50',
                                 'key': 'armour'},

        'class armour': {'param': 'class',
                         'key': 'armour'},

        'class Dn50 armour': {'param': 'class Dn50',
                              'key': 'armour'},

        'computed Dn50 underlayer': {'param': 'computed Dn50',
                                     'key': 'underlayer'},

        'class underlayer': {'param': 'rock class',
                             'key': 'underlayer'},

        'class Dn50 underlayer': {'param': 'class Dn50',
                                  'key': 'underlayer'},

        'UC': {'param': 'geotechnical',
               'key': 'attr',
               'nested': False}
    }

    return vparams

def _DE_params(
        args, variant, row, concept, structure, slopes, change_CRM_class):
    """ Get the parameters for the design explorer

    Parameters
    ----------
    args : list
        list with arguments
    variant : dict
        a variant of a concept of a type of breakwater
    row : pd.Series
        current row of the df
    concept : obj
        a breakwater concept
    structure : {'RRM', 'CRM', 'CRMR', 'RC', 'CC'}
        type of structure for export
    slopes : {angles, tuples}, optional, default: angles
        how the slopes must be exported. tuples will export the
        slope as (V,H) and angle as an angle in degrees
    change_CRM_class : bool
        True if the class of the armour layer of a CRM breakwater must
        be changed to a string, False if it must remain a float

    Returns
    -------
    dict
        dictionary with the parameter and value

    Raises
    ------
    InputError
        if the arg is not a valid argument for the design explorer
    """
    # set output dict
    to_explorer = {}

    # get the valid parameters
    if structure == 'RRM' or structure == 'CRM' or structure == 'CRMR':
        vparams = _RM_vparams()

    elif structure == 'RC' or structure == 'CC':
        vparams = _C_vparams()

    else:
        raise NotImplementedError(f'{structure} is not supported')

    for arg in args:
        if arg in vparams:
            # get parameter and key for variant (from structure)
            param = vparams[arg]['param']
            key = vparams[arg]['key']

            # get the desired value for design explorer
            if key == 'df':
                # value is in the df
                value = row[param]

            elif key == 'attr':
                # value is attribute of the concept
                # check if nested key
                if vparams[arg]['nested']:
                    # nested
                    value = getattr(concept, param)[arg]
                    param = arg
                else:
                    # not nested
                    # check if param is UC
                    if param == 'geotechnical':
                        value = getattr(concept, param)['UC']
                    else:
                        value = getattr(concept, param)

            else:
                # check if the key is part of the variant
                if key in variant:
                    # key is in variant, thus get value
                    value = variant[key][param]
                else:
                    # not all variants have a filter layer
                    value = 'None'

            # check if format of armour layer CRM must be changed
            if change_CRM_class and param == 'class':
                if not isinstance(value, str):
                    SUB = str.maketrans('3', '³')
                    unit = 'm3'.translate(SUB)
                    value = f'{np.round(value,3)} {unit}'
                else:
                    pass

            # add to the dict and round value
            if param in ['slope', 'slope_toe', 'slope_foundation']:
                if slopes == 'angles':
                    angle = np.arctan(value[0]/value[1]) * 180/np.pi
                    to_explorer[arg] = round(angle, 1)
                elif slopes == 'tuples':
                    to_explorer[arg] = [(round(value[0], 2), round(value[1], 2))]
                else:
                    raise InputError(
                        (f'method {slopes} is invalid, please use angles or '
                          'tuples'))
            elif isinstance(value, str):
                to_explorer[arg] = [value]
            elif value is None:
                to_explorer[arg] = [value]
            else:
                to_explorer[arg] = [round(value, 3)]
        else:
            # invalid arg was given for this structure
            # not raising error because might be for another structure
            # set value to 0 for correct DesignExplorer export
            to_explorer[arg] = [0]

    return to_explorer
