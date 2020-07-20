import re

from .exceptions import NotSupportedError

def _convert_string(string, mode, param):
    """ Convert a string to a list or tuple

    As a tuple and list is converted to a string by xlsxwriter it must
    be converted back to the desired type.

    Parameters
    ----------
    string : str
        the string to convert
    mode : {slope, lambda, varying}
        the type of the output
    param : str
        the parameter to convert
    """
    # validate input
    if not isinstance(string, str):
        raise TypeError(f'{param} is not correctly formatted in Excel')

    # get numbers from the string
    numbers = re.findall(r"[-+]?\d*\.\d+|\d+", string)

    # convert numbers to floats
    converted = []
    for number in numbers:
        converted.append(float(number))

    if mode == 'slope':
        # convert to tuple
        return tuple(converted)

    elif mode == 'lambda':
        return list(converted)

    elif mode == 'varying':
        pass

    else:
        raise NotSupportedError(
            f'Mode {mode} is not supported, must be slope, lambda or varying')
