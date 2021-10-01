import warnings
import matplotlib.pyplot as plt

# define errors
class InputError(Exception):
    """ Exception raised for errors in the input

    Attributes
    ----------
    message
        explanation of the error
    """
    __module__ = Exception.__module__

    def __init__(self, message):
        self.message = message


class NotSupportedError(Exception):
    """ Exception raised when input is not supported

    Attributes
    ----------
    message
        explanation of the error
    """
    __module__ = Exception.__module__

    def __init__(self, message):
        self.message = message


class RockGradingError(Exception):
    """ Exception raised if there is an error in the Rock Grading

    Attributes
    ----------
    message
        explanation of the error
    """
    __module__ = Exception.__module__

    def __init__(self, message):
        self.message = message


class ArmourUnitsError(Exception):
    """ Exception raised if there is an error in the Armour Units

    Attributes
    ----------
    message
        explanation of the error
    """
    __module__ = Exception.__module__

    def __init__(self, message):
        self.message = message


# define warnings
class LimitStateWarning(Warning):
    """ Warning in the LimitState """

    def __init__(self, msg):
        self.message = msg

class EquipmentError(Warning):

    def __init__(self, msg, other, id,  areas, depth_area):

        self.message = msg

        fig, ax = plt.subplots(figsize=(10,5))
        coordinates = other._layers(id)

        for layer, lines in coordinates.items():
            plt.plot(lines["x"], lines["y"], color="k", linewidth= 3)

        for layer, area in areas.items():
            for key, value in depth_area[layer]['Area_yrange'].items():
                coords = depth_area[layer]['Area_yrange'][key]['coordinates']
                for c in coords:
                    ax.plot(c[0], c[1], color= 'k', linewidth= 1)
                    ax.fill(c[0], c[1], color= depth_area[layer]['Area_yrange'][key]['color'])


# monkeypatch warning format
def _custom_formatwarning(msg, category, *args, **kwargs):
    # ignore everything except the message
    return f'{category.__name__}: {msg} \n'

def user_warning(msg):
    """ show user warning """
    warnings.formatwarning = _custom_formatwarning
    warnings.warn(msg, category=UserWarning)

def limitstate_warning(msg):
    """ Show warning in the LimitState """
    warnings.formatwarning = _custom_formatwarning
    warnings.warn(msg, category=LimitStateWarning)

def no_warnings():
    """ ignore warnings """
    warnings.filterwarnings('ignore')
