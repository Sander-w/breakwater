import numpy as np
import matplotlib.pyplot as plt

from .exceptions import InputError, RockGradingError

def _process_cost(structure, type, cost, Grading, validate=True):
    """ Process cost input to a dict

    Parameters
    ----------
    structure : {'RRM', 'CRM', 'RC', 'CC'}
        structure for which the cost must be verified
    type : {'Material', 'C02'}
        Type of cost to process
    Grading : :py:class:`RockGrading`
        rock grading
    cost : dict
        dictionary with the cost
    validate : bool, optional, default: True
        if True the input is validated, if False the cost dict is
        returned
    """
    # check if cost have been specified

    dictvar = None


    if cost is not None:
        # cost have been added
        # check if cost have been added to the grading
        if 'cost' in Grading[list(Grading.grading.keys())[0]]:
            # pricing has been added
            pass
        else:
            # pricing has not been added, raise error
            raise RockGradingError('There is no pricing in the RockGrading')

        # check structure to check structure specific cost
        if 'RRM' in structure:
            # check if core_price is in cost
            if 'core_price' not in cost.keys() or cost['core_price'] is None:
                raise KeyError(
                    'core_price must be specified when computing the cost of RRM')

        if 'CRM' in structure or 'CC' in structure:
            # check if unit_price is in cost
            if 'unit_price' not in cost.keys() or cost['unit_price'] is None:
                raise KeyError(
                    'unit_price must be specified when computing the cost of CRM/CC')

        if 'RC' in structure or 'CC' in structure:
            # check if fill_price is in cost
            if 'fill_price' not in cost.keys() or cost['fill_price'] is None:
                raise KeyError(
                    'fill_price must be specified when computing the cost of CRM/CC')

            # check if concrete_price is in cost
            if 'concrete_price' not in cost.keys() or cost['concrete_price'] is None:
                raise KeyError(
                    'concrete_price must be specified when computing the cost of CRM/CC')


        # check if cost must be returned
        if not validate:
            # check for optional parameters
            if 'transport_price' not in cost:
                # add as None
                cost['transport_price'] = None

            if 'dry_dock' not in cost:
                # add as None
                cost['transport_cost'] = None

            if 'length' not in cost:
                cost['length'] = None

            return cost

def cost_influence(type, lines):
    """ Plot influence of varying parameters

    Parameters
    ----------
    type = {'Material', 'C02'}
        Indicates whether the material or C02 costs are analysed
    lines : dict
        dictionary with the parameters as keys and a nested dict with
        the values and cost
    """
    # check if more than one parameter has been given

    if len(lines.keys()) > 1:
        # values must be normalised
        normalise = True

        # set title
        title = 'Influence of the varying parameters on the cost (normalised)'
        xlabel = 'normalised value of the parameter, from min to max'

        # set xlim
        xmin, xmax = 0, 1

    else:
        # normalising is not needed
        normalise = False

        # set title
        title = f'Influence of the varying parameters on the {type} cost'

    # create the figure
    plt.figure()

    # iterate over the lines
    for parameter, data in lines.items():
        # set empty list for x data
        x = []

        # check if data must be normalised
        if normalise:
            # get min and max value
            min = np.min(data['values'])
            max = np.max(data['values'])

            for value in data['values']:
                # normalise data and append to list
                x.append((value - min)/(max - min))
            # add min and max to label
            label = f'{parameter} (min={min}, max={max})'

        else:
            # x equals the values
            x = data['values']


            # label is parameter
            label = parameter
            xlabel = parameter

            # set xmax and xmax
            xmin, xmax = np.min(x), np.max(x)


        # plot data
        if type == 'Material':
            plt.plot(x, data['material_cost'], label=label)
        if type == 'CO2':
            plt.plot(x, data['CO2_cost'], label=label)

    # style figure

    plt.xlim(xmin, xmax)
    plt.xlabel(xlabel)
    plt.ylabel(f'{type} cost per m')
    plt.title(title)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()
