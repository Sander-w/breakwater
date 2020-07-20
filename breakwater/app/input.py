import numpy as np
import tkinter as tk
from tkinter.messagebox import showerror
from tkinter.ttk import Separator

from .settings import settings
from .interactive import InteractiveDesign
from .tooltip import ToolTip
from .footer import footer
from ..utils._excel_utils import _convert_string
from ..utils._kwarg_validator import _process_kwargs
from ..utils.exceptions import InputError


def _convert_entry(parameter, entry, vkwargs):
    """ Convert the entry field to the correct format

    Input of an entry field is returned as a string by Tkinter, this
    method converts the input of the entry field to the correct format

    Parameters
    ----------
    parameter : str
        name of the parameter
    entry : str
        input from the entry field
    vkwargs : dict
        dictionary of the valid kwargs for the specified structure

    Returns
    -------
    converted_entry : float, int, tuple, list, str
        input of the entry converted to the correct format
    """
    # get correct format from vkwargs
    values = vkwargs[parameter]

    if not entry:
        converted_entry = None

    elif 'float' in values['Correct']:
        converted_entry = float(entry)

    elif 'int' in values['Correct']:
        converted_entry = int(entry)

    elif 'tuple' in values['Correct'] or 'list' in values['Correct']:
        if 'slope' in parameter:
            mode = 'slope'
        elif 'lambda' in parameter:
            mode = 'lambda'
        else:
            # error
            showerror(
                'Error',
                f'input for {parameter} cannot be given as a tuple or list')

        converted_entry = _convert_string(entry, mode, parameter)

    elif 'str' in values['Correct']:
        converted_entry = entry

    else:
        # pop up with error message
        showerror('Error', f'input for {parameter} is not recognized')

    return converted_entry


class ParameterInput(tk.Frame):
    """ Input page of the app """

    def __init__(self, parent, controller):
        """ See help(ParameterInput) for more info """
        # create the frame
        tk.Frame.__init__(self, parent)

        # set controller as attribute
        self.controller = controller

        # get width and height
        width = settings('lay-out')['width']
        height = settings('lay-out')['height']

        # set height of frames
        top_height = settings('lay-out')['top_height']
        bottom_height = settings('lay-out')['bottom_height']
        center_height = height - top_height - bottom_height

        # set left space
        left_space = width/4

        # create lay out
        top_frame = tk.Frame(self, width=width, height=top_height)
        center_frame = tk.Frame(self, width=width, height=center_height)
        bottom_frame = tk.Frame(self, width=width, height=bottom_height)

        top_frame.pack(side=tk.TOP, fill=tk.X, pady=5)
        center_frame.pack(fill=tk.Y, expand=True)
        bottom_frame.pack(side=tk.BOTTOM, fill=tk.X)

        # set title
        label = tk.Label(
            top_frame, text='Parameter Input',  **settings('title'))
        label.pack()

        # set dict to store the entry field
        entries = {}
        checkboxes = {}

        # first row of the input table
        start_row = 0

        headers = ['Parameter', 'Unit', 'Input', 'Varying']
        for j, header in enumerate(headers):
            label = tk.Label(center_frame, text=header, **settings('label'))
            label.grid(row=start_row, column=j)

        start_row += 1

        # add an entry field for each parameter
        for i, (parameter, value) in enumerate(self.controller.vkwargs.items()):
            if value['Required']:
                required = 'required'
            else:
                required = 'optional'

            format = value['Correct']

            # create label for the parameter
            label = tk.Label(center_frame, text=parameter, **settings('label'))
            label.grid(row=start_row+i, column=0)
            label_tooltip = ToolTip(
                label,
                text=f'Parameter is {required}. Must be formatted as {format}')

            # create label for the unit
            unit = self.controller.vkwargs[parameter]['unit']
            label = tk.Label(center_frame, text=unit, **settings('label'))
            label.grid(row=start_row+i, column=1)

            # create entry fields
            entry = tk.Entry(center_frame)
            entry.grid(row=start_row+i, column=2)
            entry_tooltip = ToolTip(
                entry,
                text=f'Parameter is {required}. Must be formatted as {format}')

            # set default if there is one
            default = value['Default']
            if default is not None:
                entry.insert(i, str(default))

            if not value['Constant']:
                # add checkbox to initiate varying parameter
                var = tk.BooleanVar()
                checkbox = tk.Checkbutton(center_frame, variable=var)
                checkbox.grid(row=start_row+i, column=3)

                # add checkbox to dict
                checkboxes[parameter] = {'var': var, 'check': checkbox}

            # add entry to dict
            entries[parameter] = entry

        # add two empty columns to center input table
        center_frame.grid_columnconfigure(0, weight=1)
        center_frame.grid_columnconfigure(4, weight=1)

        # add lines between the columns
        Separator(center_frame, orient=tk.VERTICAL).grid(
            column=0, row=start_row-1, rowspan=i+2, sticky='nse')
        Separator(center_frame, orient=tk.VERTICAL).grid(
            column=1, row=start_row-1, rowspan=i+2, sticky='nse')

        # add line below the headers
        Separator(center_frame, orient=tk.HORIZONTAL).grid(
            column=0, row=start_row-1, columnspan=4, sticky='ews')

        # add footer
        footer(bottom_frame)

        # add navigation buttons
        nav_settings = settings('navigation')

        start_interactive = tk.Button(
            bottom_frame, text='Start Design', height=nav_settings['height'],
            width=nav_settings['width'],
            command=lambda: self.next_page(entries, checkboxes))
        start_interactive.pack(
            side=tk.RIGHT, padx=nav_settings['padx'], pady=nav_settings['pady'])

        back = tk.Button(
            bottom_frame, text='Back', height=nav_settings['height'],
            width=nav_settings['width'],
            command=lambda: self.controller.home())
        back.pack(
            side=tk.LEFT, padx=nav_settings['padx'], pady=nav_settings['pady'])

    def process_entry(self, entries, varying_parameters):
        """ Process the given entry """
        # set bool to validate the processing
        abort = False

        # output
        output = {}

        # iterate over the entry fields
        for parameter, entry in entries.items():
            if parameter in varying_parameters:
                # given as varying parameter thus pass
                continue
            else:
                # process parameter
                converted_entry = _convert_entry(
                    parameter, entry.get(), self.controller.vkwargs)

                # check if input is valid
                validator = self.controller.vkwargs[parameter]['Validator']
                valid = validator(converted_entry)

                # check if parameter is required
                required = self.controller.vkwargs[parameter]['Required']

                if valid or not required:
                    # add to dict
                    output[parameter] = converted_entry

                else:
                    # error and abort
                    abort = True
                    correct = self.controller.vkwargs[parameter]['Correct']
                    showerror(
                        'InputError',
                        f'{parameter} is incorrectly formatted, '
                        f'must be {correct}')

        # add output parameters to the controller
        self.controller.parameters.update(output)

        return abort

    def process_checkboxes(self, checkboxes):
        """ Process the input from the checkboxes """
        varying_parameters = []

        # iterate over the checkboxes to see if one has been checked
        for parameter, checkbox in checkboxes.items():
            # check if the checkbox was checked
            if checkbox['var'].get():
                # checkbox was checked
                varying_parameters.append(parameter)
            else:
                # checkbox was not checked
                pass

        # return the varying parameters
        return varying_parameters

    def next_page(self, entries, checkboxes):
        """ Process the input and go to the next page """
        # process varying input
        varying_parameters = self.process_checkboxes(checkboxes)

        # process the entry
        abort = self.process_entry(entries, varying_parameters)

        if not abort:
            # check if varying input must be given
            if any(varying_parameters):
                # update attribute for the previous page of InteractiveDesign
                self.controller.previous = 'VaryingInput'

                # go to varying input page
                self.controller.frames[VaryingInput].create_frame(
                    varying_parameters, self.controller.vkwargs)

                self.controller.show_frame(VaryingInput)
            else:
                # update attribute for the previous page of InteractiveDesign
                self.controller.previous = 'ParameterInput'

                # design breakwater for next frame
                self.controller.frames[InteractiveDesign].create_frame()

                # show next frame
                self.controller.show_frame(InteractiveDesign, default=False)


class VaryingInput(tk.Frame):
    """ Optional input frame for giving varying input """

    def __init__(self, parent, controller):
        """ See help(VaryingInput) for more info """
        # create the frame
        tk.Frame.__init__(self, parent)

        # set controller as attribute
        self.controller = controller
        self.first_generation = True

        # get width and height
        width = settings('lay-out')['width']
        height = settings('lay-out')['height']

        # set height of frames
        top_height = settings('lay-out')['top_height']
        bottom_height = settings('lay-out')['bottom_height']
        center_height = height - top_height - bottom_height

        # create lay out
        top_frame = tk.Frame(self, width=width, height=top_height)
        center_frame = tk.Frame(self, width=width, height=center_height)
        bottom_frame = tk.Frame(self, width=width, height=bottom_height)

        top_frame.pack(side=tk.TOP, fill=tk.X, pady=5)
        center_frame.pack(fill=tk.Y, expand=True)
        bottom_frame.pack(side=tk.BOTTOM, fill=tk.X)

        # set frames as attributes
        self.top = top_frame
        self.center = center_frame
        self.bottom = bottom_frame

        # set title of the page
        label = tk.Label(
            self.top, text='Varying Parameter Input', **settings('title'))
        label.pack()

        # add the footer
        footer(self.bottom)

        # create back button
        nav_settings = settings('navigation')

        back = tk.Button(
            self.bottom, text='Back', height=nav_settings['height'],
            width=nav_settings['width'],
            command=lambda: self.controller.show_frame(ParameterInput))
        back.pack(
            side=tk.LEFT, padx=nav_settings['padx'], pady=nav_settings['pady'])

    def create_frame(self, parameters, vkwargs):
        """ Create the frame for giving varying input

        Method to create the frame for when the user wants to specify
        varying parameters.

        Parameters
        ----------
        parameters : list
            list of parameters to variate
        vkwargs : dict
            dict of the valid kwargs for the specified structures
        """
        # check if this is the first time the frame is generated
        if self.first_generation:
            # first time the frame is generated
            # generate entry fields
            self.generate_input(parameters, vkwargs)

            # add start design button
            self.start_interactive_btn()

            # no longer first time the frame is generated
            self.first_generation = False

        else:
            # not the first time
            # clear the widgets of the frame
            for widget in self.center.winfo_children():
                widget.destroy()

            # generate entry fields
            self.generate_input(parameters, vkwargs)

    def start_interactive_btn(self):
        """ Create the start design button """
        nav_settings = settings('navigation')

        # add button to start designing
        start_interactive = tk.Button(
            self.bottom, text='Start Design', height=nav_settings['height'],
            width=nav_settings['width'],
            command=lambda: self.next_page())
        start_interactive.pack(
            side=tk.RIGHT, padx=nav_settings['padx'], pady=nav_settings['pady'])

    def generate_input(self, parameters, vkwargs):
        """ Generate input fields for the varying input

        Method to create the entry fields for the parameters the user
        wants to variate.

        Parameters
        ----------
        parameters : list
            list of parameters to variate
        vkwargs : dict
            dict of the valid kwargs for the specified structures
        """
        # set dict to store the entry field
        entries = {}

        # first row of the input table
        start_row = 0

        headers = ['Parameter', 'Unit', 'Minimum value', 'Maximum value']
        for j, header in enumerate(headers):
            label = tk.Label(self.center, text=header, **settings('label'))
            label.grid(row=start_row, column=j)

        start_row += 1

        # generate the input fields
        for i, parameter in enumerate(parameters):
            # add label for the parameter
            label = tk.Label(self.center, text=parameter, **settings('label'))
            label.grid(row=start_row+i)

            # add label for the unit
            unit = self.controller.vkwargs[parameter]['unit']
            label = tk.Label(self.center, text=unit, **settings('label'))
            label.grid(row=start_row+i, column=1)

            # create entry fields
            entry_min = tk.Entry(self.center)
            entry_max = tk.Entry(self.center)
            entry_min.grid(row=start_row+i, column=2)
            entry_max.grid(row=start_row+i, column=3)

            # add entry to dict
            entries[parameter] = {'min': entry_min, 'max': entry_max}

        # add lines between the columns
        Separator(self.center, orient=tk.VERTICAL).grid(
            column=0, row=start_row-1, rowspan=i+2, sticky='nse')
        Separator(self.center, orient=tk.VERTICAL).grid(
            column=1, row=start_row-1, rowspan=i+2, sticky='nse')

        # add line below the headers
        Separator(self.center, orient=tk.HORIZONTAL).grid(
            column=0, row=start_row-1, columnspan=4, sticky='ews')

        # add two empty columns to center input table
        self.center.grid_columnconfigure(0, weight=1)
        self.center.grid_columnconfigure(4, weight=1)

        # set entries as attribute
        self.entries = entries

    def get_input(self):
        """ Get the input from the entry fields

        Returns
        -------
        dict
            dictionary of the varying parameters with the specified
            minimum and maximum value
        """
        # set bool to validate the processing
        abort = False

        # create dict to store converted entries
        output = {}

        # iterate over the entry fields of the parameters
        for parameter, entry in self.entries.items():
            # create temp dict to store min and max
            temp = {}

            # create variable to store the number of decimals
            num_decimals = 0

            # get the validator and required bool
            validator = self.controller.vkwargs[parameter]['Validator']

            # iterate over the entry field min and max
            for name, value in entry.items():
                # determine the number of decimals in the input
                split = value.get().split('.')

                # check if there are decimals
                if len(split) > 1:
                    # determine number of decimals
                    computed_num_decimals = len(split[-1])
                else:
                    # no decimals
                    computed_num_decimals = 0

                # check if larger than previous num_decimals
                if computed_num_decimals >= num_decimals:
                    num_decimals = computed_num_decimals

                # convert entry field to correct type and store in temp
                converted_entry = _convert_entry(
                    parameter, value.get(), self.controller.vkwargs)

                # check if input is valid
                validator = self.controller.vkwargs[parameter]['Validator']
                valid = validator(converted_entry)

                # check if parameter is valid
                if valid:
                    # check if parameter is a slope parameter
                    if 'slope' in parameter:
                        # set custom number of decimals for degrees
                        num_decimals = 0

                        # convert slope to angles
                        try:
                            converted_entry = np.arctan(
                                converted_entry[0]/converted_entry[1])
                            converted_entry = np.round(
                                converted_entry * 180/np.pi)
                        except IndexError:
                            # parameter not given as a tuple
                            msg = (f'{parameter} must be specified as a '
                                    'tuple (V, H)')
                            showerror('Error', msg)
                            abort = True

                    # add to dict
                    temp[name] = converted_entry

                else:
                    # error and abort
                    abort = True
                    correct = self.controller.vkwargs[parameter]['Correct']
                    showerror(
                        'InputError',
                        f'{name} of {parameter} is incorrectly formatted, '
                        f'must be {correct}')

            # add number of decimals to the dict and temp to output
            temp['resolution'] = num_decimals
            output[parameter] = temp

        return output, abort

    def next_page(self):
        """ Process the input and go to the next page """
        # process varying input
        varying, abort = self.get_input()

        # check if getting input was successful
        if not abort:
            # update parameters with values from the sliders
            for parameter, values in varying.items():
                # get the minimum value
                value = values['min']

                # check if parameter is a slope parameter
                if 'slope' in parameter:
                    # convert angle back to tuple
                    value = np.tan(value*np.pi/180)
                    value = (1, 1/value)

                # add minimum slider value to the parameters
                self.controller.parameters[parameter] = value


            # design breakwater for next frame
            self.controller.frames[InteractiveDesign].create_frame(
                                                            varying=varying)

            # show next frame
            self.controller.show_frame(InteractiveDesign, default=False)
