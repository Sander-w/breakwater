import matplotlib
import numpy as np
import tkinter as tk
from tkinter.ttk import Treeview
from tkinter.messagebox import showerror, showinfo
from warnings import catch_warnings
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

from .settings import settings, size
from .footer import footer
from ..rubble import RockRubbleMound, ConcreteRubbleMound
from ..caisson import Caisson
from ..utils.exceptions import RockGradingError, ArmourUnitsError
from ..utils._kwarg_validator import _RM_vkwargs, _C_vkwargs
from ..utils.cost import cost_influence


class InteractiveDesign(tk.Frame):
    """ Interactive design page of the app """

    def __init__(self, parent, controller):
        """ See help(ParameterInput) for more info """
        tk.Frame.__init__(self, parent)

        # set controller as attribute
        self.controller = controller

        # set attribute for the sliders
        self.sliders = {}

        # create attribute to store varying parameters of the bw_types
        self.varying = {'RM': {}, 'C': {}}

        # set attribute to track button generation
        self.buttons = {'update': False, 'influence': False}

        # set attribute for the limits of the xaxis
        self.xlim = None

        # set variables
        width, height = size(
            self.controller.to_design, 'interactive',
            self.controller.num_concepts)

        # determine dimensions of the frames
        top_height = settings('lay-out')['top_height']
        bottom_height = settings('lay-out')['bottom_height']

        # set top and bottom frame
        top_frame = tk.Frame(self, width=width, height=top_height)
        self.bottom_frame = tk.Frame(self, width=width, height=bottom_height)

        # pack frames
        top_frame.pack(side=tk.TOP, fill=tk.X, pady=5)
        self.bottom_frame.pack(side=tk.BOTTOM, fill=tk.X)

        # generate design frames
        self._generate_design_frames(width, height)

        # add the title of the page
        label = tk.Label(
            top_frame, text='Interactive Design', **settings('title'))
        label.pack()

        # add the footer
        footer(self.bottom_frame)

        # add back button
        nav_settings = settings('navigation')

        back = tk.Button(
            self.bottom_frame, text='Back', height=nav_settings['height'],
            width=nav_settings['width'], command=lambda: self.back())
        back.pack(
            side=tk.LEFT, padx=nav_settings['padx'], pady=nav_settings['pady'])

    def _generate_design_frames(self, width, height):
        """ Generate the frames for the designs

        Parameters
        ----------
        width, height : int
            width and height of the super frame in pixels
        """
        # set attribute for the design frames and sliders
        self.design_frames = {}
        self.slider_frames = {}

        # frame variables
        design_width, design_height = size(
            self.controller.to_design, 'design', self.controller.num_concepts)
        info_width = settings('interactive')['info_width']
        slider_width = width - design_width - info_width

        # generate RM design frames
        for bw_type, values in self.controller.to_design.items():
            # check the number of variants for this type
            num = values['num']
            if num > 0:
                # create parent frame
                parent = tk.Frame(self, width=width, height=design_height*num)
                parent.pack_propagate(False)
                parent.pack()

                # add frame for the sliders
                sliders = tk.Frame(
                    parent, width=slider_width, height=design_height*num)
                sliders.pack_propagate(False)
                sliders.pack(side=tk.LEFT)

                # create parent for the design frames
                design_parent = tk.Frame(
                    parent, width=design_width+info_width,
                    height=design_height*num)
                design_parent.pack_propagate(False)
                design_parent.pack(side=tk.RIGHT)

                # create dict to store the design frames of this type
                design_frames = {}

                # remove num from values
                del values['num']

                # create frames for the bw_type
                for structure, active in values.items():
                    # check if structure is active
                    if active:
                        # structure is active thus create the required frames
                        # create parent for the types
                        type_parent = tk.Frame(
                            design_parent, width=design_width+info_width,
                            height=design_height)
                        type_parent.pack_propagate(False)
                        type_parent.pack()

                        # create plot frame
                        plot_frame = tk.Frame(
                            type_parent, width=design_width, height=design_height)
                        plot_frame.pack_propagate(False)
                        plot_frame.pack(side=tk.LEFT, padx=0, pady=0)

                        # create info frame
                        info_frame = tk.Frame(
                            type_parent, width=info_width, height=design_height)
                        info_frame.pack(side=tk.RIGHT, padx=5)

                        # get the name of the structure and add frames to dict
                        design_frames[structure] = {
                            'plot': plot_frame, 'info': info_frame}

                # add frames to attribute
                self.design_frames[bw_type] = design_frames
                self.slider_frames[bw_type] = sliders

    def _display_warning(self, w, structure, show_warnings):
        """ display the warning messages encountered during the design

        Parameters
        ----------
        w : list
            list of warnings
        structure : str
            name of the structure to design
        show_warnings : bool
            True if warnings must be displayed, False if not
        """
        # check if there are warnings
        if len(w) > 0 and show_warnings:
            # there are warnings
            msg = 'Warnings: \n'
            for warning in w:
                msg += f'- {warning.message} \n'

            showinfo(f'{structure} Warning', msg)

        else:
            # no warnings
            pass

    def _format_input(self, structure):
        """ Format the input for the design class

        Function to format the parameters to the correct format for the
        design classes. This means that for the RC the BermMaterial
        argument is set to the Grading.

        Parameters
        ----------
        structure : str
            name of the structure to design

        Returns
        -------
        dict
            with the parameters correctly formatted
        """
        # check structure
        if structure == 'RRM' or structure == 'CRM':
            # no formatting required
            return self.controller.parameters

        elif structure == 'RC' or structure == 'CC':
            # get the parameters
            output = self.controller.parameters.copy()

            # caisson structure, remove layers rock and units
            if 'layers_rock' in output.keys():
                del output['layers_rock']

            if 'layers_rock' in output.keys():
                del output['layers_units']

            if structure == 'RC':
                if 'BermMaterial' in output.keys():
                    del output['BermMaterial']

            return output

    def back(self):
        """ Back to the previous page """
        # remove the xlims of the previous plot
        self.xlim = None

        # clear the buttons in the bottom frame
        to_delete = ['Update Design', 'Cost influence']
        for widget in self.bottom_frame.winfo_children():
            if widget.cget('text') in to_delete:
                widget.destroy()

        self.buttons = {'update': False, 'influence': False}

        # clear the figures
        for bw_type, structures in self.design_frames.items():
            for structure, frames in structures.items():
                # destroy all widgets
                for widget in frames['plot'].winfo_children():
                    widget.destroy()

                for widget in frames['info'].winfo_children():
                    widget.destroy()

                for widget in self.slider_frames[bw_type].winfo_children():
                    widget.destroy()

                # delete the designs
                if hasattr(self, f'{structure}_old'):
                    delattr(self, f'{structure}_old')

                if hasattr(self, f'{structure}_new'):
                    delattr(self, f'{structure}_new')

        # destroy the sliders
        for bw_type, frame in self.slider_frames.items():
            for widget in frame.winfo_children():
                widget.destroy()

        # go to the previous page
        self.controller.back(self.controller.previous)

    def create_frame(self, varying=None):
        """ Create frame

        Method to create the page when the ParameterInput, and
        VaryingInput (if this is specified), is given.

        Parameters
        ----------
        varying : dict, optional, default: None
            dictionary of the varying parameters with the specified
            minimum and maximum value
        """
        # iterate over the design frames to plot the design and show info
        for bw_type, structures in self.design_frames.items():
            # check if varying parameters have been given
            if varying is not None:
                # varying parameters, thus generate sliders
                self.generate_sliders(
                    varying, self.slider_frames[bw_type], bw_type)

                # check if cost have been specified
                if self.controller.cost is not None:
                    # create cost button
                    nav_settings = settings('navigation')

                    if not self.buttons['influence']:
                        cost_button = tk.Button(
                            self.bottom_frame, text='Cost influence',
                            height=nav_settings['height'],
                            width=nav_settings['width'],
                            command=lambda: self.cost_influence())
                        cost_button.pack(
                            side=tk.RIGHT, padx=nav_settings['padx'],
                            pady=nav_settings['pady'])

                        self.buttons['influence'] = True

            else:
                # no varying input
                no_sliders(self.slider_frames[bw_type])

            # iterate over the structures of the type
            for structure, frames in structures.items():
                # show the breakwater
                self.show_breakwater(structure, frames)

    def generate_sliders(self, varying, frame, bw_type):
        """ Add sliders for the varying parameters

        Method to create the sliders for interactive design

        Parameters
        ----------
        varying : dict
            dictionary of the varying parameters with the specified
            minimum and maximum value
        frame : tk.Frame
            frame on which the sliders must be placed
        bw_type : str
            type of breakwater
        """
        # determine the width of a slider
        width = settings('interactive')['width']
        design_width = settings('interactive')['design_width']
        info_width = settings('interactive')['info_width']
        frame_width = width - design_width - info_width

        width = _slider_width(
            parameters=varying.keys(), frame_width=frame_width)

        # create bool to check if a slider has been placed
        placed_slider = False

        # get the vkwargs of the structure
        if bw_type == 'RM':
            vkwargs = _RM_vkwargs(type='both')
        elif bw_type == 'C':
            vkwargs = _C_vkwargs(type='both')
        else:
            raise NotImplementedError(f'{bw_type} has not been implemented')

        # create frame for the sliders
        parent = tk.Frame(
            frame, width=frame.winfo_reqwidth(), height=frame.winfo_reqheight())
        parent.place(anchor=tk.CENTER, relx=.5, rely=.5)

        # iterate over the varying parameters
        for parameter, values in varying.items():
            # check if parameter is valid for the given structure
            if parameter in vkwargs.keys():
                # make slider and add slider to attribute
                self.sliders[parameter] = make_slider(
                    parameter, values, parent, width)

                # change bool
                placed_slider = True

                # add to dict
                self.varying[bw_type][parameter] = varying[parameter]

        # check if a slider has been placed
        if not placed_slider:
            # no slider, add label
            no_sliders(frame)

        # create update button
        nav_settings = settings('navigation')

        if not self.buttons['update']:
            update_button = tk.Button(
                self.bottom_frame, text='Update Design',
                height=nav_settings['height'], width=nav_settings['width'],
                command=lambda: self.update())
            update_button.pack(
                side=tk.RIGHT, padx=nav_settings['padx'],
                pady=nav_settings['pady'])

            self.buttons['update'] = True

    def update(self):
        """ Update design with new varying parameters

        Method to update the design with the new parameters from the
        sliders
        """
        # iterate over the sliders
        for parameter, value in self.sliders.items():
            # get the value of the slider
            value = value.get()

            # check if parameter is a slope parameter
            if 'slope' in parameter:
                # convert angle back to tuple
                value = np.tan(value*np.pi/180)
                value = (1, 1/value)

            self.controller.parameters[parameter] = value

        # remove the xlims of the previous plot
        self.xlim = None

        # iterate over the design frames to plot the design and show info
        for bw_type, structures in self.design_frames.items():
            # iterate over the structures of the type
            for structure, frames in structures.items():
                # delete the old plot
                for widget in frames['plot'].winfo_children():
                    widget.destroy()
                for widget in frames['info'].winfo_children():
                    widget.destroy()

                # show the breakwater
                self.show_breakwater(structure, frames)

    def show_breakwater(self, structure, frames):
        """ Method to generate the plot and info

        Parameters
        ----------
        structure : str
            name of the structure to design
        frames : dict
            dictionary of tk.Frames
        """
        # design the breakwater
        self.design(structure, show_warnings=self.controller.display_warnings)

        # compute the cost and get cheapest variant
        cost = self.compute_cost(structure)

        if cost is not None:
            min_variant = min(cost, key=cost.get)
            price = cost[min_variant]
        else:
            min_variant = 'a'
            price = None

        # add the plot the frame
        self.add_plot(structure, frames['plot'], min_variant)

        # add the info of the design
        self.add_info(structure, frames['info'], min_variant, price)

    def design(self, structure, show_warnings, show_error=True):
        """ Design a breakwater

        Parameters
        ----------
        structure : str
            name of the structure to design
        show_warnings : bool
            True if warnings must be displayed, False if not
        show_error : bool, optional, default: True
            True if errors must be shown, False if not
        """
        # check if a design has already been made
        if hasattr(self, f'{structure}_old'):
            # not the first design
            # set previous design as old design
            previous_design = getattr(self, f'{structure}_new')
            setattr(self, f'{structure}_old', previous_design)

        else:
            # first design, set attributes for old design
            setattr(self, f'{structure}_old', None)

        # format input parameters
        input = self._format_input(structure)

        if structure == 'RRM':
            # design rubble mound breakwater with rock
            try:
                with catch_warnings(record=True) as w:
                    self.RRM_new = RockRubbleMound(**input)

                self._display_warning(w, 'RRM', show_warnings)

            except RockGradingError as e:
                if show_error:
                    showerror('RockGradingError', f'Computed {e}')
                delattr(self, 'RRM_old')

        elif structure == 'CRM':
            # design rubble mound breakwater with armour units
            try:
                with catch_warnings(record=True) as w:
                    self.CRM_new = ConcreteRubbleMound(**input)

                self._display_warning(w, 'CRM', show_warnings)

            except ArmourUnitsError as e:
                if show_error:
                    showerror('ArmourUnitsError', f'Computed {e}')
                delattr(self, 'CRM_old')

        elif structure == 'RC':
            # design caisson breakwater with rock as armour for the foundation
            try:
                with catch_warnings(record=True) as w:
                    self.RC_new = Caisson(
                        layers=self.controller.parameters['layers_rock'],
                        BermMaterial=self.controller.parameters['Grading'],
                        **input)

                self._display_warning(w, 'RC', show_warnings)

            except RockGradingError as e:
                if show_error:
                    showerror('RockGradingError', f'Computed {e}')
                delattr(self, 'RC_old')

        elif structure == 'CC':
            # design caisson breakwater with armour units as armour
            try:
                with catch_warnings(record=True) as w:
                    self.CC_new = Caisson(
                        layers=self.controller.parameters['layers_units'],
                        **input)

                self._display_warning(w, 'CC', show_warnings)

            except ArmourUnitsError as e:
                if show_error:
                    showerror('ArmourUnitsError', f'Computed {e}')
                delattr(self, 'CC_old')

        else:
            # invalid structure
            raise NotImplementedError(f'{structure} has not been implemented')

    def compute_cost(self, structure, output='variant'):
        """ Determine the most expensive variant

        Parameters
        ----------
        structure : str
            name of the structure to design
        output : {variant, layer, average}, optional, default: variant
            format of the output dict, variant returns the total cost
            of each variant, layer the cost of each layer for each
            variant and average returns the average cost.
        Returns
        -------
        dict
            the cost of each variant
        """
        # get the attribute
        design = getattr(self, f'{structure}_new')

        # check if cost have been specified
        if self.controller.cost is not None:
            # cost have been added
            if structure == 'RRM':
                # rubble mound structure with rock as armour
                cost = design.cost(
                    *design.variantIDs,
                    core_price=self.controller.cost['core_price'],
                    transport_cost=self.controller.cost['transport_price'],
                    output=output)

            elif structure == 'CRM':
                # rubble mound structure with ArmourUnits as armour
                cost = design.cost(
                    *design.variantIDs,
                    core_price=self.controller.cost['core_price'],
                    unit_price=self.controller.cost['unit_price'],
                    transport_cost=self.controller.cost['transport_price'],
                    output=output)

            elif structure == 'RC' or structure == 'CC':
                # get rent of dry dock and length of the breakwater
                rent = self.controller.cost['dry_dock']
                length = self.controller.cost['length']

                # check if dry dock investment must be added
                if rent is not None and length is not None:
                    # add rent of dry dock
                    design.dry_dock(rent, length)

                # check if rock armour layer
                if structure == 'RC':
                    # caisson structure with rock as armour of the foundation
                    cost = design.cost(
                        *design.variantIDs,
                        concrete_price=self.controller.cost['concrete_price'],
                        fill_price=self.controller.cost['fill_price'],
                        output=output)

                else:
                    # caisson structure with units as armour of the foundation
                    cost = design.cost(
                        *design.variantIDs,
                        concrete_price=self.controller.cost['concrete_price'],
                        fill_price=self.controller.cost['fill_price'],
                        unit_price=self.controller.cost['unit_price'],
                        output=output)

            else:
                # invalid structure
                raise NotImplementedError(
                    f'{structure} has not been implemented')

            # return the computed cost
            return cost

        else:
            return None

    def cost_influence(self, num=50):
        """ Plot the influence of the varying parameters on the cost """
        # make dict to store designs
        store_designs = {}

        # make variable to store current parameters
        parameters = self.controller.parameters.copy()

        # make dict to store the cost data and values
        lines = {}

        # iterate over the structures to design
        for bw_type, structures in self.controller.to_design.items():
            # get varying parameters of this breakwater type
            varying = self.varying[bw_type]
            # iterate over the structures
            for structure, active in structures.items():
                # check if structure is active
                if active:
                    # get current old and new design, and store in dict
                    store_designs[structure] = {
                        'old': getattr(self, f'{structure}_old'),
                        'new': getattr(self, f'{structure}_new')}

                    # compute influence of varying parameter on the cost
                    # iterate over the varying parameters
                    for parameter, values in varying.items():
                        # reset parameters back to original parameters
                        self.controller.parameters = parameters.copy()

                        # make linspace
                        varying = np.linspace(
                            values['min'], values['max'], num)

                        # add to line data
                        lines[parameter] = {'values': varying, 'cost': []}

                        # iterate over the values
                        for value in varying:
                            # check if parameter is a slope parameter
                            if 'slope' in parameter:
                                # convert angle back to tuple
                                value = np.tan(value*np.pi/180)
                                value = (1, 1/value)

                            self.controller.parameters[parameter] = value

                            # update design and get cost
                            self.design(
                                structure, show_warnings=False,
                                show_error=False)

                            with catch_warnings(record=True) as w:
                                cost = self.compute_cost(
                                    structure, output='average')

                            # add to line
                            lines[parameter]['cost'].append(
                                list(cost.values())[0])

        # create figure
        cost_influence(lines)

        # set design and parameters back to original ones
        for structure, concepts in store_designs.items():
            for key, concept in concepts.items():
                setattr(self, f'{structure}_{key}', concepts[key])

        self.controller.parameters = parameters

    def add_info(self, structure, frame, variantID, price):
        """ Add info of the plotted structure

        Parameters
        ----------
        structure : str
            name of the structure to design
        frame : tk.Frame
            frame to add the info to
        variantID : str
            identifier of the variant of which the info must be displayed
        price : float
            price of the specified variant
        """
        # get the dimensions of the frame
        param_width = 150
        value_width = settings('interactive')['info_width'] - param_width - 10

        # get the font
        font = settings('font')['info']

        # get the concept breakwater and variant
        concept = getattr(self, f'{structure}_new')
        variant = concept.get_variant(variantID)

        # create the treeview
        tree = Treeview(frame)

        # add columns
        tree['columns'] = ('#1')
        tree.column('#0', width=param_width, stretch=tk.NO)
        tree.heading('#0',text='Parameter', anchor=tk.W)
        tree.column('#1', width=value_width, stretch=tk.NO)
        tree.heading('#1', text='Value', anchor=tk.W)

        # check if a price has been specified
        if price is not None:
            # add cost
            tree.insert('', 1, iid=1, text='Cost', values=(np.round(price,2)))

            # set row
            row = 2

        else:
            # no cost, set row to first one
            row = 1

        # check if RubbleMound structure
        if structure == 'RRM' or structure == 'CRM':
            # parameters to show
            includes = ['computed Dn50', 'class']

            # add Rc
            tree.insert(
                '', row, iid=row, text='Rc',
                values=(f'{np.round(concept.Rc, 3)}m'))

            # increment row
            row += 1

            # add the info
            for layer, dimensions in variant.items():
                # add to first column
                header = tree.insert(
                    '', row, iid=row, text=layer.capitalize(), values=(''))

                # increment row
                row += 1

                # iterate over the dimensions of the layer
                for parameter, value in dimensions.items():
                    # check if value is an int or float
                    if isinstance(value, int) or isinstance(value, float):
                        if layer == 'armour' and parameter == 'class':
                            if structure == 'CRM':
                                SUB = str.maketrans('3', '³')
                                unit = 'm3'.translate(SUB)
                            else:
                                unit = 'm'
                        else:
                            unit = 'm'
                        # round the dimension to 3 decimals
                        value = f'{np.round(value, 3)}{unit}'

                    if parameter in includes:
                        # parameter and value
                        tree.insert(
                            header, 'end', iid=row,
                            text=parameter.capitalize(), values=(value))

                        # increment the row
                        row += 1

            # pack the tree to the frame
            tree.pack(side=tk.TOP, fill=tk.Y)

        # check if Caisson breakwater
        elif structure == 'RC' or structure == 'CC':
            # parameters to show
            includes = ['computed Dn50', 'class', 'd', 'B', 'h_acc']
            no_capital = ['h_acc', 'd']

            # add Rc
            Rc = variant['caisson']['Rc']
            tree.insert(
                '', row, iid=row, text='Rc', values=(f'{np.round(Rc, 3)}m'))

            # increment row
            row += 1

            # add the info
            for layer, dimensions in variant.items():
                # add to first column
                header = tree.insert(
                    '', row, iid=row, text=layer.capitalize(), values=(''))

                # increment row
                row += 1

                # iterate over the dimensions of the layer
                for parameter, value in dimensions.items():
                    # check if value is an int or float
                    if isinstance(value, int) or isinstance(value, float):
                        if layer == 'armour' and parameter == 'class':
                            if structure == 'CC':
                                SUB = str.maketrans('3', '³')
                                unit = 'm3'.translate(SUB)
                            else:
                                unit = 'm'
                        else:
                            unit = 'm'

                        # round the dimension to 3 decimals
                        value = f'{np.round(value, 3)}{unit}'

                    if parameter in includes:
                        if parameter not in no_capital:
                            parameter = parameter.capitalize()
                        # parameter and value
                        tree.insert(
                            header, 'end', iid=row, text=parameter,
                            values=(value))

                        # increment the row
                        row += 1

            # add geotechnical parameters if computed
            if self.controller.parameters['Soil'] is not None:
                # define the stresses
                stresses = ['p', 'p_cap']

                # add header of the geotechnical parameters
                geo = tree.insert(
                    '', row, iid=row, text='Geotechnical', values=(''))

                # increment row
                row += 1

                # add values
                for parameter, value in concept.geotechnical.items():
                    # check if parameter is a stess
                    if parameter in stresses:
                        # stress, round value and add unit
                        value = f'{np.round(value, 2)}kPa'
                    else:
                        # only round value
                        value = f'{np.round(value, 2)}'

                    # add to the tree and increment row
                    tree.insert(
                        geo, 'end', iid=row, text=parameter, values=(value))
                    row += 1

            # pack the tree to the frame
            tree.pack(side=tk.TOP, fill=tk.Y)

        else:
            # other type
            raise NotImplementedError(f'{structure} has not been implemented')

    def add_plot(self, structure, frame, variantID):
        """ Plot the cross section of the specified breakwater

        Parameters
        ----------
        structure : str
            name of the structure to design
        frame : tk.Frame
            frame to add the plot to
        variantID : str
            identifier of the variant to plot
        """
        # plot variables
        dpi = settings('interactive')['dpi']
        design_width, design_height = size(
            self.controller.to_design, 'design', self.controller.num_concepts)
        tick_fontsize = 6

        # create the figure
        fig = Figure(
            figsize=(design_width/dpi, design_height/dpi), dpi=dpi,
            facecolor='#f0f0f0')

        # add subplot to the figure
        ax = fig.add_subplot(111)

        # get the previous design
        old_design = getattr(self, f'{structure}_old')

        # check if old design must be plotted
        if old_design is not None:
            # get coordinates of the old design
            coordinates = old_design._layers(variantID)

            # iterate over the layers
            for layer, lines in coordinates.items():
                # plot each layer
                ax.plot(
                    lines['x'], lines['y'], color='grey', ls='--', lw=1,
                    zorder=1)

        # get the new design
        new_design = getattr(self, f'{structure}_new')

        # get coordinates of the new design
        coordinates = new_design._layers(variantID)

        # iterate over the layers
        for layer, lines in coordinates.items():
            # plot each layer
            ax.plot(lines['x'], lines['y'], color='k', zorder=20, lw=1)

        # set xlim
        xmin = ax.get_xlim()[0]
        xmax = ax.get_xlim()[1]

        if self.xlim is None:
            self.xlim = {'xmin': xmin, 'xmax': xmax}
        else:
            if np.abs(self.xlim['xmin']) > np.abs(xmin):
                xmin = self.xlim['xmin']

            if np.abs(self.xlim['xmax']) > np.abs(xmax):
                xmax = self.xlim['xmax']

        ax.set_xlim(xmin, xmax)

        # add bottom and wlev
        if structure == 'RRM' or structure == 'CRM':
            # rubble mound structure
            x_wlev = (0.5*new_design._input_arguments['B']
                      + new_design._input_arguments['slope'][1]*new_design.Rc
                            /new_design._input_arguments['slope'][0])

            wlev = new_design._LimitStates[new_design._state_overtopping].h

        elif structure == 'RC' or structure == 'CC':
            # vertical structure
            x_wlev = np.max(coordinates['caisson']['x'])
            wlev = new_design._LimitStates[
                        new_design.structure['caisson']['state_overtop']].h

        ax.hlines(
            y=wlev, xmin=xmin, xmax=-x_wlev, color='dodgerblue', lw=1)
        ax.hlines(
            y=wlev, xmin=x_wlev, xmax=xmax, color='dodgerblue', lw=1)
        ax.hlines(y=0, xmin=xmin, xmax=xmax, color='peru', lw=1, zorder=5)

        # format the ticks
        for xtick in ax.get_xticklabels():
            xtick.set_fontsize(tick_fontsize)

        for ytick in ax.get_yticklabels():
            ytick.set_fontsize(tick_fontsize)

        # format the plot with equal aspect and removing whitespace
        ax.set_aspect('equal', adjustable='box')
        ax.grid()
        fig.subplots_adjust(left=0.08, bottom=0.08, right=0.95, top=1)

        # create canvas for the figure
        canvas = FigureCanvasTkAgg(fig, frame)

        # show canvas and pack canvas
        canvas.draw()
        canvas.get_tk_widget().pack(expand=False)

def no_sliders(frame):
    """ Place a label for when no sliders have been added

    Parameters
    ----------
    frame : tk.Frame
        frame on which to pack the label
    """
    label = tk.Label(
        frame,
        text='No varying input has been specified\n Sliders not available')
    label.pack(fill=tk.Y, expand=True)

def _slider_width(parameters, frame_width):
    """ Determine the width of a slider

    Parameters
    ----------
    parameters : list
        list of parameters
    frame_width : int
        width of the frame in pixels

    Returns
    -------
    int
        required width of the slider to fill the frame
    """
    # get the longest parameter
    max_len = 0

    for parameter in parameters:
        length = len(parameter)
        if length >= max_len:
            max_len = length

    # width of a character is pixels
    char_width = 7

    # compute width of the slider
    width = frame_width - max_len*char_width - 5

    return width

def make_slider(parameter, values, frame, slider_width):
    """ Make a slider and label

    Function to make a slider and label and place them on the specified
    frame

    Parameters
    ----------
    parameter : str
        name of the parameter
    values : dict
        dictionary with the minimum and maximum value of the slider
    frame : tk.Frame
        frame on which to place the label and slider
    slider_width : int
        required width of a slider
    """
    # height of a slider
    slider_height = 42

    # create frame for the slider and label
    parent = tk.Frame(
        frame, width=frame.winfo_reqwidth(), height=slider_height)
    parent.pack_propagate(False)
    parent.pack(side=tk.TOP)

    # create frame for the label
    label_frame = tk.Frame(
        parent, width=frame.winfo_reqwidth()-slider_width,
        height=slider_height)
    label_frame.pack_propagate(False)
    label_frame.pack(side=tk.LEFT)

    # create label
    label = tk.Label(label_frame, text=parameter, padx=0)
    label.pack(side=tk.BOTTOM)

    # create frame for the slider
    slider_frame = tk.Frame(parent, width=slider_width, height=slider_height)
    slider_frame.pack_propagate(False)
    slider_frame.pack(side=tk.RIGHT)

    # determine the resolution of the slider
    resolution = 1/(10**(values['resolution']+1))

    # create slider
    slider = tk.Scale(
        slider_frame, from_=values['min'], to=values['max'],
        orient=tk.HORIZONTAL, resolution=resolution, length=slider_width)
    slider.pack()

    # return the slider
    return slider
