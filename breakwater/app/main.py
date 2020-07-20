import tkinter as tk
import os

from .settings import settings, size
from .input import ParameterInput, VaryingInput
from .interactive import InteractiveDesign
from .footer import footer


class BreakwaterDesign(tk.Tk):
    """ Tkinter app for interactive design

    Tkinter application to design breakwaters

    Parameters
    ----------
    vkwargs : dict
        dictionary of the valid arguments for the specified structures
    input : dict
        dictionary with the Python input, this is at least a LimitState
        and Grading.
    to_design : dict
        dictionary with the structures to design, key is the structure,
        value is a bool
    cost : dict
        dictionary with the price of concrete and sand
    display_warnings : bool
        True if warnings must be displayed when designing, False if not
    """

    def __init__(
            self, vkwargs, input, to_design, cost, display_warnings,
            *args, **kwargs):
        """ See help(BreakwaterDesign) for more info """
        # initialize tkinter
        tk.Tk.__init__(self, *args, **kwargs)

        # set attribute
        self.vkwargs = vkwargs
        self.parameters = input
        self.to_design = to_design
        self.cost = cost
        self.display_warnings = display_warnings

        # set attribute for the previous page of InteractiveDesign
        self.previous = 'ParameterInput'

        # determine number of concepts and set as attribute
        self.num_concepts = to_design['RM']['num'] + to_design['C']['num']

        # set icon of the app
        file = os.path.dirname(os.path.abspath(__file__))
        icon_path = os.path.join(file, 'bw-icon.ico')
        try:
            self.iconbitmap(icon_path)
        except:
            print(f'Could not load the icon from {icon_path}')

        # create a containter
        container = tk.Frame(self)

        # format main lay-out
        container.pack(side='top', fill='both', expand=True)
        container.rowconfigure(0, weight=1)
        container.columnconfigure(0, weight=1)

        # add title
        self.title('Interactive Breakwater Design')

        # create dict to store frames in
        self.frames = {}

        # add frames
        for Frame in (StartPage, ParameterInput, InteractiveDesign, VaryingInput):
            frame = Frame(container, self)
            self.frames[Frame] = frame
            frame.grid(row=0, column=0, sticky='nsew')

        # set and show home page (=start page)
        self.home_page = StartPage
        self.home()

    def show_frame(self, controller, default=True):
        """ Load a new frame

        Parameters
        ----------
        controller : tk.Frame
            The frame to show
        default : bool, optional, default: True
            True if the default geometry settings must be used, False
            if the adapted geometry settings must be showed
        """
        if default:
            # get width and height
            width = settings('lay-out')['width']
            height = settings('lay-out')['height']
        else:
            # position app in the center of the screen
            width, height = size(
                self.to_design, 'interactive', self.num_concepts)

        ws = self.winfo_screenwidth()
        hs = self.winfo_screenheight()
        x = (ws/2) - (width/2)
        y = (hs/2) - (height/2)
        self.geometry(f'{width}x{height}+{round(x)}+{round(y)}')

        # get frame
        frame = self.frames[controller]

        # display the frame
        frame.tkraise()

    def home(self):
        """ Go to the home page """
        self.show_frame(self.home_page)

    def back(self, page):
        """ Go back to the previous page

        Method to go back to the previous page from InteractiveDesign

        Parameters
        ----------
        page : str
            name of the page
        """
        if page == 'VaryingInput':
            self.show_frame(VaryingInput)

        elif page == 'ParameterInput':
            self.show_frame(ParameterInput)

        else:
            raise NotImplementedError(f'{page} is an invalid page name')


class StartPage(tk.Frame):
    """ Start page of the app """

    def __init__(self, parent, controller):
        """ see help(StartPage) for more info """
        # create the frame
        tk.Frame.__init__(self, parent)

        # get width and height
        width = settings('lay-out')['width']
        height = settings('lay-out')['height']

        # set height of frames
        top_height = settings('lay-out')['top_height']
        bottom_height = settings('lay-out')['bottom_height']
        center_height = height - top_height - bottom_height - 50

        # create lay out
        top_frame = tk.Frame(self, width=width, height=top_height)
        top_frame.pack_propagate(False)
        center_frame = tk.Frame(self, width=width, height=center_height)
        center_frame.pack_propagate(False)
        bottom_frame = tk.Frame(self, width=width, height=bottom_height)
        bottom_frame.pack_propagate(False)

        top_frame.pack(side=tk.TOP, pady=5)
        center_frame.pack()
        bottom_frame.pack(side=tk.BOTTOM)

        # add title of the page
        label = tk.Label(
            center_frame,
            text='Welcome to the Interactive Design Application of breakwater',
            **settings('title'))
        label.place(relx=0.5, rely=0.25, anchor=tk.CENTER)

        # welcome message
        label = tk.Label(
            center_frame,
            text=(''
            'All parameters required for the design can be given on the '
            'next page. On this page you can also specify \n which parameters'
            ' are allowed to vary, for each of these parameters a slider is '
            'generated. The minimum and \n maximum values for these sliders '
            'can be specified on the Varying Input page, which is '
            'accessed by checking \n the checkbox behind the parameter and '
            'pressing Start Design. \n'
            '\n'
            'When all required input has been given the Interactive Design '
            'page can be accessed, on this page a \n cross-section of the '
            'specified breakwater structure(s) is presented together with '
            'sliders if parameters \n are allowed to vary. The design can be '
            'updated by changing the sliders and pressing the \n Update '
            'Design button. A new cross-section is then plotted, with the '
            'previous design in grey.')
        )
        label.place(relx=0.5, rely=0.45, anchor=tk.CENTER)

        # add the footer
        footer(bottom_frame)

        # add next button
        nav_settings = settings('navigation')

        next_button = tk.Button(
            bottom_frame, text='Next', height=nav_settings['height'],
            width=nav_settings['width'],
            command=lambda: controller.show_frame(ParameterInput))
        next_button.pack(
            side=tk.RIGHT, padx=nav_settings['padx'], pady=nav_settings['pady'])
