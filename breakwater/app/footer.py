import tkinter as tk
from webbrowser import open_new

from .settings import settings

def footer(frame):
    """ Add footer

    Function to add a footer with text and a link to a LinkedIn page

    Parameters
    ----------
    frame : tk.Frame
        frame to which the footer must be added
    """
    # add footer
    footer = tk.Label(
        frame, text='Developed by S. Winkel. Breakwater is licensed under CC BY-NC-SA 4.0', bd=1, relief=tk.SUNKEN,
        anchor='w')
    footer.pack(side=tk.BOTTOM, fill=tk.X)
    footer.bind(
        "<Button-1>",
        lambda e: open_new('https://www.linkedin.com/in/sander-winkel/'))
