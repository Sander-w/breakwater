"""
Breakwater Design with Python
"""

__version__ = '1.1'
__author__ = 'S. Winkel'

# import hydraulic conditions
from .core.battjes import BattjesGroenendijk
from .core.goda import goda_wave_heights
from .utils.wave import shoaling_coefficient, dispersion
from .conditions import LimitState

# import materials
from .material import RockGrading, Xbloc, XblocPlus, ConcreteArmour

# import breakwaters
from .design import Configurations
from .rubble import RubbleMound, ConcreteRubbleMound, RockRubbleMound
from .caisson import Caisson

# file loader
from .design import read_breakwaters

# excel input
from .material import read_grading, read_units
from .design import read_configurations
from .utils.input_generator import generate_excel

# interactive design tool (tkinter app)
from .interactive import interactive_design

# import the soil
from .core.soil import Soil
