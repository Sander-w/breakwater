import sys
import os

# get path of parent directory
path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# insert path to the sys path
sys.path.insert(0, path)
