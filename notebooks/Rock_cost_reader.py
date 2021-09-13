import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import filedialog

root = tk.Tk()
root.withdraw()

file_path = filedialog.askopenfilename()

if len(file_path) != 0:
    file = pd.read_excel(file_path)
