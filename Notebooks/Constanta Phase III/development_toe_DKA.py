# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 11:30:15 2022

ADAPTED FROM breakwater.rubble2D.py

@author: JY6
"""

from breakwater.core.toe import toe_stability
from development_material_DKA import get_class, get_Dn50
import pandas as pd
import matplotlib.pyplot as plt

#%% Input some data for testing purposes
gradings_data = pd.read_excel(r"C:\Users\JY6\github\breakwater\Notebooks\Constanta Phase III\Input data\test_data_phase_II.xlsx", 
                                      sheet_name='input_rock_gradings',
                                      index_col = 0,
                                      skiprows = 2)

# Values from excel
Delta = 1.574 # To be calculated in notebook
toe_layers = 1.82 # Number of layers in the toe
t_underlayer = 0.7  #To be determined once for the project.
                    # Maybe not necessary when constructing on rock
t_filter = 0 #0.8 #For 60-300kg filter. 1.3 for 300-1000kg filter
rho_a = 2600

# get values from the LimitState
Hs = 3.63 # To be changed in notebook
h = 6 # To be made case dependent in notebook
Nod = 0.5 # To be made case dependent in notebook


#%% Actual calculation


# use while loop since ht depends on Dn50 of the toe
# first set temporary values for the while loop
Dn50_toe = 0
Dn50_toe_temp = 0
compute_toe = True
counter = 0
Dn50_toe_selected_log = []
Dn50_toe_computed_log = []

toe_layer_thickness = t_underlayer + t_filter + Dn50_toe*toe_layers
ht_estimate = h - toe_layer_thickness

while compute_toe:
    Dn50_toe_computed = toe_stability(Hs, 
                                      h, 
                                      ht_estimate, 
                                      Delta, 
                                      Nod)    
    Dn50_toe_computed_log.append(Dn50_toe_computed)
    
    
    # check for convergence
    if Dn50_toe_computed - Dn50_toe_temp < 0 or counter > 50:
        # value has converged, so break loop
        compute_toe = False

    counter += 1
    
    # check if computed Dn50 is larger than current normative Dn50
    if Dn50_toe_computed > Dn50_toe:
        # if larger the normative Dn50 must be changed
        Dn50_toe = Dn50_toe_computed
        
    # Set Dn50 toe class average
    class_toe = get_class(Dn50_toe_computed, rho_a, gradings_data)
    Dn50_toe_computed = get_Dn50(class_toe, gradings_data)
    
    # Check if filter layer needs to be applied
    if gradings_data.at[class_toe, 'Toe underlayer'] in gradings_data.index:
        class_filter = gradings_data.at[class_toe, 'Toe underlayer']
        Dn50_filter = get_Dn50(class_filter, gradings_data)
        t_filter = Dn50_filter*toe_layers
    
    # replace old values with the new ones
    Dn50_toe_temp = Dn50_toe_computed
    toe_layer_thickness = t_underlayer + t_filter + Dn50_toe_temp*toe_layers
    ht_estimate = h - toe_layer_thickness
    
    Dn50_toe_selected_log.append(Dn50_toe_temp)
    

#%% Plot example of figures

X = list(range(0,len(Dn50_toe_selected_log)))

plt.clf()
fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.plot(X, Dn50_toe_selected_log, label = 'Dn50_toe_selected')
ax.plot(X, Dn50_toe_computed_log, label = 'Dn50_toe_compute')
ax.set(xlabel = 'Iteration',
            ylabel = 'Dn50')
ax.legend()
ax.grid('on')
        

