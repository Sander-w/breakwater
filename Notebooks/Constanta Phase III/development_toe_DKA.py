# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 11:30:15 2022

ADAPTED FROM breakwater.rubble2D.py

@author: JY6
"""

from breakwater.core.toe import toe_stability

Delta = 1.574 # To be calculated in notebook
toe_layers = 1.82 # Number of layers in the toe

t_underlayer = 0.7 #To be determined once for the project.
# Maybe not necessary when constructing on rock

t_filter = 0 #0.8 #For 60-300kg filter. 1.3 for 300-1000kg filter


# get values from the LimitState
Hs = 3.23 # To be changed in notebook
h = 6.46 # To be made case dependent in notebook
Nod = 0.5 # To be made case dependent in notebook


# make first estimate for the water level above the toe
ht_estimate = h - t_underlayer-t_filter

# use while loop since ht depends on dn50 of the toe
# first set temporary values for the while loop
Dn50_toe = 0
Dn50_toe_temp = 0
compute_toe = True
counter = 0
Dn50_toe_log = []

while compute_toe:
    Dn50_toe_computed = toe_stability(Hs, 
                                      h, 
                                      ht_estimate, 
                                      Delta, 
                                      Nod)

    # check for convergence
    if abs(Dn50_toe_computed - Dn50_toe_temp) < 0.05 or counter > 50:
        # value has converged, so break loop
        compute_toe = False

    # replace old value with the new one
    Dn50_toe_temp = Dn50_toe_computed

    # make new estimate for the water level above the toe
    toe_thickness = t_underlayer + t_filter + Dn50_toe_temp*toe_layers
    ht_estimate = h - toe_thickness
    
    counter += 1
    
    # check if computed Dn50 is larger than current normative Dn50
    if Dn50_toe_temp > Dn50_toe:
        # if larger the normative Dn50 must be changed
        Dn50_toe = Dn50_toe_temp
        Dn50_toe_log.append(Dn50_toe_temp)
        
        
        
#        state_toe = i
    
#     class_toe = Grading.get_class(Dn50_toe)
#     class_Dn50 = Grading.get_class_dn50(class_toe)


# def _estimate_htoe(self, Dn50=0):
#     """Method to estimate the height of the toe"""
#     # set ht variable
#     htoe = 0

#     # get normative layer thickness
#     for i, id in enumerate(self.variantIDs):
#         # get structure
#         structure = self.get_variant(id)

#         # get thickness of the layer
#         t_armour = self._layer_thickness(
#             "armour", self._input_arguments["armour"], structure
#         )
#         t_underlayer = self._layer_thickness("underlayer", "Rock", structure)
#         t_filter = self._layer_thickness("filter layer", "Rock", structure)

#         htoe_est = t_underlayer + t_filter

#         if Dn50 != 0:
#             htoe_est += np.ceil(t_armour / Dn50) * Dn50

#         # check if larger than previous estimate
#         if htoe_est > htoe:
#             htoe = htoe_est

#     return htoe



