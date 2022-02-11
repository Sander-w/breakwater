#%%
# Load packages
from math import nan
from unittest import result
from numpy import NaN
from shapely.geometry.polygon import Polygon
from shapely.geometry import LineString, box

from copy import deepcopy

import datetime

import geolib
from pathlib import Path

import breakwater as bw

from geo_dsuite import (
    AnalysisMethod,
    DStabilityResults,
    Headlines,
    Reflines,
    SettledStabilityAnalysis,
    SettledStabilityStage,
    Materials,
    SoilLayers,
    SoilBoundaries,
    analysis_method,
    settled_stability
)

from typing import List
from geo_dsuite import soil_layers

from geo_dsuite.soil_layers import LayerDict


#--------------------------------------------------------------------------
#%%
# Breakwater package input
from demo_2D_executable import breakwater_design

bw_configs = breakwater_design()

# Plot variant from Breakwater package
bw_configs.df.loc[len(bw_configs.df)-1].concept.plot('all')

#%%
# Define config and variant
config_idx = len(bw_configs.df)-1
print(f'Configuration no. {config_idx}')
config = bw_configs.df.loc[config_idx].concept
variant = 'a'

# Get coordinates for specific idx from breakwater package input
variant_polygons = {}
for variant_idx in config.variantIDs:
    variant_polygons[variant_idx] = config.to_polygon_coordinates(variant_idx)

# Create polygons from coordinates
for variant_idx,variant_dict in variant_polygons.items():
    for key,value in variant_dict.items():
        variant_polygons[variant_idx][key] = Polygon(value).buffer(0)

#--------------------------------------------------------------------------
#%%
# Soil parameters input
materials = Materials.from_excel(r"C:\Users\QPQ\Python\Projects\breakwater\Notebooks\Geotech/Materials_with_defaults.xlsx")
materials.index.name=None
materials_input = materials.to_soils()

#--------------------------------------------------------------------------
# %%
# Geometry setup
# Base geometry example 
geom_list = []
materials_list = []

geom_list.append(Polygon([(-200, -30), (-200, -20), (200, -20), (200, -30)]))
materials_list.append(materials_input[9].name)
geom_list.append(Polygon([(-200, -20), (-200, -10), (200, -10), (200, -20)]))
materials_list.append(materials_input[5].name)
geom_list.append(Polygon([(-200, -10), (-200, -5), (200, -5), (200, -10)]))
materials_list.append(materials_input[4].name)
geom_list.append(Polygon([(-200, -5), (-200, 0), (200, 0), (200, -5)]))
materials_list.append(materials_input[2].name)

# Append breakwater geometry and material to the respective lists
for variant_idx,variant_dict in variant_polygons.items():
    for key,value in variant_dict.items():
        if variant_idx == variant:
            geom_list.append(value)
            materials_list.append(materials.name.loc['Zand_Dijk_Za']) 

# SoilLayers setup
my_soil_layers = SoilLayers.from_base_input_list(geom_list, materials_list)
my_soil_layers = my_soil_layers.fix_geometry_consistency()

# Plot soillayers
my_soil_layers.plot(materials=materials)        

#--------------------------------------------------------------------------
#%%
# Headlines setup
hl1 = LineString([(-50, 10), ((-10/2), 10), ((10/2 + 2*10), 1), (50, 1)])
my_headlines = Headlines.from_base_input("Freatische lijn (PL1)", hl1, True)
my_headlines = my_headlines.set_index("name").sort_index()

#%%
# Reflines setup
rf1 = LineString([(-50, -5), (50, -5)])
my_reflines = Reflines.from_base_input("Referentie PL1", rf1, "Freatische lijn (PL1)", "Freatische lijn (PL1)")
my_reflines = my_reflines.set_index("name").sort_index()

#--------------------------------------------------------------------------
#%%
# Analysis method setup
# Search grid
search_grid = geolib.models.dstability.analysis.DStabilitySearchGrid(
    bottom_left=geolib.geometry.one.Point(
        label='', id=None, x=-60.0, y=-999.0, z=30, tolerance=0.0001),
    number_of_points_in_x= 30, 
    number_of_points_in_z= 20, 
    space= 1.0
    )

# Analysis method
my_analysis_method = geolib.models.dstability.analysis.DStabilityBishopBruteForceAnalysisMethod(
    search_grid=search_grid, 
    bottom_tangent_line_z= -10, 
    number_of_tangent_lines= 25, 
    space_tangent_lines= 0.5
    )

#--------------------------------------------------------------------------
#%%
# SettledStabilityStage setup
initial = SettledStabilityStage(
        name="Initial",
        headlines=my_headlines,
        reflines=my_reflines,
        soil_layers=my_soil_layers,
        analysis_method=my_analysis_method,
    )

#%%
# Stability analysis
analysis = SettledStabilityAnalysis(
        materials=materials,
        stages=[initial],
        path=Path(f"results"),
    )

#%%
# To settlement model 
m = initial.to_settlement_model(materials=materials)
m.serialize(Path(f"results/ScratchSettlementModel_ConfigNr-{config_idx}_Variant-{variant}.sli"))

#%%
# To stability model and calculate
model, mappings = analysis.to_stability_model()
model.serialize(Path(f"results/ScratchStabilityModel_ConfigNr-{config_idx}_Variant-{variant}.stix"))
model.execute()

#--------------------------------------------------------------------------
# %%
# Phased execution section

# Geometry list from breakwater package
new_geom_list = [variant_polygons[variant][x] for x in variant_polygons[variant]]

# Create soillayers to enable .fix_geometry_consistency()
soil_layer_phased = SoilLayers.from_base_input_list(new_geom_list, materials_list)
soil_layer_phased = soil_layer_phased.fix_geometry_consistency()

# Create soillayers
soil_layer_phased = SoilLayers.from_base_input_geom(soil_layer_phased.unary_union, materials.name.loc['Zand_Dijk_Za'], "20")

# Define bounds and step size
(xmin, ymin, xmax, ymax) = soil_layer_phased.unary_union.bounds
steps = (ymax-ymin)
steps_list = [1 for x in range(int(round(steps,0)))]
steps_list[-1] = -1 # Last step always to the top of the geometry with "-1"

# Split layer
new_soil_layer_phased = soil_layer_phased.layer_in_stages("20", steps_list)

# Iterate over split layers to check for MultiPolygons
new_soil_layer_phased = SoilLayers.check_multipolygons(new_soil_layer_phased)

# Resort layers
new_soil_layer_dataframe = SoilLayers(new_soil_layer_phased).vertical_sort()
soil_layer_loading = deepcopy(new_soil_layer_dataframe.iloc[:-1])
# Plot soillayers
new_soil_layer_dataframe.plot(materials=materials)

#%%
# Initialize phased stages
initial = SettledStabilityStage(
        name="Initial",
        headlines=my_headlines,
        reflines=my_reflines,
        soil_layers=soil_layer_loading,
        analysis_method=my_analysis_method,
    )

analysis = SettledStabilityAnalysis(
        materials=materials,
        stages=[initial],
        path=Path(f"results"),
    )

# To phased stability file
model, mappings = analysis.to_stability_model()
model.serialize(Path(f"results/ScratchStabilityModel_ConfigNr-{config_idx}_Variant-{variant}_Phased.stix"))
# %%

# Settlement model setup
settle_model = initial.to_settlement_model(materials=materials)

# %%
# Apply load to settlement model from last layer in phased execution

# Create last layer as load
new_load = deepcopy(new_soil_layer_dataframe.geometry[-1])
new_load.geometry = soil_layers.top_boundary_from_polygon(new_load)
new_load_points = soil_layers._skip_consecutive_duplicates([(geolib.geometry.one.Point(x=x,z=z)) for x,z in new_load.geometry.coords])

# Add non uniform load to settlement model
settle_model.add_non_uniform_load('Last load', new_load_points, datetime.timedelta(days=100), 18, 20)

# To settlement file
settle_model.serialize(Path(f"results/ScratchSettlementModel_ConfigNr-{config_idx}_Variant-{variant}_Phased.sli"))


# %%
#--------------------------------------------------------------------------
# Miscellaneous methods
# def skip_duplicates(elements):
#         prev = None
#         for elem in elements:
#             if elem != prev:
#                 yield elem
#                 prev = elem


# def sort_clockwise(coords: List, clockwise=True):
#     """Sort coordinate list clockwise"""

#     if clockwise:
#         # compute centroid
#         cent=(sum([p[0] for p in coords])/len(coords),sum([p[1] for p in coords])/len(coords))
#         # sort by polar angle
#         coords.sort(key=lambda p: -math.atan2(p[1]-cent[1],p[0]-cent[0]))

#     else:
#         # compute centroid
#         cent=(sum([p[0] for p in coords])/len(coords),sum([p[1] for p in coords])/len(coords))
#         # sort by polar angle
#         coords.sort(key=lambda p: math.atan2(p[1]-cent[1],p[0]-cent[0]))
    
#     return coords
