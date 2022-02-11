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

# Define variant to be calculated
variant = 'b'

#%%
# Miscellaneous methods
def skip_duplicates(elements):
        prev = None
        for elem in elements:
            if elem != prev:
                yield elem
                prev = elem

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

#--------------------------------------------------------------------------

# Breakwater package input
battjes = bw.BattjesGroenendijk(Hm0=4.4, h=15, slope_foreshore=(1,100))
H2_per = battjes.get_Hp(0.02)
ULS = bw.LimitState(
    h=15, Hs=4.5, Hm0=4.4, H2_per=H2_per, Tp=9.4, Tm=8.8, T_m_min_1=9.7,
    Sd=5, Nod=2, q=20, label='ULS')

rho_rock = 2650

NEN = bw.RockGrading(rho=rho_rock)
M50_core = (1+1000)/2
Dn50_core = (M50_core / rho_rock)**(1/3)

RRM = bw.RockRubbleMound(
    slope=(1,2), slope_foreshore=(1,100), rho_w=1025, B=5.5, N=2100,
    LimitState=ULS, Grading=NEN, Dn50_core=0.4, core_material={'class': 'QR 1-1000kg', 'Dn50': Dn50_core})

variant_polygons = {}
for variant_idx in RRM.variantIDs:
    variant_polygons[variant_idx] = RRM.to_polygon_coordinates(variant_idx)

for variant_idx,variant_dict in variant_polygons.items():
    for key,value in variant_dict.items():
        variant_polygons[variant_idx][key] = Polygon(value)

# Plot variant from Breakwater package
RRM.plot(variant)

#%%
# Soil parameters input
materials = Materials.from_excel("Materials_with_defaults.xlsx")
materials.index.name=None
materials_input = materials.to_soils()
materials_list = []

#--------------------------------------------------------------------------
# %%
# Geometry setup
# Base geometry example 
geom_list = []

poly1_point1 = -100, -30
poly1_point2 = -100, -20
poly1_point3 = 50, -20
poly1_point4 = 50, -30

poly2_point1 = -100, -20
poly2_point2 = -100, -10
poly2_point3 = 50, -10
poly2_point4 = 50, -20

poly3_point1 = -100, -10
poly3_point2 = -100, -5
poly3_point3 = 50, -5
poly3_point4 = 50, -10

poly4_point1 = -100, -5
poly4_point2 = -100, 0
poly4_point3 = 50, 0
poly4_point4 = 50, -5

poly1_pointlist = poly1_point1, poly1_point2, poly1_point3, poly1_point4
poly2_pointlist = poly2_point1, poly2_point2, poly2_point3, poly2_point4
poly3_pointlist = poly3_point1, poly3_point2, poly3_point3, poly3_point4
poly4_pointlist = poly4_point1, poly4_point2, poly4_point3, poly4_point4

poly1 = Polygon([[p[0], p[1]] for p in poly1_pointlist])
poly2 = Polygon([[p[0], p[1]] for p in poly2_pointlist])
poly3 = Polygon([[p[0], p[1]] for p in poly3_pointlist])
poly4 = Polygon([[p[0], p[1]] for p in poly4_pointlist])

geom_list.append(poly1)
materials_list.append(materials_input[9].name)

geom_list.append(poly2)
materials_list.append(materials_input[5].name)

geom_list.append(poly3)
materials_list.append(materials_input[4].name)

geom_list.append(poly4)
materials_list.append(materials_input[2].name)

# Breakwater geometry from breakwater package
for variant_idx,variant_dict in variant_polygons.items():
    for key,value in variant_dict.items():
        if variant_idx == variant:
            geom_list.append(value)
            materials_list.append(materials_input[1].name)
            
#--------------------------------------------------------------------------
#%%
# Headlines setup
hl1_point1 = -50, 10
hl1_point2 = -10/2, 10
hl1_point3 = (10/2 + 2*10), 1
hl1_point4 = 50, 1

line1_pointlist = hl1_point1, hl1_point2, hl1_point3, hl1_point4
hl1 = LineString([[p[0], p[1]] for p in line1_pointlist])

head_line_data = [{
    "name": "Freatische lijn (PL1)",  
    "geometry": hl1,
    "is_phreatic": True,
    }
]

my_headlines = Headlines.from_base_headlines(head_line_data)
my_headlines = my_headlines.set_index("name").sort_index()

#%%
# Reflines setup
rf1_point1 = -50, -5
rf1_point2 = 50, -5

rfline1_pointlist = rf1_point1, rf1_point2
rf1 = LineString([[p[0], p[1]] for p in rfline1_pointlist])

rf_line_data = [{
    "name": "Referentie PL1",
    "geometry": rf1,
    "headline_top": "Freatische lijn (PL1)",
    "headline_bottom": "Freatische lijn (PL1)"
    }
]

my_reflines = Reflines.from_base_reflines(rf_line_data)
my_reflines = my_reflines.set_index("name").sort_index()


#--------------------------------------------------------------------------
#%%
# Analysis method setup
# Search grid
bottom_left = geolib.geometry.one.Point(label='', id=None, x=-45.0, y=-999.0, z=12, tolerance=0.0001)
number_of_points_in_x = 30
number_of_points_in_z = 20
search_grid = geolib.models.dstability.analysis.DStabilitySearchGrid(bottom_left=bottom_left,number_of_points_in_x=number_of_points_in_x, number_of_points_in_z=number_of_points_in_z, space=1.0)

# Tangent lines
bottom_tangent_line_z = -10
number_of_tangent_lines = 25
space_tangent_lines = 0.5

# Analysis method
my_analysis_method = geolib.models.dstability.analysis.DStabilityBishopBruteForceAnalysisMethod(search_grid=search_grid, bottom_tangent_line_z=bottom_tangent_line_z, number_of_tangent_lines=number_of_tangent_lines, space_tangent_lines=space_tangent_lines)

#--------------------------------------------------------------------------


# TODO Create standard breakwater input method/class
# TODO Create breakwater import method


# SoilLayers setup
layer_data = []

for idx, geom in enumerate(geom_list):
    layer_data.append({
            "name": f"{idx}",
            "notes": "",
            "geometry": geom,
            "material": materials_list[idx]
        }
    )

my_soil_layers = SoilLayers.from_base_layers(layer_data)
my_soil_layers = my_soil_layers.fix_geometry_consistency()

# Plot soillayers
my_soil_layers.plot(materials=materials)

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
# To settlement model 
m = initial.to_settlement_model(materials=materials)
m.serialize(Path(f"results/ScratchSettlementModelVariant{variant}.sli"))


#%%
# Stability analysis
analysis = SettledStabilityAnalysis(
        materials=materials,
        stages=[initial],
        path=Path(f"results"),
    )

#%%
# To stability model
model, mappings = analysis.to_stability_model()
model.serialize(Path(f"results/ScratchStabilityModelVariant{variant}.stix"))
#--------------------------------------------------------------------------

# %%
# Calculate stability model
model.execute()

#--------------------------------------------------------------------------
# %%
# Phased execution section

# Geometry list from breakwater package
new_geom_list = [variant_polygons[variant][x] for x in variant_polygons[variant]]

# Soil layer setup
# Geometry dict setup
phased_layer_data = []
for idx, geom in enumerate(new_geom_list):
    phased_layer_data.append({
            "name": f"{idx}",
            "notes": "",
            "geometry": geom,
            "material": materials_list[idx]
        }
    )

# Create soillayers to enable .fix_geometry_consistency()
soil_layer_phased = SoilLayers.from_base_layers(phased_layer_data)
soil_layer_phased = soil_layer_phased.fix_geometry_consistency()

# TODO Fix underlying layer splitting in a single method

# Merge all applicable layers
combined = [x for x in soil_layer_phased.geometry]
combined_geo = combined[0]
[combined_geo := combined_geo.union(x) for x in combined]

# Geometry list setup (2)
phased_layer_data = [{
            "name": "20",
            "notes": "",
            "geometry": combined_geo,
            "material": materials_list[-1]
        }]

# Create soillayers
soil_layer_phased = SoilLayers.from_base_layers(phased_layer_data)
soil_layer_phased = soil_layer_phased.fix_geometry_consistency()

# Split merged geometry into phased execution
#%%
# Split merged geometry into phased execution

# Define bounds and step size
(xmin, ymin, xmax, ymax) = combined_geo.bounds
steps = (ymax-ymin)
steps_list = [1 for x in range(int(round(steps,0)))]
steps_list[-1] = -1 # Last step always to the top of the geometry with "-1"

# Split layer
new_soil_layer_phased = soil_layer_phased.layer_in_stages("20", steps_list)

# Iterate over split layers to check for MultiPolygons
pop_list = []
for idx, geom in enumerate(new_soil_layer_phased):
    if type(new_soil_layer_phased[idx].geometry) != Polygon:
        # When found, split Multipolygons into Polygons with .buffer(0)
        multipoly = new_soil_layer_phased[idx].geometry.buffer(0) 
        for i, subgeo in enumerate(multipoly):
            
            stage = deepcopy(new_soil_layer_phased[idx])
            stage.name = new_soil_layer_phased[idx].name + f".{i+1}"
            stage.notes = ""
            stage.geometry = subgeo
            stage.material = materials_list[-1]

            # Append to geometry
            new_soil_layer_phased.append(stage)

        # Keep track of MultiPolygon layers
        pop_list.append(idx)

# Remove Multipolygon layers
[new_soil_layer_phased.pop(x) for x in pop_list]

# Resort layers
new_soil_layer_dataframe = SoilLayers(new_soil_layer_phased).vertical_sort()
# Plot soillayers
new_soil_layer_dataframe.plot(materials=materials)

#%%
# Initialize phased stages
initial = SettledStabilityStage(
        name="Initial",
        headlines=my_headlines,
        reflines=my_reflines,
        soil_layers=new_soil_layer_dataframe,
        analysis_method=my_analysis_method,
    )

analysis = SettledStabilityAnalysis(
        materials=materials,
        stages=[initial],
        path=Path(f"results"),
    )

# To phased stability file
model, mappings = analysis.to_stability_model()
model.serialize(Path(f"results/ScratchStabilityModelVariant{variant}_Phased.stix"))
# %%

# %%
# Apply load to settlement model from last layer in phased execution
# Layers setup without last layer
new_load = deepcopy(new_soil_layer_dataframe.geometry[-1])
soil_layer_loading = deepcopy(new_soil_layer_dataframe.iloc[:-1])


# Stage setup without last layer
initial = SettledStabilityStage(
        name="Initial",
        headlines=my_headlines,
        reflines=my_reflines,
        soil_layers=soil_layer_loading,
        analysis_method=my_analysis_method,
    )

# Settlement model setup
settle_model = initial.to_settlement_model(materials=materials)

# Top boundary for load
new_load.geometry = soil_layers.top_boundary_from_polygon(new_load)
new_load_points = skip_duplicates([(geolib.geometry.one.Point(x=x,z=z)) for x,z in new_load.geometry.coords])

# Add non uniform load to settlement model
settle_model.add_non_uniform_load('Last load', new_load_points, datetime.timedelta(days=100), 18, 20)

# To settlement file
settle_model.serialize(Path(f"results/ScratchSettlementModelVariant{variant}_Phased.sli"))

# %%


# polys = configs.df.loc[99].concept.to_polygon_coordinates('a')

# from shapely.geometry import Polygon
# variant_polygons = {}
# for variant_idx in CRMR.variantIDs:
#     variant_polygons[variant_idx] = configs.df.loc[99].concept.to_polygon_coordinates('a')

# for variant_idx,variant_dict in variant_polygons.items():
#     for key,value in variant_dict.items():
#         variant_polygons[variant_idx][key] = Polygon(value)

# from shapely.geometry import Polygon
# variant_polygons = {}
# for variant_idx in CRMR.variantIDs:
#     variant_polygons[variant_idx] = configs.df.loc[99].concept.to_polygon_coordinates('a')

# for variant_idx,variant_dict in variant_polygons.items():
#     for key,value in variant_dict.items():
#         variant_polygons[variant_idx][key] = Polygon(value)