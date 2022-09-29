from pykml import parser
import numpy as np
from shapely.geometry import Polygon
from shapely.geometry.polygon import orient
import math
from math import radians, cos, sin, asin, sqrt
import string

from breakwater.conditions import LimitState
from breakwater.utils.exceptions import InputError
from breakwater.rubble_2D import RockRubbleMound, ConcreteRubbleMound, ConcreteRubbleMoundRevetment


def coords_from_kml(kml_file):

    """

    Parameters
    ----------
    kml_file: str
        string with the path of the kml_file
    root_location: loc
        location within the kml file with the coordinates
    Returns
    -------
    tuple
    """

    root = parser.fromstring(open(kml_file, 'rb').read())
    coords = str(root.Document.Folder.Placemark.LineString.coordinates).strip()

    coords = coords.split(' ')

    coord_tuples = []

    for c in coords:
        c_split = c.split(',')[:-1]
        coord_tuples.append((float(c_split[0]), float(c_split[1])))



    return coord_tuples

def create_shape(kml_file):

    """

    Parameters
    ----------
    coordinates: list
        list containing tuples with (lon, lat) of the LinearRing
    Returns
    -------
    shapely.geometry
    """

    coords = coords_from_kml(kml_file)

    if len(coords) >= 2:
        shape = Polygon(coords)
        #counterclock wise orientated
        shape_counterclock = orient(polygon= shape, sign= 1)

    else:
        xmed = (coords[0][0] + coords[1][0]) / 2
        ymed = (coords[0][1] + coords[1][1]) / 2
        coords.append([xmed, ymed])
        shape = Polygon(coords)
        #counterclock wise orientated
        shape_counterclock = orient(polygon= shape, sign= 1)

    return shape_counterclock

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance in meters between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles. Determines return value units.
    return (c * r) * 1000

def structure_orientation(counterclock_coords, wave_direction= 'right'):
    """
    Coordinates should be given in counter clockwise direction.
    This function gives the orientation of the different parts of the structure
    with respect to the North.

    Parameters
    ----------
    coords: list
        coordinates in counter clockwise direciton
    wave_direction: str
        Either left or right. Default is right. When walking in
        counter clockwise direction from point to point, from where
        are the waves coming (which direction is un-sheltered).

    Returns
    -------

    """

    orientation_dicts = {}
    r = [0]
    s = [0]
    c = [0]
    for i in range(len(counterclock_coords)-1):
        coords_anticlock = counterclock_coords[i:i+2]
        long2, long1 = coords_anticlock[1][0], coords_anticlock[0][0]
        lat2, lat1 = coords_anticlock[1][1], coords_anticlock[0][1]

        dist = math.ceil(haversine(lon1= long1, lat1= lat1, lon2= long2, lat2= lat2))

        dLon = (long2 - long1)
        x = math.cos(math.radians(lat2)) * math.sin(math.radians(dLon))
        y = math.cos(math.radians(lat1)) * math.sin(math.radians(lat2)) - math.sin(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.cos(math.radians(dLon))
        bearing = math.atan2(x,y)   # use atan2 to determine the quadrant
        bearing = math.degrees(bearing)


        if bearing < 0:
            bearing = 360 + bearing

        if 0 <= bearing <= 90 or 180 <= bearing <= 270:
                orientation = bearing + 90
        else:
            if bearing > 270:
                orientation = bearing - 270
            else:
                orientation = bearing + 90

        if wave_direction != 'right':
            if orientation >= 180:
                orientation -= 180
            else:
                orientation += 180

        orientation = int(orientation)
        letter = list(string.ascii_uppercase)

        if dist < 10:
            if f'section {letter[r[-1]]}' not in orientation_dicts.keys():
                orientation_dicts[f'section {letter[r[-1]]}'] = {'coordinates': [[(long1, lat1), (long2, lat2)]], 'orientation': [orientation], 'distance': [dist]}
                s.append(r[-1]+1)
            else:
                orientation_dicts[f'section {letter[r[-1]]}']['coordinates'].append([(long1, lat1), (long2, lat2)])
                orientation_dicts[f'section {letter[r[-1]]}']['orientation'].append(orientation)
                orientation_dicts[f'section {letter[r[-1]]}']['distance'].append(dist)
        else:
            orientation_dicts[f'section {letter[s[-1]]}'] = {'coordinates': [[(long1, lat1), (long2, lat2)]], 'orientation': [orientation], 'distance': [dist]}
            r.append(s[-1]+1)
            s.append(r[-1])

    return orientation_dicts

def wave_angles_structure(kml_path, wave_conditions, LimitStates= None, wave_direction = 'right', shape = 'LinearRing'):

    """
    Determines all the angles of the wave with respect to the orientation of all
    the parts of the design

    Parameters
    ----------
    kml_path: str
        location of the kml file
    wave_conditions: dict
        key is the wave direction and the items a dict with Hm0 and Tp
    LimitStates: lst
        list with LimitState objects
    wave_directions: dict
        per wave direction the significant wave height and wave period is provided
    shape: {'LinearRing', 'Polygon'}
        can be either LinearRing or Polygon. Use Polygon when the shape is closed
    Returns
    -------
    dict
    """

    original_coords = coords_from_kml(kml_path)
    shape_counterclock = create_shape(kml_path)

    counterclock_coords = shape_counterclock.exterior.coords

    if shape != 'Polygon':
        # delete duplicates created by polygon when the shape is a Linestring
        counterclock_coords = dict.fromkeys(counterclock_coords)
        # if there where only two coordinates, a fictive one was added to create a polygon. Delete this one
        counterclock_coords = [c for c in counterclock_coords if c in original_coords]


    orientation_dict = structure_orientation(counterclock_coords= counterclock_coords, wave_direction= wave_direction)

    for key, items in orientation_dict.items():
        orientation_dict[key]['wave_conditions'] = []
        for wave_angle, wavespecs in wave_conditions.items():
            beta = np.array(items['orientation']) - wave_angle
            beta = abs(beta)

            beta = np.where(beta >= 270, 360 - beta, beta)
            beta = np.where(beta <= 90, beta, None)
            if LimitStates == None:
                orientation_dict[key]['wave_conditions'].append({'Hm0': wavespecs['Hm0'], 'Tp': wavespecs['Tp'], 'beta': beta})
            else:
                orientation_dict[key]['wave_conditions'].append({'Hm0': wavespecs['Hm0'], 'Tp': wavespecs['Tp'], 'beta': beta, 'LimitState': LimitStates[wave_angle]})

    return orientation_dict

def all_designs(
        kml_path,
        wave_conditions,
        shape,
        slope,
        slope_foreshore,
        B,
        N,
        rho_w,
        Grading,
        core_material,
        LimitState,
        structure_type= 'breakwater',
        ArmourUnit= None,
        wave_direction= 'right',
        safety=1,
        slope_toe=(2, 3),
        B_toe=None,
        layers=1,
        layers_underlayer=2,
        vdm= 'max',
        filter_rule=None,
        Soil=None,
        phi=40,
        id=None,
        **kwargs):
    """

    Parameters
    ----------
    kml_path: str
            path with the kml file
    shape: str
        shape of the structure. default is LinearRing
    wave_conditions: dict
        dictionary with orientation the key and item is a dict with Hm0 and Tp
    wave_direction: str, default= 'right'
        When walking in counter-clockwise direction along the coordinates, do the waves
        come from the right or the left
    structure_type: str, default= 'breakwater'
        type of structure, either breakwater or revetment
    h: float
        Water depth at the structure
    Sd: int
        Allowed damage level of the armour
    Nod: int
        Allowed damage level at the toe
    q: int
        Allowed number of overtopping in q/l/m
    slope : tuple
        Slope of the armour layer (V, H). For example a slope of 3V:4H
        is defined as (3, 4)
    slope_foreshore : tuple
        slope of the foreshore (V, H). For example a slope of 1:100 is
        defined as (1, 100)
    rho_w : float
        density of water [kg/m³]
    B : float
        Crest width [m]
    N : int
        Number of incident waves at the toe of the structure [-]
    LimitState : :py:class:`LimitState` or list of :py:class:`LimitState`
        ULS, SLS or another limit state defined with
        :py:class:`LimitState`
    Grading : :py:class:`RockGrading`
        standard rock grading defined in the NEN-EN 13383-1 or a user
        defined rock grading
    core_material : dict
        core_class and Dn50
    safety : float, optional, default: 1
        safety factor of design (number of standard deviations from the
        mean)
    slope_toe : tuple, optional, default: (2,3)
        slope of the toe
    B_toe : float, optional, default: None
        width of the top of the toe in meters. By default the width of
        toe is taken as 3 * Dn50_toe.
    beta : float, optional, default: 0
        angle between direction of wave approach and a line normal to
        the breakwater (degrees).
    layers : int, optional, default: 2
        number of layers in the armour layer
    layers_underlayer : int, optional, default: 2
        number of layers in the underlayer
    vdm : {min, max, avg}, optional, default: max
        value to return in case both the deep and shallow water formula
        are valid. min for the lowest value, max for the highest value
        and avg for the average value, default is max.
    Soil : :py:class:`Soil`, optional, default: None
        by default Soil is None, which means that the geotechnical checks
        are not performed. By specifying a Soil object, the geotechnical
        checks are automatically performed.
    phi : float, optional, default: 40
        internal friction angle of rock [degrees]
    id : int, optional, default: None
        add a unique id to the breakwater

    Returns
    -------
    dict
    """
    orientation_dict = wave_angles_structure(kml_path= kml_path, wave_conditions= wave_conditions, LimitStates= LimitState,
                                             wave_direction = wave_direction, shape = shape)
    for part, specs in orientation_dict.items():
        for i in range(len(specs['wave_conditions'])):
            condition = specs['wave_conditions'][i]
            orientation_dict[part]['wave_conditions'][i]['design'] = []
            for beta in condition['beta']:
                if beta != None:
                    if structure_type == 'breakwater' and ArmourUnit == None:
                        design = RockRubbleMound(
                                                slope= slope,
                                                slope_foreshore= slope_foreshore,
                                                rho_w= rho_w,
                                                B= B,
                                                N= N,
                                                LimitState= condition['LimitState'],
                                                Grading= Grading,
                                                core_material= core_material,
                                                safety= safety,
                                                slope_toe= slope_toe,
                                                structure_type= "breakwater",
                                                B_toe= B,
                                                beta= beta,
                                                layers= layers,
                                                layers_underlayer= layers_underlayer,
                                                vdm= vdm,
                                                Soil= Soil,
                                                phi= phi,
                                                id= id,
                                                **kwargs,
                                                )

                    elif structure_type == 'breakwater' and ArmourUnit != None:
                        design = ConcreteRubbleMound(
                                                    slope= slope,
                                                    slope_foreshore= slope,
                                                    B= B,
                                                    rho_w= rho_w,
                                                    LimitState= condition['LimitState'],
                                                    ArmourUnit= ArmourUnit,
                                                    Grading= Grading,
                                                    core_material= core_material,
                                                    safety= safety,
                                                    slope_toe= slope_toe,
                                                    structure_type= "breakwater",
                                                    B_toe= B_toe,
                                                    beta= beta,
                                                    layers= layers,
                                                    layers_underlayer= layers_underlayer,
                                                    filter_rule= filter_rule,
                                                    Soil= Soil,
                                                    phi= phi,
                                                    id= id,
                                                    **kwargs,
                                                    )

                    elif structure_type == 'revetment' and ArmourUnit != None:
                        design = ConcreteRubbleMoundRevetment(
                                                            slope= slope,
                                                            slope_foreshore= slope_foreshore,
                                                            B= B,
                                                            rho_w= rho_w,
                                                            LimitState= condition['LimitState'],
                                                            ArmourUnit= ArmourUnit,
                                                            Grading= Grading,
                                                            core_material= core_material,
                                                            safety= safety,
                                                            slope_toe= slope_toe,
                                                            structure_type= 'revetment',
                                                            B_toe= B_toe,
                                                            beta= beta,
                                                            layers= layers,
                                                            layers_underlayer= layers_underlayer,
                                                            filter_rule= filter_rule,
                                                            Soil= Soil,
                                                            phi= phi,
                                                            id= id
                                                            )
                    else:
                        raise InputError(f'The combination of {structure_type} with the material is not possible.')

                    orientation_dict[part]['wave_conditions'][i]['design'].append(design)

                else:
                    orientation_dict[part]['wave_conditions'][i]['design'].append(None)

    return orientation_dict

