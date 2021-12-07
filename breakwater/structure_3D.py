from breakwater.shape_3D.shape import create_shape, coords_from_kml
import matplotlib.pyplot as plt
import math

def shape_3D(kml_path, shape = 'LinearRing'):
    """

    Parameters
    ----------
    kml_path: str
        path with the kml file
    shape: str
        shape of the structure. default is LinearRing

    Returns
    -------
    Shapely.shape

    """

    return create_shape(kml_file= kml_path, shape= shape)

def governing_wave_direction(kml_path, wave_directions, shape= 'LinearRing'):

    """

    Parameters
    ----------
    wave_directions: dict
        per wave direction the significant wave height and wave period is provided

    Returns
    -------

    """

    structure_coords = coords_from_kml(kml_file= kml_path)

    for i in range(len(structure_coords)-1):
        long2, long1 = structure_coords[i+1][0], structure_coords[i][0]
        lat2, lat1 = structure_coords[i+1][1], structure_coords[i][1]
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
            orientation = bearing - 90

        for key, items in wave_directions.items():
            angle = key - orientation
            if angle < 0:
                continue
            else:
                pass #work further here

        plt.figure()
        plt.plot([long1, long2], [lat1, lat2])
        plt.title(f'{orientation}')
