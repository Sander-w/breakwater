from pykml import parser
from shapely.geometry import LinearRing

def coords_from_kml(kml_file, root_location):

    """

    Parameters
    ----------
    kml_file: str
        string with the path of the kml_file
    root_location: loc
        location within the kml file with the coordinates
    Returns
    -------
    list
    """

    root = parser.fromstring(open(kml_file, 'rb').read())
    #coords = str(root.Document.Placemark.Polygon.outerBoundaryIs.LinearRing.coordinates).strip()
    coords = str(root.root_location).strip()

    coords = coords.split(' ')

    coord_tuples = []

    for c in coords:
        c_split = c.split(',')[:-1]
        coord_tuples.append((float(c_split[0]), float(c_split[1])))

    return coord_tuples

def create_shape(kml_file, shape = 'LinearRing'):

    """

    Parameters
    ----------
    coordinates: list
        list containing tuples with (lon, lat) of the LinearRing
    shape: str
        shape of the structure
    Returns
    -------

    """

    coords = coords_from_kml(kml_file)

    if shape == 'LinearRing':
        shape_structure = LinearRing(coords)

    return shape_structure




