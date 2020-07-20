def settings(key):
    """ Lay-out settings """
    settings = {

        'lay-out': {
            'width': 600,
            'height': 750,
            'top_height': 30,
            'bottom_height': 80,
        },

        'font': {
            'large': ('Verdana', 12),
            'info': ('Verdana', 10)
        },

        'navigation': {
            'height': 3,
            'width': 20,
            'padx': 5,
            'pady': 5,
        },

        'interactive': {
            'width': 1600,
            'height': 900,
            'design_width': 1000,
            'design_height': 400,
            'info_width': 300,
            'dpi': 150
        },

        'label': {
            'padx': 5,
            'pady': 1,
        },

        'title': {
            'font': ('Verdana', 12, 'bold'),
            'pady': 5
        },

    }

    return settings[key]

def size(to_design, frame, num_concepts):
    """ Returns required frame size based on the structures to design

    Parameters
    ----------
    to_design : dict
        dictionary with the structures to design, key is the structure,
        value is a bool
    frame : {interactive, design}
        for which frame the geometry must be returned
    num_concepts : int
        number of concepts

    Returns
    float, float
        the required width and height of the specified frame
    """
    # width is fixed
    width = 1600

    # check number of structures
    if num_concepts == 1:
        # set frame height
        # check if rubble mound or vertical
        if to_design['RM']['RRM'] or to_design['RM']['CRM']:
            # RubbleMound
            height = 500
            design_height = 400
            design_width = 1000

        elif to_design['C']['RC'] or to_design['C']['CC']:
            # vertical
            height = 600
            design_height = 500
            design_width = 1000

    elif num_concepts == 2:
        # set frame heights
        height = 860
        if to_design['C']['RC'] and to_design['C']['CC']:
            width = 1400
            design_height = 380
            design_width = 800

        else:
            design_height = 380
            design_width = 1000

    if frame == 'interactive':
        return width, height

    elif frame == 'design':
        return design_width, design_height

    else:
        raise NotImplementedError(f'{frame} is not implemented')
