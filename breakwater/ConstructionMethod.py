class Equipment:
    """
    General equipment class

    Parameters
    ----------
    grading_limit: str
        the maximum grading limit the equipment can handle
    design_type: dict
        Which layers can be designed with what grading at what cost, {layer: {grading: price}}
    operation_type: list of strings
        How is the material placed. e.g. ['bulk', 'individual']
    """

    def __init__(self, name, grading_limit, design_type, operation_type):

        self.name = name
        self.grading_limit = grading_limit
        self.design_type = design_type
        self.operation_type = operation_type

    def get_price(self, layer, grading_layer):
        return self.design_type[layer][grading_layer]

class Vessel(Equipment):
    """
    Vessel class (marine based equipment)

    Parameters
    ----------
    grading_limit: str
        the maximum grading limit the equipment can handle
    design_type: list of strings
        which layers can the equipment design. e.g. ['core', 'filter','underlayer']
    operation_type: list of strings
        How is the material placed. e.g. ['bulk', 'individual']
    draught: int
        draught of the vessel
    installation_depth: float
        At which water depth is the vessel used
    """

    def __init__(self,name, grading_limit, design_type, operation_type, draught, installation_depth):

        super().__init__(name, grading_limit, design_type, operation_type)
        self.draught = draught
        self.installation_depth = installation_depth

    def install(self, end, layer, grading_layer):
        install = False

        if layer in self.design_type.keys() and grading_layer in self.design_type[layer].keys():
            if (self.installation_depth - self.draught) >= end:
                install = True
        return install


class Crane(Equipment):

    """
    Crane class (Land based equipment)

    Parameters
    ----------
    grading_limit: str
        the maximum grading limit the equipment can handle
    design_type: list of strings
        which layers can the equipment design. e.g. ['core', 'filter','underlayer']
    operation_type: list of strings
        How is the material placed. e.g. ['bulk', 'individual']
    reach_y: int
        What is the deepest part the equipment can reach
    reach_x: int
        What is the furthest the equipment can rach
    installation_depth: float
        At which water depth is the crane used
    """

    def __init__(self, name, grading_limit, design_type, operation_type, reach_y, reach_x, installation_depth):

        super().__init__(name, grading_limit, design_type, operation_type)
        self.installation_depth = installation_depth
        self.reach_y = reach_y
        self.reach_x = reach_x

    def install(self, layer, grading_layer, start, max_h, xbot, xtop):
        install = False
        if layer in self.design_type.keys() and grading_layer in self.design_type[layer].keys():
            if (start > self.installation_depth) or \
                ((max_h - start) <= self.reach_y and abs(xbot - xtop) <= self.reach_x):
                install = True
        return install

class CraneVessel(Equipment):

    """
    Vessel class (Land based equipment)

    Parameters
    ----------
    grading_limit: str
        the maximum grading limit the equipment can handle
    design_type: list of strings
        which layers can the equipment design. e.g. ['core', 'filter','underlayer']
    operation_type: list of strings
        How is the material placed. e.g. ['bulk', 'individual']
    reach_x: int
        What is the furthest the equipment can rach
    draught: int
        draught of the vessel
    installation_depth: float
        At which water depth is the crane used
    """

    def __init__(self, name, grading_limit, design_type, operation_type, reach_x, draught, installation_depth):

        super().__init__(name, grading_limit, design_type, operation_type)
        self.reach_x = reach_x
        self.draught = draught
        self.installation_depth = installation_depth

    def install(self, layer, grading_layer, slope, start, x_end):
        V, H = slope
        install = False
        if layer in self.design_type.keys() and grading_layer in self.design_type[layer].keys():
            if self.reach_x - (self.draught * V / H + (start - self.installation_depth) * V / H) >= x_end:
                install = True
        return install


