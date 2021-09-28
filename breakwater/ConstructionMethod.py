class Equipment:
    """
    General equipment class

    Parameters
    ----------
    grading_limit: str
        the maximum grading limit the equipment can handle
    design_type: list of strings
        which layers can the equipment design. e.g. ['core', 'filter','underlayer']
    operation_type: list of strings
        How is the material placed. e.g. ['bulk', 'individual']
    reach_y: int
        From/to which depth can stones be placed. From or to dependent on equipment type (vessel or crane)
    """

    def __init__(self, name, grading_limit, design_type, operation_type, reach_y):

        self.name = name
        self.grading_limit = grading_limit
        self.design_type = design_type
        self.operation_type = operation_type
        self.reach_y = reach_y

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
    reach_y: int
        From how deep can the layer be placed
    installation_depth: float
        At which water depth is the vessel used
    """

    def __init__(self,name, grading_limit, design_type, operation_type, reach_y, draught, installation_depth):

        super().__init__(name, grading_limit, design_type, operation_type, reach_y)
        self.draught = draught
        self.installation_depth = installation_depth
        self.min_depth = installation_depth - reach_y

class Crane(Equipment):

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
    reach_y: int
        What is the deepest part the equiment can reach
    installation_depth: float
        At which water depth is the crane used
    """

    def __init__(self, grading_limit, design_type, operation_type, reach_y, installation_depth):

        super().__init__(grading_limit, design_type, operation_type, reach_y)
        self.installation_depth = installation_depth
        self.max_depth = installation_depth - reach_y



