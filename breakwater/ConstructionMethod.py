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
    reach: tuple
        What is the shallowest and deepest part the equiment can reach
    """

    def __init__(self, name, grading_limit, design_type, operation_type, reach):

        self.name = name
        self.grading_limit = grading_limit
        self.design_type = design_type
        self.operation_type = operation_type
        self.reach = reach

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
    reach: tuple, ints (min, max)
        What is the shallowest and deepest part the equiment can reach
    draught: float
        Maximum draught of the vessel
    installation_depth: float
        At which water depth is the vessel used
    """

    def __init__(self,name, grading_limit, design_type, operation_type, reach, draught, installation_depth):

        super().__init__(name, grading_limit, design_type, operation_type, reach)
        self.draught = draught
        self.installation_depth = installation_depth
        self.max_depth = installation_depth - draught
        self.shallow_rech = installation_depth - reach[0]
        self.deep_reach = installation_depth - reach[1]

class Crane(Equipment):

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
    reach: tuple
        What is the shallowest and deepest part the equiment can reach
    draught: float
        Maximum draught of the vessel
    installation_depth: float
        At which water depth is the vessel used
    """

    def __init__(self, grading_limit, design_type, operation_type, reach):

        super().__init__(grading_limit, design_type, operation_type, reach)



