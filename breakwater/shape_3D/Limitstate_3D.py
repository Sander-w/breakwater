from breakwater.conditions import LimitState

class LimitState_3D():

    def __init__(self, h, wave_conditions, H2_per, label, **kwargs):
        """

        Parameters
        ----------
        h
        wave_conditions
        H2_per
        label
        kwargs
        """

        self.Limit_states = {}
        for orientation, wavespecs in wave_conditions.items():
            H2 = {'H2_per': H2_per[orientation]}
            rest = {**wavespecs, **H2, **kwargs}
            LS = LimitState(h= h, label= f'{label}_{orientation}', **rest)
            self.Limit_states[orientation] = LS

    def transform_periods(self, scalar):

        for orientation, ls in self.Limit_states.items():
            ls.transform_periods(scalar= scalar)


