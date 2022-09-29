from breakwater.core.battjes import BattjesGroenendijk

class BattjesGroenendijk_3D:

    """
    Apply the battjesgroenendijk to the entire wave rose assuming the water depth and slope of the fore shore is
    the same everywhere.

    Parameters
    ----------
    h: float
        The water depth
    wave_conditions: dict
        dictionary with the orientation of the waves as the key and the value is a dictionary with possible keys:
        Hs, Hm0, Tp, Tm, T_m_min.
    waveheight: str
        which wave height from the wave_conditions should be used. Default is Hm0
    slope_foreshore: tuple
        slope of the foreshore

    """

    def __init__(self, h, wave_conditions, slope_foreshore, waveheight = 'Hm0'):

        self.battjes_dict = {}

        for orientation, wavespecs in wave_conditions.items():
            H = wavespecs[waveheight]
            battjes = BattjesGroenendijk(Hm0= H, h= h, slope_foreshore= slope_foreshore)
            self.battjes_dict[orientation] = battjes

    def get_Hp(self, P):

        Hp = {}

        for orientation, battjes in self.battjes_dict.items():
            Hp[orientation] = battjes.get_Hp(P= P)

        return  Hp
