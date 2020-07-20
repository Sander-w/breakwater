import context
import unittest
import warnings
import numpy as np

from breakwater.core.goda import goda_wave_heights
from breakwater.conditions import LimitState
from breakwater.material import RockGrading, Xbloc
from breakwater.caisson import Caisson


class TestCaisson(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        warnings.filterwarnings('ignore')

        H13, Hmax = goda_wave_heights(
            h=15, d=11, Ho=5.7, T=7.28, slope_foreshore=(1,100))
        H13_SLS, Hmax_SLS = goda_wave_heights(
            h=13, d=8, Ho=4, T=7.18, slope_foreshore=(1,100))

        ULS = LimitState(
            deep_water=False, h=15, Nod=0.5, Sd=2, Ho=5, q=10, Hm0=5.22,
            H13=H13, Hmax=Hmax, T_m_min_1=8.66, Tp=7.87, Tm=6.3, T13=7.28,
            label='ULS')
        SLS = LimitState(
            deep_water=False, h=13, Nod=0.5, Sd=2, Ho=4, q=20, Hm0=4.22,
            H13=H13_SLS, Hmax=Hmax_SLS, T_m_min_1=8.66, Tp=7.87, Tm=6.3,
            T13=7.18, label='SLS')

        NEN = RockGrading()

        cls.bw = Caisson(
            Pc=0.6, rho_c=2400, rho_fill=1600, rho_w=1000, Bm=8, hb=4, mu=0.5,
            layers=2, BermMaterial=NEN, LimitState=[SLS, ULS], safety=1,
            slope_foreshore=(1,100), SF_sliding=1.2, SF_turning=1.2,
            beta=0, slope_foundation=(2,3))

    def test_structure_caisson(self):
        caisson_structure = self.bw.structure['caisson']
        computed = {
            'hb': 4, 'h_acc': 11, 'Pc': 0.6, 'd': 8.333,
            'Rc': 14.222, 'state_overtop': 1,
            'B': 18.888, 'state_goda': 1, 'Bm': 8}
        for i, (param, val) in enumerate(caisson_structure.items()):
            self.assertAlmostEqual(val, list(computed.values())[i], 3)

    def test_structure_armour(self):
        armour = self.bw.structure['armour']
        computed = [1.334, 'HMA_6000/10000', 1.445, 1, 2]
        for i, (param, val) in enumerate(armour.items()):
            self.assertAlmostEqual(val, computed[i], 3)

    def test_structure_underlayer(self):
        underlayer = self.bw.structure['foundation']
        computed = [[0.586, 0.671],
                    ['HMA_300/1000', 'HMA_1000/3000'],
                    [0.615, 0.895],
                    'see armour']
        for i, (param, vals) in enumerate(underlayer.items()):
            for j, val in enumerate(vals):
                self.assertAlmostEqual(val, computed[i][j], 3)

    def test_effective_width_foundation(self):
        B_eff = self.bw._effective_width_foundation(
            rho_w=1000, hb=4, B=18.895, m=991289)
        self.assertAlmostEqual(B_eff, 22.309, 3)


if __name__ == '__main__':
    unittest.main()
