import context
import unittest
import warnings
import numpy as np

from breakwater.conditions import LimitState
from breakwater.utils.exceptions import LimitStateWarning


class TestLimitState(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.ULS = LimitState(
            h=8.5, Nod=0.5, q=10, Hm0=3.22, H13=3.27, T_m_min_1=7.66, Tp=6.96,
            Tm=6.33, label='ULS')
        cls.SLS = LimitState(
            h=8.5, Sd=10, Hm0=2.82, T_m_min_1=8.58, label='SLS')

    def test_deep_water(self):
        self.assertFalse(self.ULS.deep_water)
        self.assertFalse(self.SLS.deep_water)

    def test_transform_periods(self):
        with warnings.catch_warnings(record=True) as w:
            self.SLS.transform_periods(scalar=0.5)

        self.assertEqual(len(w), 1)
        self.assertEqual(w[0].category, LimitStateWarning)

        parameter = ['Sd', 'Hm0', 'T_m_min_1', 'Tp', 'Tm', 'T13']
        computed = [10, 2.82, 8.58, 7.8, 6.24, 7.215]
        for i, (arg, val) in enumerate(self.SLS.conditions.items()):
            self.assertEqual(arg, parameter[i])
            self.assertEqual(val, computed[i])

    def test_get_Hs(self):
        with warnings.catch_warnings(record=True) as w:
            self.assertEqual(self.ULS.get_Hs('H13'), 3.27)
            self.assertEqual(self.ULS.get_Hs('Hm0'), 3.22)
            self.assertEqual(self.SLS.get_Hs('H13'), 2.82)

        self.assertEqual(len(w), 1)
        self.assertEqual(w[0].category, LimitStateWarning)

    def test_L(self):
        L0 = self.ULS.L(period='T_m_min_1', deep_water=True)
        L = self.ULS.L(period='T_m_min_1', deep_water=False)
        self.assertAlmostEqual(L0, 91.611, 3)
        self.assertAlmostEqual(L, 63.122, 3)

    def test_s(self):
        self.assertAlmostEqual(self.ULS.s(number='mean'), 0.05227, 5)
        self.assertAlmostEqual(self.ULS.s(number='spectral'), 0.03515, 5)
        self.assertAlmostEqual(self.ULS.s(H='Hm0', T='T_m_min_1'), 0.03515, 5)

    def test_surf_similarity(self):
        alpha = np.arctan(3/4)
        xi_mean = self.ULS.surf_similarity(alpha=alpha, number='mean')
        xi = self.ULS.surf_similarity(alpha=alpha, H='Hm0', T='T_m_min_1')
        self.assertAlmostEqual(xi_mean, 3.280, 3)
        self.assertAlmostEqual(xi, 4.000, 3)

    def test_Sd(self):
        self.ULS.Sd(G=1, nv=0.45)
        self.assertAlmostEqual(self.ULS.conditions['Sd'], 0.909, 3)

    def test_Nod(self):
        LS = LimitState(
            h=8.5, Sd=10, Hm0=2.82, T_m_min_1=8.58, label='LS')

        LS.Nod(G=1, nv=0.45)
        self.assertEqual(LS.conditions['Nod'], 5.5)


if __name__ == '__main__':
    unittest.main()
