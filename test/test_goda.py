import context
import unittest
import warnings
import numpy as np

from breakwater.core.goda import Goda, goda_wave_heights
from breakwater.utils.wave import shoaling_coefficient


class TestWaveHeights(unittest.TestCase):

    def test_shoaling_coefficient(self):
        # example from Goda (2000), page 78
        Ks = shoaling_coefficient(h=8, T=12, H0=4.5)
        self.assertAlmostEqual(Ks, 1.24, 2)

    def test_wave_heights(self):
        waves = goda_wave_heights(
            h=10.1, d=5.6, Ho=6.3, T=11.4, slope_foreshore=(1,100))
        self.assertAlmostEqual(waves[0], 5.796, 3)
        self.assertAlmostEqual(waves[1], 8.050, 3)


class TestGoda(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # example from Goda (2000), page 140
        with warnings.catch_warnings(record=True) as cls.w:
            cls.goda = Goda(
                Hs=5.796, Hmax=8.050, h=10.1, d=5.6, h_acc=7.1, hc=3.4, Bm=8,
                T=11.4, beta=15*np.pi/180, slope_foreshore=np.arctan(1/100),
                rho=1030, B=15)

    def test_warnings(self):
        self.assertEqual(len(self.w), 1)
        self.assertEqual(self.w[0].category, UserWarning)

    def test_impulsive(self):
        alpha_star = self.goda._impulsive_pressure(Bm=8, Hmax=8.050)
        self.assertAlmostEqual(alpha_star, 0.3223, 4)

    def test_pressures(self):
        self.assertAlmostEqual(self.goda.p1, 97644.84, 2)
        self.assertAlmostEqual(self.goda.p3, 87211.93, 2)
        self.assertAlmostEqual(self.goda.p4, 69674.10, 2)
        self.assertAlmostEqual(self.goda.pu, 65739.31, 2)

    def test_required_mass(self):
        mass = self.goda.required_mass(
            SF_sliding=1.16, SF_turning=2.56, mu=0.6)
        self.assertAlmostEqual(mass/1000, 345.4477, 4)

    def test_bearing_pressure(self):
        pe = self.goda.bearing_pressure(
            Pc=0.72142857, rho_c=2400, rho_fill=1600)
        self.assertAlmostEqual(pe/1000, 291.040, 3)

    def test_mass(self):
        mass = self.goda.mass(Pc=0.8, rho_c=2400, rho_fill=1600)
        self.assertEqual(mass/1000, 352.8)

    def test_centroid(self):
        y = self.goda._pressure_centroid()
        self.assertAlmostEqual(y, 5.155, 3)

    def test_eccentricity(self):
        mass = self.goda.mass(Pc=0.8, rho_c=2400, rho_fill=1600)
        e = self.goda.eccentricity(mass)
        self.assertAlmostEqual(e, 2.049, 3)

    def test_required_width(self):
        with warnings.catch_warnings(record=True) as w:
            B = self.goda.required_width(
                Pc=0.72142857, rho_c=2400, rho_f=1600, rho_w=1030, mu=0.6,
                t=0.5, SF_sliding=1.16, SF_turning=2.56)

            self.assertAlmostEqual(B, 15.200, 3)

        self.assertEqual(len(w), 1)
        self.assertEqual(w[0].category, UserWarning)


if __name__ == '__main__':
    unittest.main()
