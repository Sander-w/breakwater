import context
import unittest
import numpy as np

from breakwater.core.battjes import BattjesGroenendijk


class TestBattjesGroenendijk(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.battjes = BattjesGroenendijk(
            Hm0=0.13, h=0.27, slope_foreshore=np.arctan(1/100))

    def test_gammainc_upper(self):
        self.assertAlmostEqual(
            self.battjes.gammainc_upper(0.01, 0.1), 1.80324, 5)
        self.assertAlmostEqual(
            self.battjes.gammainc_upper(2, 1), 0.73575888, 8)

    def test_gamma_lower(self):
        self.assertAlmostEqual(
            self.battjes.gammainc_lower(0.01, 0.1), 97.62934376, 8)
        self.assertAlmostEqual(
            self.battjes.gammainc_lower(1.5, 1), 0.37894469, 8)

    def test_solver(self):
        self.battjes.Htr_tilde = 0.1
        H1_H2_tilde = self.battjes._solver([7.003, 1.060])
        self.assertAlmostEqual(H1_H2_tilde[0], 0.000, 3)
        self.assertAlmostEqual(H1_H2_tilde[1], 0.000, 3)

    def test_get_Hp(self):
        self.assertAlmostEqual(self.battjes.get_Hp(0.01), 0.171, 3)
        self.assertAlmostEqual(self.battjes.get_Hp(0.001), 0.192, 3)
        self.assertRaises(ValueError, self.battjes.get_Hp, P=1)

    def test_get_Hn(self):
        self.assertAlmostEqual(self.battjes.get_Hn(3), 0.135, 3)
        self.assertRaises(ValueError, self.battjes.get_Hn, N=0)


if __name__ == '__main__':
    unittest.main()
