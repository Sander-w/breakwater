import context
import unittest

from breakwater.core.soil import Soil


class TestSoil(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.soil = Soil(c=20, phi=15, gamma=19, n=0.4)

    def test_brinch_hansen(self):
        # example from CTB2310 Soil Mechanics at the TU Delft
        # no inclination factors
        p = self.soil.brinch_hansen(p=0, t=0, B=10, L=10, q=19, sat=False)
        self.assertAlmostEqual(p, 462.51, 2)

        # with inclination factors
        p = self.soil.brinch_hansen(t=30, p=175, B=10, L=10, q=19, sat=False)
        self.assertAlmostEqual(p, 191.54, 2)

    def test_saturated_soil(self):
        self.soil.saturated_weight(rho_w=1000)
        self.assertEqual(self.soil.gamma_sat, 22.924)


if __name__ == '__main__':
    unittest.main()
