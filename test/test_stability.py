import context
import unittest
import numpy as np

import breakwater.core.stability as stability


class TestRock(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.slope = np.arctan(1/3)

    def test_xi_cr(self):
        # deep water model constants, Rock Manual (2007) Box 5.13
        xi_cr_deep = stability.xi_critical(
            Cpl=6.2, Cs=1, P=0.4, alpha=self.slope)
        self.assertAlmostEqual(xi_cr_deep, 3.0, 1)

    def test_vandermeer_deep(self):
        # plunging conditons, Rock Manual (2007) Box 5.13
        Dn50 = stability.vandermeer_deep(
            Hs=5, Delta=1.6, P=0.4, Sd=5, N=2100, xi_m=1.85, alpha=self.slope,
            safety=0)
        self.assertAlmostEqual(Dn50, 1.26, 2)
        # plunging conditions + safety, Rock Manual (2007) Box 5.14
        Dn50 = stability.vandermeer_deep(
            Hs=5, Delta=1.6, P=0.4, Sd=5, N=2100, xi_m=1.85, alpha=self.slope,
            safety=1.75)
        self.assertAlmostEqual(Dn50, 1.42, 2)
        # test surging conditions
        Dn50 = stability.vandermeer_deep(
            Hs=4, Delta=1.6, P=0.4, Sd=5, N=2100, xi_m=4.46,
            alpha=np.arctan(1/1.5),safety=0)
        self.assertAlmostEqual(Dn50, 1.5519, 4)

    def test_vandermeer_shallow(self):
        # plunging conditions, Rock Manual (2007) Box 5.15 (with errata)
        Dn50 = stability.vandermeer_shallow(
            Hs=4, H2=4.95, Delta=1.6, P=0.4, Sd=2, N=2273, xi_s_min_1=2.39,
            alpha=self.slope, safety=0)
        self.assertAlmostEqual(Dn50, 1.27, 2)
        # surging conditons
        Dn50 = stability.vandermeer_shallow(
            Hs=3.5, H2=4.2, Delta=1.6, P=0.4, Sd=2, N=2273, xi_s_min_1=4.01,
            alpha=np.arctan(1/2), safety=0)
        self.assertAlmostEqual(Dn50, 1.3713, 4)

    def test_hudson(self):
        # Checked with values from Xbloc guidelines, table 1
        Dn = stability.hudson(H=5.01, Kd=16, Delta=1.33, alpha=np.arctan(3/4))
        self.assertAlmostEqual(Dn**3, 2.5, 1)


if __name__ == '__main__':
    unittest.main()
