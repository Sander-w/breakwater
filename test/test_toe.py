import context
import unittest
import numpy as np

import breakwater.core.toe as toe


class TestRock(unittest.TestCase):

    def test_toe_stability(self):
        Dn50 = toe.toe_stability(Hs=4, h=8, ht=5, Delta=1.6, Nod=0.5)
        self.assertAlmostEqual(Dn50, 0.7411, 4)

    def test_berm_stability(self):
        Dn50 = toe.toe_berm_stability(
            Hs=5.796, T=11.4, d=7.1, Bm=8, Delta=1.6)
        self.assertAlmostEqual(Dn50, 1.58, 2)


if __name__ == '__main__':
    unittest.main()
