import context
import unittest
import warnings

from breakwater.utils.wave import dispersion
from breakwater.core.scour import scour_protection


class TestScourProtection(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.L = dispersion(T=8, h=15)

    def test_L(self):
        self.assertAlmostEqual(self.L, 81.7897, 4)

    def test_vertical(self):
        w = scour_protection(self.L, slope=None)
        self.assertAlmostEqual(w, 20.4474, 4)

    def test_bounds(self):
        w_lower = scour_protection(self.L, slope=(1,1.75))
        w_upper = scour_protection(self.L, slope=(1,1.2))
        self.assertAlmostEqual(w_lower, 8.1790, 4)
        self.assertAlmostEqual(w_upper, 12.2685, 4)

    def test_out_of_range(self):
        with warnings.catch_warnings(record=True) as w:
            width = scour_protection(self.L, slope=(1,2))

        self.assertAlmostEqual(width, 7.0637, 4)
        self.assertEqual(len(w), 1)
        self.assertEqual(w[0].category, UserWarning)


if __name__ == '__main__':
    unittest.main()
