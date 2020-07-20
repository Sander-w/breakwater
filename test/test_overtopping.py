import context
import unittest
import warnings
import numpy as np

import breakwater.core.overtopping as overtopping


class TestOvertopping(unittest.TestCase):

    def test_rubble_mound(self):
        # EurOtop (2018) case study 7 --> run 3
        with warnings.catch_warnings(record=True) as w:
            Rc = overtopping.rubble_mound(
                Hm0=2.72, q=1.18, xi_m_min_1=1.93, alpha=np.arctan(1/3),
                beta=0, gamma_b=1, gamma_v=1, gam_star=1, armour_layer='Rock',
                layers=2, permeability='permeable', safety=1, Gc=6.4,
                Dn50=1.28)

        self.assertAlmostEqual(Rc, 3.65, 2)
        self.assertEqual(len(w), 1)
        self.assertEqual(w[0].category, UserWarning)

        # EurOtop (2018) case study 8
        Rc = overtopping.rubble_mound(
            Hm0=4.9, q=21, xi_m_min_1=3.89, alpha=np.arctan(1/1.5), beta=0,
            gamma_b=1, gamma_v=1, gam_star=1, armour_layer='Accropode II',
            safety=1)
        self.assertAlmostEqual(Rc, 5.61, 2)

    def test_vertical_deep(self):
        # EurOtop (2018) case study 6, run 2
        # test of vertical_deep --> eq7.1
        with warnings.catch_warnings(record=True) as w:
            Rc = overtopping.vertical(
                Hm0=2.53, q=0.059, h=10.28, d=10.28, L_m_min_1=44,
                s_m_min_1=0.057, safety=0, limit=False)
            self.assertAlmostEqual(Rc, 5.94, 2)
            Rc = overtopping.vertical(
                Hm0=2.53, q=0.059, h=10.28, d=10.28, L_m_min_1=44,
                s_m_min_1=0.057, safety=0, limit=True)
            self.assertAlmostEqual(Rc, 5.56, 2)

        # test warnings
        self.assertEqual(len(w), 2)
        self.assertEqual(w[0].category, UserWarning)

    def test_composite_normal(self):
        # EurOtop (2018) case study 10, run 1
        # test of composite_normal --> eq7.14
        Rc = overtopping.vertical(
            Hm0=2.47, q=4.2, h=5.04, d=2.79, L_m_min_1=44,
            s_m_min_1=0.0554, safety=0)
        self.assertAlmostEqual(Rc, 5.60, 2)

    def test_vertical_no_breaking(self):
        # EurOtop (2018) case study 12, run 1
        # test of vertical_no_breaking --> eq7.5 (vertical wall)
        Rc = overtopping.vertical(
            Hm0=2.22, q=0.429, h=4.97, d=4.97, L_m_min_1=45,
            s_m_min_1=0.0497, safety=0)
        self.assertAlmostEqual(Rc, 5.67, 2)

    def test_vertical_normal(self):
        # EurOtop (2018) case study 12, run 4
        # test of vertical_normal --> eq7.8
        Rc = overtopping.vertical(
            Hm0=1.56, q=0.35, h=3.43, d=3.43, L_m_min_1=55.6,
            s_m_min_1=0.0278, safety=0)
        self.assertAlmostEqual(Rc, 7.21, 2)

    def test_vertical_low(self):
        # case study 12, run 4 adapted so that freeboard is low
        # test of vertical_low --> eq7.7
        Rc = overtopping.vertical(
            Hm0=1.56, q=16.2, h=3.43, d=3.43, L_m_min_1=55.6,
            s_m_min_1=0.0278, safety=0)
        self.assertAlmostEqual(Rc, 1.999, 3)

    def test_composite_low(self):
        # case study 10, run 1 adapted so that freeboard is low
        # test of composite_low --> eq7.15
        Rc = overtopping.vertical(
            Hm0=2.47, q=20.4, h=5.04, d=2.79, L_m_min_1=44,
            s_m_min_1=0.0554, safety=0)
        self.assertAlmostEqual(Rc, 3.298, 3)

    def test_composite_no_breaking(self):
        # case study 10, adapted so that there is no breaking
        # test of vertical_no_breaking --> eq7.5 (composite vertical)
        Rc = overtopping.vertical(
            Hm0=1.5, q=0.13, h=5.80, d=3.45, L_m_min_1=20,
            s_m_min_1=0.0532, safety=0)
        self.assertAlmostEqual(Rc, 4.156, 3)


if __name__ == '__main__':
    unittest.main()
