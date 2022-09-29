import context
import unittest
import numpy as np

from breakwater.rubble_2D import RockRubbleMound, ConcreteRubbleMound
from breakwater.conditions import LimitState
from breakwater.material import RockGrading, Xbloc
from breakwater.utils.exceptions import no_warnings


class TestRockRubbleMound(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        no_warnings()

        ULS = LimitState(
            deep_water=False, h=10, Nod=0.5, label='ULS', Sd=2, q=10,
            Hm0=3.22, T_m_min_1=7.66, Tp=6.96, Tm=6.33)
        NEN = RockGrading()

        cls.bw = RockRubbleMound(
            slope=(3,4), slope_foreshore=(1,100), B=6.4, rho_w=1000,
            Dn50_core=0.3, LimitState=ULS, Grading=NEN, safety=1, id=1,
            N=3000)

    def test_structure_armour(self):
        armour = self.bw.structure['armour']
        computed = [1.399, 'HMA_6000/10000', 1.445, 0, 2]
        for i, (param, val) in enumerate(armour.items()):
            self.assertAlmostEqual(val, computed[i], 3)

    def test_structure_underlayer(self):
        underlayer = self.bw.structure['underlayer']
        computed = [[0.586, 0.671],
                    ['HMA_300/1000', 'HMA_1000/3000'],
                    [0.615, 0.895],
                    'see armour',
                    2]
        for i, (param, vals) in enumerate(underlayer.items()):
            if isinstance(vals, list):
                for j, val in enumerate(vals):
                    self.assertAlmostEqual(val, computed[i][j], 3)
            else:
                self.assertEqual(vals, computed[i])

    def test_structure_filter(self):
        filter = self.bw.structure['filter layer']
        computed = [[None, 0.306, 0.415],
                    [None, 'LMA_40/200', 'LMA_60/300'],
                    [None, 0.335, 0.388],
                    'see armour',
                    2]
        for i, (param, vals) in enumerate(filter.items()):
            if isinstance(vals, list):
                for j, val in enumerate(vals):
                    self.assertAlmostEqual(val, computed[i][j], 3)
            else:
                self.assertEqual(vals, computed[i])

    def test_structure_toe(self):
        toe = self.bw.structure['toe']
        computed = [0.790, 'HMA_1000/3000', 0.895, 0]
        for i, (param, val) in enumerate(toe.items()):
            self.assertAlmostEqual(val, computed[i], 3)

    def test_area(self):
        areas = self.bw.area('a')
        computed = [109.035, 52.574, 171.323, 9.613]
        for i, (layer, area) in enumerate(areas.items()):
            self.assertAlmostEqual(area, computed[i], 3)


class TestConcreteRubbleMound(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        no_warnings()

        ULS = LimitState(
            deep_water=False, h=10, Nod=0.5, q=10, Hm0=3.22, Tp=6.96, Tm=6.33,
            T_m_min_1=7.66, label='ULS')
        SLS = LimitState(
            deep_water=False, h=9, Nod=0.5, q=10, Hm0=2.52, Tp=6.96, Tm=6.33,
            T_m_min_1=7.66, label='ULS')
        NEN = RockGrading()
        xbloc = Xbloc()

        cls.bw = ConcreteRubbleMound(
            slope=(3,4), slope_foreshore=(1,100), B=6.4, rho_w=1000,
            Dn50_core=0.3, LimitState=[SLS, ULS], id=1, ArmourUnit=xbloc,
            Grading=NEN, safety=1)

    def test_structure_armour(self):
        armour = self.bw.structure['armour']
        computed = [0.949, 1, 1, 1, 1]
        for i, (param, val) in enumerate(armour.items()):
            self.assertAlmostEqual(val, computed[i], 3)

    def test_structure_underlayer(self):
        underlayer = self.bw.structure['underlayer']
        computed = [[0.392, 0.532],
                    ['LMA_60/300', 'HMA_300/1000'],
                    [0.388, 0.615],
                    'see armour',
                    2]
        for i, (param, vals) in enumerate(underlayer.items()):
            if isinstance(vals, list):
                for j, val in enumerate(vals):
                    self.assertAlmostEqual(val, computed[i][j], 3)
            else:
                self.assertEqual(vals, computed[i])

    def test_structure_filter(self):
        self.assertRaises(KeyError, lambda: self.bw.structure['filter layer'])

    def test_structure_toe(self):
        toe = self.bw.structure['toe']
        computed = [0.450, 'HMA_300/1000', 0.615, 1]
        for i, (param, val) in enumerate(toe.items()):
            self.assertAlmostEqual(val, computed[i], 3)


if __name__ == '__main__':
    unittest.main()
