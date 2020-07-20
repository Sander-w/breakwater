import context
import unittest
import numpy as np

from breakwater.core.bishop import Bishop, SlipCircle


class TestBishop(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # example from Soil Mechanics (CTB2310) at TU Delft
        # circle = {'xy': (2.23, 8.20), 'r': 8.5}
        circle = SlipCircle(centre=(2.23, 8.20), r=8.5)
        cls.slip = Bishop(point2=(10.44,6), SlipCircle=circle)
        cls.slip.add_layer(18, 25, phi=5, name='soil')
        cls.slip.compute(num_slices=5)

        # get the slices
        cls.slices = cls.slip.circles[1].slices

    def test_angle(self):
        # slip angles from the example
        # difference in angles is due to rounding intermediate values
        example = [-8.08, 6.14, 20.71, 36.98, 59.5]

        # iterate over the computed values
        for slice, params in self.slices.items():
            self.assertAlmostEqual(
                params['alpha_s'], example[slice]*np.pi/180, 3)

    def test_heights(self):
        # height of slice from the example
        example = [0.75, 1.99, 2.68, 2.69, 1.25]

        # iterate over the computed values
        for slice, params in self.slices.items():
            h = params['h']['soil'][0] - params['h']['soil'][1]
            self.assertAlmostEqual(h, example[slice], 2)

    def test_load(self):
        # loads of each slice from the example
        example = [-1.89, 3.82, 17.04, 29.12, 19.37]

        # iterate over the slices
        for slice, coords in self.slices.items():
            # compute load and strength
            load = self.slip._load(coords['alpha_s'], coords['h'])

            self.assertAlmostEqual(load, example[slice], 2)

    def test_strenght(self):
        # strength of each slice from the example
        # difference in values is again due to intermediate rounding
        # addionally the values of the example for the last two slices
        # do not yield the same result as the given result, the computed
        # result is closer to the result from the Bishop class
        example = [26.77, 28.03, 30.23, 34.33, 46.26]

        # iterate over the slices
        for slice, coords in self.slices.items():
            # compute load and strength
            strength = self.slip._strength(
                alpha_s=coords['alpha_s'], heights=coords['h'], F=1,
                pressure=False)
            self.assertAlmostEqual(strength, example[slice], 2)

    def test_F(self):
        self.assertAlmostEqual(self.slip.circles[1].F, 2.542, 3)


class TestBishopTwoLayers(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # example from Soil Mechanics (CTB2310) at TU Delft
        # circle = {'xy': (2.23, 8.20), 'r': 8.5}
        circle = SlipCircle(centre=(2.23, 8.20), r=8.5)

        slip = Bishop(point2=(10.44,6), SlipCircle=circle)
        slip.add_layer(gamma=18, c=25, phi=5, name='layer1', ymin=-20, ymax=0)
        slip.add_layer(gamma=18, c=25, phi=5, name='layer2', ymin=0, ymax=6)
        slip.compute(num_slices=5)

        # set as attribute of the class
        cls.slip = slip

        # get the slices
        cls.slices = cls.slip.circles[1].slices

    def test_angle(self):
        # slip angles from the example
        # difference in angles is due to rounding intermediate values
        example = [-8.08, 6.14, 20.71, 36.98, 59.5]

        # iterate over the computed values
        for slice, params in self.slices.items():
            self.assertAlmostEqual(
                params['alpha_s'], example[slice]*np.pi/180, 3)

    def test_heights(self):
        # height of slice from the example
        example = [0.75, 1.99, 2.68, 2.69, 1.25]

        # iterate over the computed values
        for slice, params in self.slices.items():
            # compute the height of the layer
            h = 0
            for layer, height in params['h'].items():
                h += height[0] - height[1]

            self.assertAlmostEqual(h, example[slice], 2)

    def test_load(self):
        # loads of each slice from the example
        example = [-1.89, 3.82, 17.04, 29.12, 19.37]

        # iterate over the slices
        for slice, coords in self.slices.items():
            # compute load and strength
            load = self.slip._load(coords['alpha_s'], coords['h'])
            self.assertAlmostEqual(load, example[slice], 2)

    def test_strenght(self):
        # strength of each slice from the example
        # difference in values is again due to intermediate rounding
        # addionally the values of the example for the last two slices
        # do not yield the same result as the given result, the computed
        # result is closer to the result from the Bishop class
        example = [26.77, 28.03, 30.23, 34.33, 46.26]

        # iterate over the slices
        for slice, coords in self.slices.items():
            # compute load and strength
            strength = self.slip._strength(
                alpha_s=coords['alpha_s'], heights=coords['h'], F=1,
                pressure=False)
            self.assertAlmostEqual(strength, example[slice], 2)

    def test_F(self):
        self.assertAlmostEqual(self.slip.circles[1].F, 2.542, 3)


if __name__ == '__main__':
    unittest.main()
