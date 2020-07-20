import context
import unittest
import warnings
import numpy as np

from breakwater.design import Configurations
from breakwater.conditions import LimitState
from breakwater.material import RockGrading, Xbloc


class TestCombinations(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # define the config file, with the input
        # of this input all possible combinations are generated
        cls.config = {
            'A': None,
            'B': np.linspace(1, 5, 5),
            'C': np.linspace(11, 15, 5),
            'D': 18,
            'E': np.linspace(21, 25, 5)}

        # define vkwargs to set varying and fixed parameters
        vkwargs = {
            'A': {'Constant': True},
            'B': {'Constant': False},
            'C': {'Constant': False},
            'D': {'Constant': True},
            'E': {'Constant': False}}

        # generate combinations
        cls.varying, cls.max_combinations = Configurations._get_combinations(
            vkwargs=vkwargs, config=cls.config)

    def test_max_combinations(self):
        # check maximum number of combinations
        self.assertEqual(self.max_combinations, 125)

    def test_check_combinations(self):
        # check if A and D are not in the varying dict
        self.assertNotIn('A', self.varying)
        self.assertNotIn('D', self.varying)

        # define index for B, C and E
        i_B, i_C, i_E = 0, 0, 0

        # check combinations
        for i in range(self.max_combinations):
            # get unique combination
            combination = Configurations._get_concept_set(
                configs=self.varying, id=i+1)

            # check combination
            self.assertEqual(combination['B'], self.config['B'][i_B])
            self.assertEqual(combination['C'], self.config['C'][i_C])
            self.assertEqual(combination['E'], self.config['E'][i_E])

            # increment E index every combination
            i_E += 1

            # increment C index every 5 combinations
            if (i+1) % 5 == 0:
                i_C += 1

                # reset E index
                i_E = 0

            # increment B index every 25 combinations
            if (i+1) % 25 == 0:
                i_B += 1

                # reset C and E index
                i_C, i_E = 0, 0


class TestConfigurations(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        ULS = LimitState(
            h=15, Nod=0.5, Sd=2, Ho=5, q=10, Hm0=3.22, T_m_min_1=8.66,
            label='ULS', Hmax=6)
        SLS = LimitState(
            h=10, Nod=0.5, Sd=2, Ho=4, q=20, Hm0=2.52, T_m_min_1=7.66,
            label='SLS')

        with warnings.catch_warnings(record=True) as w:
            ULS.transform_periods(0.5)
            SLS.transform_periods(0.5)
            ULS.check_deep_water()
            SLS.check_deep_water()

        xbloc = Xbloc()
        NEN = RockGrading()

        cls.configs = Configurations(
            structure=['RRM', 'CRM', 'RC', 'CC'], LimitState=ULS,
            slope=((1,6), (3,4), 6), slope_foreshore=(1,100), rho_w=1000,
            B=(3,8,5), N=3000, P=0.4, Grading=NEN, Dn50_core=(0.1, 0.4, 5),
            safety=1, slope_toe=(2,3), ArmourUnit=xbloc, rho_c=2400,
            rho_fill=1600, Pc=(0.6, 0.8, 5), Bm=(4,8,5),  hb=(1.5,5,5),
            mu=0.5, SF_sliding=1.2, SF_turning=1.2, beta=0,
            slope_foundation=(2,3), BermMaterial=xbloc)

    def test_values(self):
        # test parameters and values in the df
        # define list with columns to exclude as these are not parameters
        exclude = ['concept', 'id', 'type', 'warnings']

        # drop columns in exclude and define new df
        df = self.configs.df.drop(columns=exclude)

        # defined values
        defined = {
            'B': [3, 4.25, 5.5, 6.75, 8, np.nan],
            'B_toe': [None, np.nan],
            'Bm': [np.nan, 4, 5, 6, 7, 8],
            'Dn50_core': [0.1, 0.175, 0.25, 0.325, 0.4, np.nan],
            'Pc': [np.nan, 0.6, 0.65, 0.7, 0.75, 0.8],
            'hb': [np.nan, 1.5, 2.375, 3.25, 4.125, 5],
            'slope': [9.46, 15.82, 21.8, 27.32, 32.35, 36.87, np.nan],
            'slope_foundation': [np.nan, 33.69],
            'slope_toe': [33.69, np.nan]}

        # iterate over the columns of the df of the parameters
        for parameter in df.columns:
            # check if parameter is slope
            if 'slope' in parameter:
                # check if not all nan
                if not df[parameter].isnull().all():
                    # convert tuples to floats
                    df[parameter] = df[parameter].apply(
                        lambda x: convert(x))

            # get unique and correct values of the parameter
            unique = df[parameter].unique()
            correct = defined[parameter]

            # iterate over the unique values to check them
            for value, correct_value in zip(unique, correct):
                if value is None:
                    self.assertIsNone(correct_value)
                elif np.isnan(value):
                    np.testing.assert_equal(value, correct_value)
                else:
                    self.assertAlmostEqual(value, correct_value)

    def test_concepts(self):
        # set correct number of concepts and variants
        correct_num_variants = 1005
        correct_num_concepts = 550

        # set number of concepts and variants
        num_variants = 0
        num_concepts = 0

        # compute total number of generated concepts + variants
        for _, row in self.configs.df.iterrows():
            # get variantIDs of the concept
            variantIDs = row.concept.variantIDs

            # add to number of variants and concepts
            num_concepts += 1
            num_variants += len(variantIDs)

        # check numbers
        self.assertEqual(num_concepts, correct_num_concepts)
        self.assertEqual(num_variants, correct_num_variants)

    def test_types(self):
        # correct types
        specified = ['RRM', 'CRM', 'RC', 'CC']

        # get types
        types = list(self.configs.df.type.unique())

        self.assertListEqual(types, specified)


def convert(slope):
    if (isinstance(slope, tuple) or isinstance(slope, list)
            or isinstance(slope, np.void)):
        return np.round(np.arctan(slope[0]/slope[1])*180/np.pi, 2)
    else:
        return np.nan


if __name__ == '__main__':
    unittest.main()
