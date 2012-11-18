from spade import spade
import unittest

class TestSPADE(unittest.TestCase):
    '''Functional tests for SPADE.'''

    def test_cardinality_eq_1(self):
        '''Test identification of frequent one-element sequences.'''

        sequences = [
            ('A',),
            ('B',),
            ('A',),
            ('C',),
            ('A',),
            ('B',)
            ]

        self.assertEqual(spade(sequences,2),set([('A',),('B',)]))

    def test_cardinality_eq_2(self):
        '''Test identification of frequent two-element sequences.'''

        sequences = [
            ('A','B',),
            ('B','A',),
            ('A','B',),
            ('B',)
            ]

        self.assertEqual(spade(sequences,2),set([('A',),('B',),('A','B',)]))

    def test_cardinality_gt_2(self):
        '''Test identification of frequent sequences with more than two elements.'''

        sequences = [
            ('A','B','C','D',),
            ('A','B','C','D',),
            ]

        self.assertEqual(spade(sequences,2),set([
            ('A',),('B',),('C',),('D',),
            ('A','B',),('B','D',),('B','C',),('A','C',),('A','D',),('C','D',),
            ('A','B','C',),('B','C','D',),('A','B','D',),('A','C','D',),
            ('A','B','C','D',)
            ]))


if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(TestSPADE)
    unittest.TextTestRunner(verbosity=2).run(suite)
