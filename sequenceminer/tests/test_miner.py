from sequenceminer.miner import mine
import unittest

class TestMiner(unittest.TestCase):
    '''Functional tests for Miner.'''

    def test_cardinality_eq_1(self):
        '''Test identification of frequent one-element sequences.'''

        sequences = [
            (0,5,('A',)),
            (1,2,('B',)),
            (2,0,('A',)),
            (3,5,('C',)),
            (4,15,('A',)),
            (5,15,('B',))
            ]

        self.assertEqual(set(mine(sequences,2).keys()),set([
            ('A',),('B',)
            ]))

    def test_cardinality_eq_2(self):
        '''Test identification of frequent two-element sequences.'''

        sequences = [
            (0,0,('A',)),
            (0,1,('B',)),
            (1,0,('B',)),
            (1,1,('A',)),
            (2,0,('A',)),
            (2,1,('B',)),
            (3,0,('B',))
            ]

        self.assertEqual(set(mine(sequences,2).keys()),set([
            ('A',),('B',),('A','B',)
            ]))

    def test_cardinality_gt_2(self):
        '''Test identification of frequent sequences with more than two elements.'''

        sequences = [
            (0,0,('A',)),
            (0,1,('B',)),
            (0,2,('C',)),
            (0,3,('D',)),
            (1,0,('A',)),
            (1,1,('B',)),
            (1,2,('C',)),
            (1,3,('D',)),
            ]

        self.assertEqual(set(mine(sequences,2).keys()),set([
            ('A',),('B',),('C',),('D',),
            ('A','B',),('B','D',),('B','C',),('A','C',),('A','D',),('C','D',),
            ('A','B','C',),('B','C','D',),('A','B','D',),('A','C','D',),
            ('A','B','C','D',)
            ]))

if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(TestMiner)
    unittest.TextTestRunner(verbosity=2).run(suite)
