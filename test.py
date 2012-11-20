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


    def test_temporal_join(self):
        '''Test temporal joins of sequence ID lists.'''

        from spade import temporal_join

        # simple join of disjoint sequences
        id_list_i = [{'item':('A',),'sid':1,'eid':0}]
        id_list_j = [{'item':('B',),'sid':1,'eid':1}]
        join_result = temporal_join(id_list_i,id_list_j)

        self.assertEqual(join_result,
            { ('A','B',): [{'item':('A','B',),'sid':1,'eid':1}] })

        # join with one overlapping element
        id_list_i = [{'item':('A','B',),'sid':1,'eid':1}]
        id_list_j = [{'item':('B','C',),'sid':1,'eid':2}]
        join_result = temporal_join(id_list_i,id_list_j)

        self.assertEqual(join_result,
            { ('A','B','C',): [{'item':('A','B','C',),'sid':1,'eid':2}] })

        # join with two overlapping elements
        id_list_i = [{'item':('A','B','C',),'sid':1,'eid':2}]
        id_list_j = [{'item':('B','C','D',),'sid':1,'eid':3}]
        join_result = temporal_join(id_list_i,id_list_j)

        self.assertEqual(join_result,
            { ('A','B','C','D',): [{'item':('A','B','C','D',),'sid':1,'eid':3}] })


    def test_subset_to_support(self):
        '''Test subsetting a list of sequences to those that meet or exceed a given support threshold.'''
        
        from spade import subset_to_support,IdList

        # ensure that multiple occurrences in the same sequence do not inflate support
        sequences = [
            ('A',),
            ('B',),
            ('A',),
            ('C','C',)
            ]
        
        id_list = IdList(sequences)

        id_list_subset = subset_to_support(id_list,2)

        self.assertEqual(id_list_subset,
            { ('A',): [{'item':('A',),'sid':0,'eid':0},{'item':('A',),'sid':2,'eid':0}] })



if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(TestSPADE)
    unittest.TextTestRunner(verbosity=2).run(suite)
