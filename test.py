from spade import spade
import unittest

class TestSPADE(unittest.TestCase):
    '''Functional tests for SPADE.'''

    def test_cardinality_eq_1(self):
        '''Test identification of frequent one-element sequences.'''

        sequences = [
            (0,('A',)),
            (1,('B',)),
            (2,('A',)),
            (3,('C',)),
            (4,('A',)),
            (5,('B',))
            ]

        self.assertEqual(spade(sequences,2),set([('A',),('B',)]))

    def test_cardinality_eq_2(self):
        '''Test identification of frequent two-element sequences.'''

        sequences = [
            (0,('A','B',)),
            (1,('B','A',)),
            (2,('A','B',)),
            (3,('B',))
            ]

        self.assertEqual(spade(sequences,2),set([('A',),('B',),('A','B',)]))

    def test_cardinality_gt_2(self):
        '''Test identification of frequent sequences with more than two elements.'''

        sequences = [
            (0,('A','B','C','D',)),
            (1,('A','B','C','D',)),
            ]

        self.assertEqual(spade(sequences,2),set([
            ('A',),('B',),('C',),('D',),
            ('A','B',),('B','D',),('B','C',),('A','C',),('A','D',),('C','D',),
            ('A','B','C',),('B','C','D',),('A','B','D',),('A','C','D',),
            ('A','B','C','D',)
            ]))


    def test_temporal_join(self):
        '''Test temporal joins of sequence ID lists.'''

        from spade import temporal_join,Element,Event

        # simple join of disjoint sequences
        element_i = Element(('A',),Event(sid=1,eid=0))
        element_j = Element(('B',),Event(sid=1,eid=1))
        join_result = temporal_join(element_i,element_j)

        self.assertEqual(join_result.values(),[
            Element(('A','B',),Event(sid=1,eid=1))
            ])

        # join with one overlapping element
        element_i = Element(('A','B',),Event(sid=1,eid=1))
        element_j = Element(('B','C',),Event(sid=1,eid=2))
        join_result = temporal_join(element_i,element_j)

        self.assertEqual(join_result.values(),[ 
            Element(('A','B','C',),Event(sid=1,eid=2))
            ])

        # join with two overlapping elements
        element_i = Element(('A','B','C',),Event(sid=1,eid=2))
        element_j = Element(('B','C','D',),Event(sid=1,eid=3))
        join_result = temporal_join(element_i,element_j)

        self.assertEqual(join_result.values(),[ 
            Element(('A','B','C','D',),Event(sid=1,eid=3))
            ])


    def test_subset_to_support(self):
        '''Test subsetting a list of sequences to those that meet or exceed a given support threshold.'''
        
        from spade import subset_to_support,Element,Event
        from keydefaultdict import _KeyDefaultDict

        # ensure that multiple occurrences in the same sequence do not inflate support
        sequences = [
            (0,('A',)),
            (1,('B',)),
            (2,('A',)),
            (3,('C','C',))
            ]
        
        # parse input sequences into individual item Elements
        elements = _KeyDefaultDict(Element) 

        for sid,sequence in sequences:
            for eid,item in enumerate(sequence):
                elements[tuple(item)] |= Element(tuple(item),Event(sid=sid,eid=eid))

        
        # identify frequent single elements
        elements = subset_to_support(elements,2)

        self.assertEqual(elements.values(),[
            Element(('A',),Event(sid=0,eid=0),Event(sid=2,eid=0)),
            ])


if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(TestSPADE)
    unittest.TextTestRunner(verbosity=2).run(suite)
