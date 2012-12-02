from sequenceminer.miner import temporal_join,Element,Event
import unittest

class TestTemporalJoin(unittest.TestCase):
    '''Test temporal joins of sequence ID lists.'''

    def test_join_disjoint(self):
        '''simple join of disjoint sequences'''
        element_i = Element(('A',),Event(sid=1,eid=0))
        element_j = Element(('B',),Event(sid=1,eid=1))
        join_result = temporal_join(element_i,element_j)

        self.assertEqual(join_result.values(),[
            Element(('A','B',),Event(sid=1,eid=1))
            ])

    def test_join_overlap_1(self):
        '''join with one overlapping element'''
        element_i = Element(('A','B',),Event(sid=1,eid=1))
        element_j = Element(('B','C',),Event(sid=1,eid=2))
        join_result = temporal_join(element_i,element_j)

        self.assertEqual(join_result.values(),[ 
            Element(('A','B','C',),Event(sid=1,eid=2))
            ])

    def test_join_overlap_2(self):
        '''join with two overlapping elements'''
        element_i = Element(('A','B','C',),Event(sid=1,eid=2))
        element_j = Element(('B','C','D',),Event(sid=1,eid=3))
        join_result = temporal_join(element_i,element_j)

        self.assertEqual(join_result.values(),[ 
            Element(('A','B','C','D',),Event(sid=1,eid=3))
            ])

    def test_join_coincident_items_1(self):
        '''join coincident items''' 
        element_i = Element(('A',),Event(sid=1,eid=1))
        element_j = Element(('B',),Event(sid=1,eid=1))
        join_result = temporal_join(element_i,element_j)

        self.assertEqual(join_result.values(),[
            Element(('AB',),Event(sid=1,eid=1))
            ])

    def test_join_coincident_items_2(self):
        '''join coincident items as parts of larger, overlapping atoms''' 
        element_i = Element(('P','A',),Event(sid=1,eid=1))
        element_j = Element(('P','B',),Event(sid=1,eid=1))
        join_result = temporal_join(element_i,element_j)

        self.assertEqual(join_result.values(),[
            Element(('P','AB',),Event(sid=1,eid=1))
            ])

    def test_join_coincident_items_3(self): 
        '''join coincident items as parts of larger, otherwise disjoint atoms'''
        element_i = Element(('Z','A',),Event(sid=1,eid=1))
        element_j = Element(('P','B',),Event(sid=1,eid=1))
        join_result = temporal_join(element_i,element_j)

        self.assertEqual(join_result.values(),[
            Element(('P','AB',),Event(sid=1,eid=1)),
            Element(('Z','AB',),Event(sid=1,eid=1)),
            ])

    def test_join_coincident_items_4(self):
        '''join coincident items as parts of itemsets'''

        element_i = Element(('A','B',),Event(sid=1,eid=1))
        element_j = Element(('BC',),Event(sid=1,eid=1))
        join_result = temporal_join(element_i,element_j)

        self.assertEqual(join_result.values(),[
            Element(('BC',),Event(sid=1,eid=1)),
            Element(('A','BC',),Event(sid=1,eid=1))
            ])

if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(TestTemporalJoin)
    unittest.TextTestRunner(verbosity=2).run(suite)
