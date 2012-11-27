from sequenceminer.miner import subset_to_support,Element,Event
import unittest

class TestSubsetToSupport(unittest.TestCase):
    '''Test temporal joins of sequence ID lists.'''

    def test_subset_items(self):
        '''Test subsetting a list of sequences to individual items that meet or exceed a given support threshold.'''
        
        from sequenceminer.keydefaultdict import _KeyDefaultDict

        # ensure that multiple occurrences in the same sequence do not inflate support
        sequences = [
            (0,0,('A',)),
            (1,0,('B',)),
            (2,0,('A',)),
            (3,0,('C',)),
            (3,1,('C',)),
            ]
        
        # parse input sequences into individual item Elements
        elements = _KeyDefaultDict(Element) 

        for sid,eid,itemset in sequences:
            for item in itemset:
                elements[tuple(item)] |= Element(tuple(item),Event(sid=sid,eid=eid))
        
        # identify frequent single elements
        elements = subset_to_support(elements,2)

        self.assertEqual(elements.values(),[
            Element(('A',),Event(sid=0,eid=0),Event(sid=2,eid=0)),
            ])

if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(TestTemporalJoin)
    unittest.TextTestRunner(verbosity=2).run(suite)
