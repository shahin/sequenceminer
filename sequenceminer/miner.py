import collections
from keydefaultdict import _KeyDefaultDict 

Event = collections.namedtuple('Event','sid eid')

class Element(object):
    '''An element of the set of all possible subsequences, and a description of
    where that element occurs in the input sequences.
    '''

    def __init__(self,seq,*events):

        self.seq = seq
        self.events = set()

        for event in events:
            self.events.add(event)

    def __ior__(self,other_element):
        '''Implements the assignment operator |= by returning an Element whose
        events attribute is a union of the events of both input Elements.
        '''

        self.events |= other_element.events
        return self

    def __repr__(self):

        return self.__dict__.__repr__()

    def __eq__(self,other):

        return (self.seq == other.seq and self.events == other.events)

      
def subset_to_support(elements,support_threshold):
    '''Given an IdList, return an IdList containing only those atoms which
    meet the support threshold.
    '''

    subsetted = _KeyDefaultDict(Element)

    for element_name,element in elements.items():
        support = len(set([event.sid for event in element.events]))
        if support >= support_threshold:
            subsetted[element_name] = element
                    
    return subsetted


def count_frequent_two_seq(elements,support_threshold):
    '''Given an IdList of atoms, return a dictionary of two-sequences as keys with
    the frequency of each two-sequence as the value.
    '''

    # Given an dictionary of Elements, convert it to a horizontal ID list in order to
    # count the frequency of each two-sequence of atoms.
    horizontal_db = {} 

    for element_name,element in elements.items():
                
        for event in element.events:

            if event.sid not in horizontal_db:
                 horizontal_db[event.sid] = []

            horizontal_db[event.sid].append((element_name,event.eid))

    # create counts using horizontal_db
    counts = collections.defaultdict(int)
    
    for sid,seq in horizontal_db.iteritems():
        
        for event_index_i,event_i in enumerate(seq):
            for event_index_j,event_j in enumerate(seq[event_index_i+1:]):
                        
                if event_i[1] <= event_j[1]:
                    two_seq = event_i[0]+event_j[0]
                else:
                    two_seq = event_j[0]+event_i[0]

                counts[two_seq] += 1

    # this is followed by temporal joins between atoms in pairs, so
    # include only unique combinations
    return {tuple(sorted(two_seq)) for two_seq,count in counts.iteritems() if count >= support_threshold}


def temporal_join(element_i,element_j):
    '''Given two elements, return a dictionary of new elements indexed by
    their corresponding item sequences.
    '''

    join_results = _KeyDefaultDict(Element)
    
    for event_index_i,event_i in enumerate(element_i.events):
        for event_index_j,event_j in enumerate(element_j.events):
    
            if event_i.sid == event_j.sid:
                                        
                sid = event_i.sid
                superseqs = tuple()
                superseqs_events = tuple()
            
                # these two atoms occur in the same sequence
                # if they occur at different times (different eids), then
                # their combination atom has the later eid by Corollary 1 (Zaki 2001)
                if event_i.eid > event_j.eid:
                    superseq = element_j.seq + tuple(element_i.seq[-1])
                    superseq_event = Event(sid=sid,eid=event_i.eid)
                    join_results[superseq] |= Element(superseq,superseq_event)

                elif event_i.eid < event_j.eid:
                    superseq = element_i.seq + tuple(element_j.seq[-1])
                    superseq_event = Event(sid=sid,eid=event_j.eid)
                    join_results[superseq] |= Element(superseq,superseq_event)

                elif element_i.seq[-1] != element_j.seq[-1]:

                    superseq_event = Event(sid=sid,eid=event_j.eid)

                    # for coincident atoms, join the last element of one atom to the other
                    # ensure that the itemset is sorted
                    superseq_i = element_i.seq[:-1] + tuple([
                        ''.join(sorted(set(element_i.seq[-1] + element_j.seq[-1])))
                        ])
                    join_results[superseq_i] |= Element(superseq_i,superseq_event)

                    superseq_j = element_j.seq[:-1] + tuple([
                        ''.join(sorted(set(element_i.seq[-1] + element_j.seq[-1])))
                        ])

                    # if both resulting atoms are identical, only add it once
                    if superseq_j != superseq_i:
                        join_results[superseq_j] |= Element(superseq_j,superseq_event)
                
    return join_results


def enumerate_frequent_seq(elements,support_threshold):
    '''Recursively traverse the sequence lattice, generating frequent n+1-length
    sequences from n-length sequences provided in the id_list parameter.'''

    frequent_elements = _KeyDefaultDict(Element)

    for element_index_i,seq_i in enumerate(elements.keys()):

        frequent_elements_inner = _KeyDefaultDict(Element)
            
        for element_index_j,seq_j in enumerate(elements.keys()[element_index_i+1:]):

            R = temporal_join(elements[seq_i],elements[seq_j])

            for seq,element in R.items():
                support = len(set([event.sid for event in element.events]))
                if support >= support_threshold:
                    frequent_elements_inner[seq] |= element


        for seq,element in frequent_elements_inner.items():
            frequent_elements[seq] |= element

        for seq,element in enumerate_frequent_seq(frequent_elements_inner,support_threshold).items():
            frequent_elements[seq] |= element

    return frequent_elements


def mine(sequences,support_threshold):
    '''SPADE (Zaki 2001) is performed in three distinct steps:
    1. Identify frequent single elements.
    2. Identify frequent two-element sequences.
    3. Identify all remaining sequences of three elements or more.
    '''

    # parse input sequences into individual item Elements
    elements = _KeyDefaultDict(Element) 

    for sid,eid,itemset in sequences:
        for item in itemset:
            elements[tuple(item)] |= Element(tuple(item),Event(sid=sid,eid=eid))

    # identify frequent single elements
    elements = subset_to_support(elements,support_threshold)

    # identify frequent two-element sequences using a horizontal database
    freq_elements_len_eq_2 = count_frequent_two_seq(elements,support_threshold)

    # generate ID lists for frequent two-element sequences discovered above
    elements_len_eq_2 = _KeyDefaultDict(Element)

    for two_seq in freq_elements_len_eq_2:

        R = temporal_join(elements[tuple(two_seq[0])],elements[tuple(two_seq[1])])

        for seq,element in R.items():
            support = len(set([event.sid for event in element.events]))
            if support >= support_threshold:
                elements_len_eq_2[seq] |= element

    # identify and generate ID lists for all remaining sequences
    freq = enumerate_frequent_seq(elements_len_eq_2,support_threshold)

    # collect all identified sequences of any length
    for seq,element in elements_len_eq_2.items():
        freq[seq] |= element

    for seq,element in elements.items():
        freq[seq] |= element

    # return frequent sequences
    return freq


def read_sequences(filename):
    '''Read sequences from a CSV.

    The CSV contains one line per sequence with columns defined as follows:
    - First column is a unique integer as sequence ID (sid)
    - Second column is a sequence-unique integer as event ID (eid)
    - Each remaining column contains an item as a character string with columns
      arranged in sequence order
    '''

    import csv

    sequences = []

    with open(filename) as f:
        seqreader = csv.reader(f,delimiter=',')
        for seqline in seqreader:
            sequences.append(
                    tuple([ seqline[0],seqline[1],tuple(seqline[2:]) ])
                    )

    return sequences


def main(argv):

    import pprint as pp
    import argparse

    parser = argparse.ArgumentParser(description=
            'Generate frequent subsequences meeting the minimum support threshold.'
            )

    parser.add_argument('--file',dest='input_sequence_file',
            help='A comma-delimited text file containing input sequences.',
            required=True)
    parser.add_argument('--support',dest='support_threshold',type=int,
            help='The minimum number of occurrences of a frequent sequence.',
            required=True)

    args = parser.parse_args(argv)
    sequences = read_sequences(args.input_sequence_file)

    frequent_sequences = mine(sequences,args.support_threshold)
    pp.pprint(frequent_sequences.keys())


if __name__ == "__main__":

    import sys

    main(sys.argv[1:])

