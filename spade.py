import collections

class IdList(collections.Mapping):
    '''Dictionary implementation for ID Lists as described in Zaki (2001) that
    helps with bookkeeping and provides descriptive semantics.'''

    def __init__(self,sequences):
        '''Given a list of sequences, create a dictionary of (sid,eid) pairs
        for each atom in the sequence.
        '''
            
        self.events = {}

        for seq_index in range(len(sequences)):
            seq = sequences[seq_index]

            for element_index in range(len(seq)):
                element = seq[element_index]

                for item in element:

                    event = {'item':tuple(item),'sid':seq_index,'eid':element_index}
                    if event['item'] in self.events:
                        self.events[event['item']].append(event)
                    else:
                        self.events[event['item']] = [event]

    def add(self,atoms):

        for atom in atoms:
            if atom['item'] in self.events:

                atom_events = self.events[atom['item']]
                if (atom['sid'],atom['eid']) not in [(event['sid'],event['eid']) for event in atom_events]:
                    self.events[atom['item']].append(atom)

            else:
                self.events[atom['item']] = [atom]


    def __getitem__(self,item):

        return self.events[item]

    def __iter__(self):

        return iter(self.events)

    def __len__(self):

        return len(self.events)

    def __str__(self):

        return str(self.events)

    def __repr__(self):

        return self.events.__repr__()

def subset_to_support(atom_id_list,support_threshold):
    '''Given an IdList, return an IdList containing only those atoms which
    meet the support threshold.
    '''

    subsetted = IdList({})
    
    for atom in atom_id_list:
        if len(atom_id_list[atom]) >= support_threshold:
            subsetted.add(atom_id_list[atom])
                    
    return subsetted

def count_frequent_two_seq(id_list,support_threshold):
    '''Given an IdList of atoms, return a dictionary of two-sequences as keys with
    the frequency of each two-sequence as the value.
    '''

    # Given an IdList, convert it to a horizontal IdList in order to
    # count the frequency of each two-sequences of atoms.

    horizontal_db = {} 

    for atom,events in id_list.iteritems():
                
        for event in events:

            sid = event['sid']
            eid = event['eid']

            if sid not in horizontal_db:
                 horizontal_db[sid] = []

            horizontal_db[sid].append((atom,eid))

    # create counts using horizontal_db

    counts = {}
    
    for sid in horizontal_db.keys():
        seq = horizontal_db[sid]
        
        for event_index_i in range(len(seq)):
            event_i = seq[event_index_i]
                
            for event_index_j in range(event_index_i+1,len(seq)):
                event_j = seq[event_index_j]
                        
                if seq[event_index_i][1] <= seq[event_index_j][1]:
                    two_seq = event_i[0]+event_j[0]
                else:
                    two_seq = event_j[0]+event_i[0]

                if two_seq in counts:
                    counts[two_seq] += 1
                else:
                    counts[two_seq] = 1

    # this is followed by temporal joins between atoms in pairs, so
    # include only unique combinations
    return {tuple(sorted(two_seq)) for two_seq,count in counts.iteritems() if count >= support_threshold}


def temporal_join(id_list_i,id_list_j):
    '''Given two IdLists, return a dictionary of new IdLists indexed by
    tuples of the new atoms formed by a temporal join of the original
    IdLists.
    '''

    new_id_list = IdList({})
    
    atom_i = id_list_i[0]['item']
    atom_j = id_list_j[0]['item']
    
    for event_index_i in range(len(id_list_i)):
        for event_index_j in range(len(id_list_j)):
                                    
            event_i = id_list_i[event_index_i]
            event_j = id_list_j[event_index_j]
    
            if event_i['sid'] == event_j['sid']:
                                        
                sid = event_i['sid']
                new_atom = {}
            
                # these two atoms occur in the same sequence
                # if they occur at different times (different eids), then
                # their combination atom has the later eid by Corollary 1 (Zaki 2001)
                if event_i['eid'] > event_j['eid']:
                    new_atom = {'item':atom_j+tuple(atom_i[-1]),'sid':sid,'eid':event_i['eid']}
                elif event_i['eid'] < event_j['eid']:
                    new_atom = {'item':atom_i+tuple(atom_j[-1]),'sid':sid,'eid':event_j['eid']}
                elif atom_i[-1] != atom_j[-1]:
                    new_atom = {'item':(atom_i+atom_j),'sid':sid,'eid':event_j['eid']}
        
                if len(new_atom) > 0:
                    new_id_list.add([new_atom])
                
    return new_id_list

def enumerate_frequent_seq(id_list,support_threshold):
    '''Recursively traverse the sequence lattice, generating frequent n+1-length
    sequences from n-length sequences provided in the id_list parameter.'''

    frequent_seq = IdList({})

    for atom_index_i in range(len(id_list)):
        atom_i = id_list.keys()[atom_index_i]

        frequent_seq_inner = IdList({})
            
        for atom_index_j in range(atom_index_i+1,len(id_list)):
            atom_j = id_list.keys()[atom_index_j]
            #print "atom_i: " + str(atom_i) + " atom_j: " + str(atom_j)

            R = temporal_join(id_list[atom_i],id_list[atom_j])

            for seq,events in R.items():
                n_distinct_sequences = len(set([event['sid'] for event in events]))
                if n_distinct_sequences >= support_threshold and seq not in frequent_seq_inner:
                    frequent_seq_inner.add(events)

        for item,idlist in frequent_seq_inner.events.items():
            frequent_seq.add(idlist)
        for item,idlist in enumerate_frequent_seq(frequent_seq_inner,support_threshold).events.items():
            frequent_seq.add(idlist)

    return frequent_seq


def spade(sequences,support_threshold):
    '''SPADE is performed in three distinct steps:
    1. Identify frequent single elements.
    2. Identify frequent two-element sequences.
    3. Identify all remaining sequences of three elements or more.
    '''

    # identify frequent single elements
    id_list_singles = subset_to_support(IdList(sequences),support_threshold)

    # identify frequent two-element sequences using a horizontal database
    hid_list = count_frequent_two_seq(id_list_singles,support_threshold)

    # generate ID lists for frequent two-element sequences discovered above
    id_list_pairs = IdList({})

    for two_seq in hid_list:

        R = temporal_join(id_list_singles[tuple(two_seq[0])],id_list_singles[tuple(two_seq[1])])

        for seq,events in R.items():
            n_distinct_sequences = len(set([event['sid'] for event in events]))
            if n_distinct_sequences >= support_threshold and seq not in id_list_pairs:
                id_list_pairs.add(events)

    # identify and generate ID lists for all remaining sequences
    freq = enumerate_frequent_seq(id_list_pairs,support_threshold)

    # collect all identified sequences of any length
    for item,idlist in id_list_pairs.events.items():
        freq.add(idlist)

    for item,idlist in id_list_singles.events.items():
        freq.add(idlist)

    # return frequent sequences
    return set(freq.events.keys())


if __name__ == "__main__":

    import pprint as pp

    sequences = [
        ('A','B',),
        ('B','A',),
        ('A','B',),
        ('B',)
        ]

    pp.pprint(spade(sequences,2))
