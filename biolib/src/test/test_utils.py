'''
Mock objects to have outside projects independence
'''
def _complement(seq):
    '''It complements the given sequence'''
    compdict = { 'a':'t', 'c':'g', 'g':'c', 't':'a', 'u':'t',
                 'm':'k', 'r':'y', 'w':'w', 's':'s', 'y':'r',
                 'k':'m', 'v':'b', 'h':'d', 'd':'h', 'b':'v',
                 'x':'x', 'n':'n',
                 'A':'T', 'C':'G', 'G':'C', 'T':'A', 'U':'T',
                 'M':'K', 'R':'Y', 'W':'W', 'S':'S', 'Y':'R',
                 'K':'M', 'V':'B', 'H':'D', 'D':'H', 'B':'V',
                 'X':'X', 'N':'N', '*':'*', '-':'-', ' ':' '
                 }
    return ''.join([compdict[base] for base in seq])
class Seq(object):
    '''This inmutable seq class has only test porpouses'''
    # pylint: disable-msg=R0903
    #We know that there are too few public methods, we don't need more
    def __init__(self, seq):
        '''We can set the sequence using an str.'''
        self.seq = seq
    def complement(self):
        '''It isn't mutable, so we return a new sequence'''
        compseq = _complement(self.seq)
        return self.__class__(compseq)
    def __len__(self):
        '''It returns the length'''
        return len(self.seq)
    def __str__(self):
        '''The string representation'''
        return self.seq
    def __getitem__(self, index):
        '''The slicing'''
        return self.seq[index]

class SeqWithQuality(object):
    '''A simple seq with Quality class to do some tests.'''
    # pylint: disable-msg=R0903
    #We know that there are too few public methods, we don't need more
    def __init__(self, seq, qual):
        '''The init requieres a sequence and a quality.'''
        self.seq = seq
        if not isinstance(qual, list):
            self.qual = [qual]
        else:
            self.qual = qual
    def __len__(self):
        '''The length'''
        return len(self.seq)
    def __getitem__(self, index):
        '''The slicing'''
        return self.__class__(seq=self.seq[index], qual=self.qual[index])
    def complement(self):
        '''It returns a complemented sequence'''
        return self.__class__(seq=_complement(self.seq), qual=self.qual)
    def __add__(self, other):
        '''It adds another seqwithquality.'''
        seq = self.seq + other.seq
        qual = self.qual[:]
        qual.extend(other.qual)
        return self.__class__(seq=seq, qual=qual)
class SeqRecord(object):
    '''A simple seq with name do some tests.'''
    # pylint: disable-msg=R0903
    #We know that there are too few public methods, we don't need more

    def __init__(self, seq, name=''):
        '''The init requires a sequence.'''
        self.seq = seq
        self.name = name
        self.letter_annotations = {}
    def __len__(self):
        '''The length'''
        return len(self.seq)
    def __getitem__(self, index):
        '''It returns a new sliced SeqRecord.'''
        return self.__class__(seq=self.seq[index])
    def complement(self):
        '''It returns a complemented SeqRecord'''
        return self.__class__(seq=_complement(self.seq))
    def __add__(self, other):
        '''It returns a new SeqRecord with seqs added.'''
        seq = self.seq + other.seq
        return self.__class__(seq=seq)
    
class Seqmut(object):
    '''This mutable seq class has only test porpouses'''
    # pylint: disable-msg=R0903
    #We know that there are too few public methods, we don't need more
    def __init__(self, seq):
        '''We can set the seq with an str'''
        self.seq = seq
    def complement(self):
        '''It's mutable, so we change the sequence'''
        self.seq = _complement(self.seq)
    def __len__(self):
        '''It returns the length'''
        return len(self.seq)
    def __str__(self):
        '''The string representation'''
        return self.seq
    def __getitem__(self, index):
        '''The slicing'''
        return self.seq[index]
