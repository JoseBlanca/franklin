'''
Created on 2009 mar 27

@author: peio
'''

class SeqWithQuality(object):
    '''
    This class is a container of a sequence with quality atribute
    '''

    def __init__(self, name, seq, qual= None):
        '''
        To initialice this class we need name and seq arguments.
        Arguments:
            name : Name of the secuence
            seq  : Sequence in a string, it accepts Biopython seq, 
                that is used to complement the seq. 
            qual : Quality of the sequence, in a list object
        '''
        
        self.name = name
        self._seq  = seq
        self._qual = None
        self.qual  = self._set_qual(qual)
     
    def __add__(self, seq2):
        '''It return a new object with both seqs joined '''
        return self.__class__(name = self.name, \
                              seq  = self._seq + seq2.seq, \
                              qual = self._qual + seq2.qual)
        
    def __getitem__(self, index):
        ''' It returns another object but only with the sclice chosed'''
        if isinstance(index, int):
            qual = [self.qual[index]]
        elif isinstance(index,slice):
            qual = self.qual[index]
        name = '%s_%s' % (self.name, str(index))
        return self.__class__(name = name, seq  = self._seq[index], \
                              qual = qual )
    def __len__(self):       
        ''' It returns the length of the sequence'''
        return len(self._seq)
     
    def __repr__(self):
        ''' It prints the seq in a readable'''
        sprint = 'Name : ' + self.name.__repr__() + '\n'
        sprint += 'Seq  : ' + self._seq.__repr__()  + '\n'
        sprint += 'Quality : '
        for qual in self._qual:
            sprint += '%d ' % int(qual)
        sprint += '\n'
        return sprint 
    
    
            
    def complement(self):
        ''' it returns a new object with the complementary strand of the seq '''
        
        return self.__class__(name = self.name + '_complemented', \
                              seq  = self._seq.complement(), \
                              qual = self._qual) 
    def _get_qual(self):
        ''' It returns the quality of the sequence'''
        return self._qual
    
    def _set_qual(self, qual):
        '''  It sets the quality of the sequence and it checks if 
        it is of the same lenght'''
        if qual is not None:
            if (len(qual) != len(self._seq)):
                raise ValueError('Quality values must have the same \
                lenght of sequence')
            else:
                self._qual = qual
    qual = property(_get_qual, _set_qual)   
    
    def _get_seq(self):
        ''' It returns the seq of the sequence '''
        return self._seq
    seq = property(_get_seq)
    
        