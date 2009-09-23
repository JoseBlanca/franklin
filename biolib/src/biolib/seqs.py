'''
Created on 2009 mar 27

@author: peio
'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of biolib.
# biolib is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# biolib is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with biolib. If not, see <http://www.gnu.org/licenses/>.

class SeqWithQuality(object):
    'This class is a container of a sequence with quality atribute'
    #we don't need more public methods
    #pylint: disable-msg=R0903

    def __init__(self, seq=None, qual= None, name=None, length=None,
                description=None, annotations=None):
        '''To initialice this class we need name and seq arguments.

        Arguments:
            name : Name of the secuence
            seq  : Sequence in a string, it accepts Biopython seq,
                that is used to complement the seq.
            qual : Quality of the sequence, in a list object
            length : The sequence length, useful for empty seq objects
            annotations: A dict with the annotations that affect the sequence
                         as a whole, not just a fragment of the sequence
            description: an str with the sequence description
        '''
        self._length = length
        self.name = name
        self._seq = None
        self.seq  = seq
        self._qual = None
        self.qual  = qual
        self.description = description
        self.annotations = annotations

    def __add__(self, seq2):
        '''It returns a new object with both seq and qual joined '''
        if self._qual is  None or seq2.qual is None:
            qual = None
        else:
            qual = self._qual + seq2.qual

        return self.__class__(name = self.name, \
                              seq  = self._seq + seq2.seq, \
                              qual = qual)

    def __getitem__(self, index):
        ''' It returns another object but only with the slice chose'''
        qual = self.qual
        if qual is not None:
            if isinstance(index, int):
                qual = [qual[index]]
            elif isinstance(index,slice):
                qual = qual[index]
        else:
            qual = None
        return self.__class__(name=self.name, seq=self._seq[index], qual=qual,
                              description=self.description)

    def __len__(self):
        ''' It returns the length of the sequence'''
        if self._length is not None:
            return self._length
        elif self._seq is not None:
            return len(self._seq)
        else:
            return 0

    def __repr__(self):
        ''' It prints the seq in a readable'''
        sprint = 'Name : ' + self.name.__repr__() + '\n'
        sprint += 'Seq  : ' + self._seq.__repr__()  + '\n'
        sprint += 'Quality : '
        quals = self._qual
        if quals is not None:
            for qual in quals:
                sprint += '%d ' % int(qual)
            sprint += '\n'
        else:
            sprint += 'None'
        return sprint

    def __str__(self):
        'It returns just the str of the seq property'
        return str(self.seq)

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
        if self._qual is not None:
            raise AttributeError("Can't reset the qual attribute")
        length = len(self)
        if length and qual is not None and length != len(qual):
            if self.name is None:
                print_name = ''
            else:
                print_name = ' in ' + self.name
            raise ValueError('qual should have the same length as this ' +
                            self.__class__.__name__  + print_name)
        self._qual = qual
    qual = property(_get_qual, _set_qual)

    def _set_seq(self, seq):
        'It sets the seq property'
        if self._seq is not None:
            raise AttributeError("Can't reset the seq attribute")
        length = len(self)
        if length and seq is not None and length != len(seq):
            raise ValueError('seq should have the same length as this ' +
                            self.__class__.__name__)
        self._seq = seq

    def _get_seq(self):
        ''' It returns the seq of the sequence '''
        return self._seq
    seq = property(_get_seq, _set_seq)

    def copy(self, seq=None, qual=None, name=None):
        '''Given a seqrecord it returns a new seqrecord with seq or qual changed.

        This is necessary because our SeqWithQuality is inmutable
        '''
        if seq is None:
            seq = self.seq
        if qual is None:
            qual = self.qual
        if name is None:
            name = self.name
        return self.__class__(name=name, seq=seq, qual=qual,
                              annotations=self.annotations,
                              description=self.description)

class Seq(str):
    "It represents a sequence. It's basically an str with some extra methods"
    # too many public method is str's fault
    #pylint: disable-msg=R0904
    _super = str
    def __new__(cls, string, *args, **kargs):
        'It creates a new Seq instance'
        return str.__new__(cls, string)
    def __init__(self, seq):
        'It requires a sequence (str) to be initialized.'
        pass

    def complement(self):
        'It returns a complemented the sequence'
        compdict = { 'a':'t', 'c':'g', 'g':'c', 't':'a', 'u':'t',
                     'm':'k', 'r':'y', 'w':'w', 's':'s', 'y':'r',
                     'k':'m', 'v':'b', 'h':'d', 'd':'h', 'b':'v',
                     'x':'x', 'n':'n',
                     'A':'T', 'C':'G', 'G':'C', 'T':'A', 'U':'T',
                     'M':'K', 'R':'Y', 'W':'W', 'S':'S', 'Y':'R',
                     'K':'M', 'V':'B', 'H':'D', 'D':'H', 'B':'V',
                     'X':'X', 'N':'N', '*':'*', '-':'-', ' ':' '
                     }
        return self.__class__(''.join([compdict[base] for base in self]))

    def __add__(self, seq):
        'It returns an new added sequence'
        add_result = self._super.__add__(self, seq)
        return self.__class__(add_result)

