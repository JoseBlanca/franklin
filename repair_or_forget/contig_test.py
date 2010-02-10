'''It checks that the contig code is ok.'''

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

import unittest
from biolib.contig import Contig
from biolib.locatable_sequence import locate_sequence
from test.test_utils import Seq, SeqRecord, SeqRecord

# pylint: disable-msg=R0201
# disable to many public methods
# pylint: disable-msg=R0904
# disable to many statements
# pylint: disable-msg=R0915
# disable: Except doesn't do anything
# pylint: disable-msg=W0704

class ContigTests(unittest.TestCase):
    '''It checks the Contig class behaviour.'''
    def test_basic_behaviour(self):
        '''It tests the init and the basic slicing.'''
        #a basic assembly with all sequences starting at the same place
        #aligment
        #     0123456
        #seq1 ATCTGCA
        #seq2 ATGTGCA
        #seq3 ATGTGCT
        seq = 'ATCTGCA'
        alignment = Contig(sequences=[seq, 'ATGTGCA', 'ATGTGCT'])
        assert len(alignment) == 3
        assert alignment.ncols == 7
        #some slicing
        assert alignment[0, 0] == 'A'
        assert alignment[0] is seq
        assert alignment[1, 1:3] == 'TG'
        #now we get a column
        assert alignment[1:, 6] == 'AT'

        #the same with seqRecords
        seqrec1 = SeqRecord('ATCTGCA')
        seqrec2 = SeqRecord('ATGTGCA')
        alignment = Contig(sequences=[seqrec1, seqrec2, SeqRecord('ATGTGCT')])
        assert len(alignment) == 3
        assert alignment.ncols == 7
        #if we want the SeqRecod.seq we should ask for it
        assert alignment[0, 0].seq == 'A'
        #a seqrecord slice returns a seqrecord slice
        assert alignment[0] is seqrec1
        assert alignment[1, 1:3].seq == 'TG'
        #now we get a column
        assert alignment[:, 2].seq == 'CGG'

        #now a contig without consensus
        #ref   012345
        #seq1    0123
        #seq1    CTAT   Locatable
        #seq2  012345
        #seq2  tTCTag   Masked
        #seq3  43210
        #seq3  TAGAT    RevComp
        #seq4   43210
        #seq4   cGAgc   Locatable, Masked and RevComp

        contig = Contig()
        #seq0
        seq = 'TAGAT'
        contig.append(seq)
        #seq1
        seq = 'CTAT'
        contig.append_to_location(seq, start=2)
        #seq2
        seq = 'tTCTag'
        mask = (1, 3)
        contig.append_to_location(seq, mask=mask)
        #seq3
        seq = 'TAGAT'
        contig.append_to_location(seq, strand=-1)
        #seq4
        seq = 'cgAGc'
        mask = (2, 3)
        contig.append_to_location(seq, start=1, strand=-1, mask=mask)

        assert len(contig) == 5
        assert contig.ncols == 6
        assert contig[1, 2] == 'C'
        self.assertRaises(IndexError, contig.__getitem__, (1, 0))
        #if we have located a sequence inside the contig with an slice we'll
        #get a LocatableSequences
        assert contig[0, 1:3].sequence == 'AG'
        assert contig[1, 1:3].sequence == 'C'

        #now a contig composed by seqwithquality
        #ref   012345
        #seq1    ATCT fwd
        #qual1   5678
        #seq2  GCTA   revcomp
        #qual2 89123
        contig = Contig()
        seq1 = SeqRecord(Seq('ATCT'), [5, 6, 7, 8])
        contig.append_to_location(seq1, start=2)
        seq2 = SeqRecord(Seq('ATCG'), [3, 2, 1, 9])
        contig.append_to_location(seq2, strand=-1, forward=False)

        assert len(contig) == 2
        assert contig.ncols == 6
        assert contig[1, 0].seq == 'C'
        assert contig[0, 1:4].sequence.seq == 'AT'
        assert contig[0, 1:4].sequence.qual[0] == 5

        #now we get a column
        col = contig[:, 2]
        #the seqwithqual
        assert col.seq == 'AA'
        assert str(col.qual) == str([5, 2])
        #only the seq
        assert contig[:, 0].seq == 'C'
        #only the qual
        assert str(contig[:, 5].qual) == str([8])

        #let's slice an assembly
        #      0123456789
        #seq1  ATCT
        #seq2     TGACA
        #seq3         ACT
        contig2 = Contig()
        contig2.append('ATCT')
        contig2.append_to_location('TGACA', start=3)
        contig2.append_to_location('ACT', start=7)
        contig21 = contig2[:, 1:7]
        assert len(contig21) == 2
        assert contig21.ncols == 6
        assert contig21[0, 0] == 'T'
        assert contig21[1, 2] == 'T'

    def test_consensus(self):
        '''It tests that a consensus can be added.'''
        contig1 = Contig(consensus='ACGT')
        assert contig1.consensus == 'ACGT'

        #now we slice the contig
        contig2 = contig1[:, :2]
        assert contig2.consensus == 'AC'

    def test_empty_seqs_slice(self):
        '''We want also the empty seqs, could we get them?'''
        seq = [locate_sequence('ACTG', location=5), 'AAAA', 'TTTTTT']
        contig = Contig(seq)
        # Slice slice combination
        #default behaviour
        new_contig = contig[:, 0:1]
        assert len(new_contig) == 2
        #now we want the empty seqs
        contig.return_empty_seq = True
        new_contig = contig[:, 0:1]
        assert new_contig[0] is None
    
        # int slice combination
        contig.return_empty_seq = True
        new_contig = contig[0, 2:3]
        assert new_contig is None
        
    def test_empty_seqs_int(self):
        '''We want also the empty seqs, could we get them?'''
        seq = [locate_sequence('ACTG', location=5), 'AAAA', 'TTTTTT']
        contig = Contig(seq)        

        # These are with slice int combination
        #default behaviour
        contig.return_empty_seq = False
        new_contig = contig[:, 1]
        assert len(new_contig) == 2
        #now we want the empty seqs
        contig.return_empty_seq = True
        new_contig = contig[:, 0]
        assert new_contig[0] is None
        
        #These are with int int combination
        #default behaviour
        contig.return_empty_seq = False
        try:
            #pylint: disable-msg=W0104
            #the statement does has an effect
            contig[0, 1]
            self.fail("Index Error expected")
        except IndexError:
            pass
        
    @staticmethod
    def test_consensus_empty():
        ''' We test here if the contig watches what returns the contig when 
        you ask a column that doesn't exist'''
        consensus = locate_sequence('TT', location=0)
        contig = Contig(sequences=['AAAA'], consensus = consensus)
        assert  contig[:, 4:4].consensus  == None

    @staticmethod
    def xtest_contig_str():
        'It checks the str contig representation'
        consensus = locate_sequence('AA', location=1)
        contig = Contig(sequences=['xAAx'], consensus = consensus)
        assert str(contig) == '0              ->xAAx\nconsensus      -> AA\n'
 
def main():
    '''It runs the tests'''
    unittest.main()

if __name__ == '__main__':
    main()
