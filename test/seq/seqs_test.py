'''
Created on 2009 mar 27

@author: peio
'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of franklin.
# franklin is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# franklin is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with franklin. If not, see <http://www.gnu.org/licenses/>.

import unittest
from franklin.seq.seqs import SeqWithQuality, Seq, SeqFeature, get_seq_name
from Bio.SeqFeature import FeatureLocation, ExactPosition
from Bio.Alphabet import Alphabet

class SeqsTest(unittest.TestCase):
    '''Tests the seq with quality class '''
    #pylint: disable-msg=R0904

    @staticmethod
    def test_seqs_string():
        ''' Here we test it we can initialice a seq with quality
         and if we can print it. Seq is going to be a normal string'''

        #sequence1 = Seq('aaavvttt')

        # First we initialice the quality in the init
        seq1 = SeqWithQuality(name='seq1', seq=Seq('aaavvttt'),
                               qual=[2, 4 , 1, 4, 5, 6, 12, 34])
        assert seq1

        # Here we add the quality after the initialization
        seq2 = SeqWithQuality(name='seq2', seq=Seq('aaavvttt'))
        seq2.qual = [2, 4 , 1, 4, 5, 6, 12, 34]
        assert seq2

        # Let's check the sliceability of the seq
        for num in range(len(seq1)):
            seq3 = seq1[num]
            assert seq3

    @staticmethod
    def test_seq_seq():
        ''' We are going to check the same tests but with a BioPython seq
        class object inside seq'''
        sequence1 = Seq('aaaccttt')

        # First we initialice the quality in the init
        seq1 = SeqWithQuality(name='seq1', seq=sequence1, \
                               qual=[2, 4 , 1, 4, 5, 6, 12, 34])
        assert seq1

        # Here we add the quality after the initialization
        seq2 = SeqWithQuality(name='seq2', seq=sequence1)
        seq2.qual = [2, 4 , 1, 4, 5, 6, 12, 34]
        assert seq2

        # We check if the seq can be complemented
        seq3 = seq2.complement()
        assert seq3.seq == 'tttggaaa'

        # Let's check the sliceability of the seq
        assert seq2[0:2].seq == 'aa'
        assert seq2[0:2].qual == [2, 4]

    @staticmethod
    def test_add_seqs():
        ''' It checks if the seqs can be joined in a unique seq'''
        sequence1 = Seq('aaaccttt')

        # First we initialice the quality in the init
        seq1 = SeqWithQuality(name='seq1', seq=sequence1,
                              qual=[2, 4 , 1, 4, 5, 6, 12, 34])
        seq2 = SeqWithQuality(name='seq2', seq=sequence1,
                              qual=[2, 4 , 1, 4, 5, 6, 12, 34])
        seq3 = seq1 + seq2
        assert seq3.seq == 'aaacctttaaaccttt'
        assert seq3.qual == [2, 4 , 1, 4, 5, 6, 12, 34, 2, 4 , 1, 4, 5,
                             6, 12, 34]
        assert seq3.name == 'seq1+seq2'

    @staticmethod
    def test_description():
        'A SeqWithQuality can have description and annotations'
        desc = 'a short sequence'
        annots = {'type':'region', 'go':['0001', '0002'], 'database':'my'}
        seq = SeqWithQuality(seq=Seq('A'), description=desc, annotations=annots)
        assert seq.description is desc
        assert seq.annotations is annots


    @staticmethod
    def test_repr():
        'It test the __repr__ function'
        desc = 'a short sequence'
        annots = {'type':'region', 'go':['0001', '0002'], 'database':'my'}

        qual_list = [2, 4 , 1, 4, 5, 6, 12, 34]
        seq1 = SeqWithQuality(name='seq1', seq=Seq('aaaccttt'),
                              description=desc, annotations=annots,
                              qual=qual_list,)
        seqfeature = SeqFeature(location=FeatureLocation(5, 8), type='ortholog',
                                qualifiers={'arabidposys':['arab1', 'arab2']})
        seq1.features.append(seqfeature)

        a = eval(repr(seq1))
        assert  a.annotations == {'go': ['0001', '0002'], 'type': 'region',
                                  'database': 'my'}
        assert a.seq == 'aaaccttt'
        assert a.features[0].qualifiers == {'arabidposys': ['arab1', 'arab2']}
        assert a.qual == qual_list

    @staticmethod
    def test_upper():
        'It test the uppercase funtion'
        seq = SeqWithQuality(Seq('actg'))
        seq.upper()
        assert seq.seq == 'actg'

class SeqTest(unittest.TestCase):
    'It tests the Seq object.'
    @staticmethod
    def test_complement():
        'It test the Seq complement method'
        seq = Seq('ACTG')
        seq2 = seq.complement()
        assert seq2 == 'TGAC'

    @staticmethod
    def test_add():
        'It tests the add method'
        seq = Seq('AC') + Seq('TG')
        assert seq == 'ACTG'
        assert seq.complement() #is still a Seq

    @staticmethod
    def test_getitem():
        'It tests the get item method'
        seq = Seq('ACTG')
        seq2 = seq[1:3]
        assert str(seq2) == 'CT'
        assert seq2.complement() #is still a Seq
        assert seq[::-1].complement() #is still a Seq

    @staticmethod
    def test_repr():
        'It tests that the repr works ok'
        text = 'A' * 100
        seq = Seq(text)
        assert text in repr(seq)

class TestGetSeqName(unittest.TestCase):
    'It tests that we can get a sequence name'
    @staticmethod
    def test_get_seq_name():
        'It tests that we can get a sequence name'
        #with no name attribute -> uuid
        assert len(get_seq_name('AA')) > 10
        seq = SeqWithQuality(id='seqid', name='seqname', seq=Seq('ATGAT'))
        assert get_seq_name(seq) == 'seqname'
        seq = SeqWithQuality(id='seqid', seq=Seq('ATGAT'))
        assert get_seq_name(seq) == 'seqid'
        seq = SeqWithQuality(seq=Seq('ATGAT'))
        assert  len(get_seq_name(seq)) > 10

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
