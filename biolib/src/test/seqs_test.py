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

import unittest
from biolib.seqs import SeqWithQuality, Seq

class SeqsTest(unittest.TestCase):
    '''Tests the seq with quality class '''
    #pylint: disable-msg=R0904

    @staticmethod
    def test_seqs_string():
        ''' Here we test it we can initialice a seq with quality
         and if we can print it. Seq is going to be a normal string'''
        
        #sequence1 = Seq('aaavvttt')
        
        # First we initialice the quality in the init 
        seq1 = SeqWithQuality(name = 'seq1', seq='aaavvttt', \
                               qual=[2, 4 , 1, 4, 5, 6, 12, 34])
        assert seq1

        # Here we add the quality after the initialization
        seq2 = SeqWithQuality(name = 'seq2', seq = 'aaavvttt')
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
        seq1 = SeqWithQuality(name = 'seq1', seq=sequence1, \
                               qual=[2, 4 , 1, 4, 5, 6, 12, 34])
        assert seq1
        
        # Here we add the quality after the initialization
        seq2 = SeqWithQuality(name = 'seq2', seq = sequence1)
        seq2.qual = [2, 4 , 1, 4, 5, 6, 12, 34]
        assert seq2 
        
        # We check if the seq can be complemented
        seq3 = seq2.complement()
        assert seq3.seq == 'tttggaaa'
    
        # Let's check the sliceability of the seq
        seq4 = []
        for num in range(len(seq1)):
            seq4.append(seq1[num])
        assert seq4[0].seq  == 'a'
        assert seq4[0].qual == [2]

    @staticmethod
    def test_add_seqs():
        ''' It checks if the seqs can be joined in a unique seq'''
        sequence1 = Seq('aaaccttt')
        
        # First we initialice the quality in the init 
        seq1 = SeqWithQuality(name = 'seq1', seq = sequence1, \
                               qual = [2, 4 , 1, 4, 5, 6, 12, 34])
        seq2 = SeqWithQuality(name = 'seq1', seq = sequence1, \
                               qual = [2, 4 , 1, 4, 5, 6, 12, 34])
        
        seq3 = seq1 + seq2
        assert seq3.seq == 'aaacctttaaaccttt'
        assert seq3.qual == [2, 4 , 1, 4, 5, 6, 12, 34, 2, 4 , 1, 4, 5, \
                             6, 12, 34]

    def test_empty_seq(self):
        'EmptySeq basic behaviour'
        seq = SeqWithQuality(length=300)
        assert len(seq) == 300

        #We can create an empty seq and add the sequence latter
        seq = SeqWithQuality(length=4)
        seq.seq = 'ACTG'

        #we'll fail if we try to reassing
        try:
            seq.seq = 'ACTG'
            self.fail('AttributeError expected')
            #pylint: disable-msg=W0704
        except AttributeError:
            pass

    @staticmethod
    def test_description():
        'A SeqWithQuality can have description and annotations'
        desc = 'a short sequence'
        annots = {'type':'region', 'go':['0001', '0002'], 'database':'my'}
        seq = SeqWithQuality(seq='A', description=desc, annotations=annots)
        assert seq.description is desc
        assert seq.annotations is annots

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
    def getitem():
        'It tests the get item method'
        seq = Seq('ACTG')
        seq2 = seq[1:3]
        assert seq2 == 'AC'
        assert seq2.complement() #is still a Seq
        assert seq[::-1].complement() #is still a Seq

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
