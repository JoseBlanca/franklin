'''
Created on 2009 uzt 30

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

from biolib.seqvariation import  SeqVariation, calculate_pic
from biolib.snp_cleaner import (create_second_allele_number_filter,
                                create_seqvar_close_to_limit_filter,
                                create_bad_quality_allele_remover,
                                create_pic_filter,
                                create_cap_enzyme_filter)

from biolib.locatable_sequence import locate_sequence, Location
from biolib.contig import Contig
from biolib.seqs import SeqWithQuality



class SeqVariationFilteringTest(unittest.TestCase):
    'It checks the filtering methods.'
    @staticmethod
    def test_filter_by_number_reads():
        'It checks that filter the reads with not much reads.'
        #It checks the number of times that the second most abundant allele
        #has been read
        snp = SeqVariation(alleles={'A':2, 'T':3, 'C':4})
        second_allele_read_times = create_second_allele_number_filter(3)
        assert second_allele_read_times(snp) == True

        second_allele_read_times = create_second_allele_number_filter(4)
        assert second_allele_read_times(snp) == False


    def test_remove_bad_quality_alleles(self):
        'It test it we can remove the alleles with bad qualities'
        #the good reads in 454 are around a quality of 25 to 40
        #a bad letter is less than 20
        #for sanger the phred quality values are related to the probability
        #of having an error in the sequence
        #http://www.phrap.com/phred/
        #Phred quality score  Probability base is wrong   Accuracy of the base
        #       10                   1 in 10                  90%
        #       20                   1 in 100                 99%
        #       30                   1 in 1,000               99.9%
        #       40                   1 in 10,000              99.99%
        #       50                   1 in 100,000             99.999%
        seq1 = SeqWithQuality(seq='AATAA', qual=[30, 30, 30, 30, 30])
        seq2 = SeqWithQuality(seq='AACAA', qual=[30, 30, 30, 30, 30])
        seq3 = SeqWithQuality(seq='AAGAA', qual=[30, 30, 15, 30, 30])
        contig = Contig([seq1, seq1, seq2, seq2, seq3, seq3])
        snp = SeqVariation(location=2, alignment=contig,
                           alleles={'T':[0, 1], 'C':[2, 3], 'G':[4, 5]})

        remove_bad_quality_alleles = create_bad_quality_allele_remover(25)
        snp2 = remove_bad_quality_alleles(snp)
        assert len(snp2.alleles) == 2


        #The same test but with locatable sequences
        seq1 = locate_sequence(sequence=SeqWithQuality(seq='AATAA', \
                                                     qual=[30, 30, 30, 30, 30]),
                               location=0)
        seq2 = locate_sequence(sequence=SeqWithQuality(seq='AACAA', \
                                                     qual=[30, 30, 30, 30, 30]),
                               location=0)
        seq3 = SeqWithQuality(seq='AAGAA', qual=[30, 30, 15, 30, 30])
        contig = Contig([seq1, seq1, seq2, seq2, seq3, seq3])
        snp = SeqVariation(location=2, alignment=contig,
                           alleles={'T':[0, 1], 'C':[2, 3], 'G':[4, 5]})
        snp2 = remove_bad_quality_alleles(snp)
        assert len(snp2.alleles) == 2

        # One of the sequences does not have quality, We do not provide
        # default quelity. It have to give and error
        seq1 = locate_sequence(sequence=SeqWithQuality(seq='AATAA', \
                                                     qual=[30, 30, 30, 30, 30]),
                               location=0)
        seq2 = locate_sequence(sequence=SeqWithQuality(seq='AACAA', \
                                                     qual=[30, 30, 30, 30, 30]),
                               location=0)
        seq3 = SeqWithQuality(seq='AAGAA' )
        contig = Contig([seq1, seq1, seq2, seq2, seq3, seq3])
        snp = SeqVariation(location=2, alignment=contig,
                           alleles={'T':[0, 1], 'C':[2, 3], 'G':[4, 5]})
        self.failUnlessRaises(ValueError, remove_bad_quality_alleles, snp)

        # The  gaps does not hace quality , so we ignore the quelity but we use
        # them
        seq1 = SeqWithQuality(seq='AA-AA', qual=[30, 30, 0, 30, 30])
        seq2 = SeqWithQuality(seq='AACAA', qual=[30, 30, 30, 30, 30])
        seq3 = SeqWithQuality(seq='AAGAA', qual=[30, 30, 15, 30, 30])
        contig = Contig([seq1, seq1, seq2, seq2, seq3, seq3])
        snp = SeqVariation(location=2, alignment=contig,
                           alleles={'-':[0, 1], 'C':[2, 3], 'G':[4, 5]})
        snp2 = remove_bad_quality_alleles(snp)
        assert len(snp2.alleles) == 2
        #if there is no quality we get an error

        # Now we have to check if it return a seqvar with only one allele after
        # the quality testing
        seq1 = SeqWithQuality(seq='AATAA', qual=[30, 30, 15, 30, 30])
        seq2 = SeqWithQuality(seq='AACAA', qual=[30, 30, 30, 30, 30])
        seq3 = SeqWithQuality(seq='AAGAA', qual=[30, 30, 15, 30, 30])
        contig = Contig([seq1, seq1, seq2, seq2, seq3, seq3])
        snp = SeqVariation(location=2, alignment=contig,
                           alleles={'T':[0, 1], 'C':[2, 3], 'G':[4, 5]})
        snp2 = remove_bad_quality_alleles(snp)
        assert snp2 is None

        # It checks indels with more than one column
        seq1 = SeqWithQuality(seq='AA--AA', qual=[30, 30, 0, 0, 30, 30])
        seq2 = SeqWithQuality(seq='AACCAA', qual=[30, 30, 30, 30, 30, 30])
        seq3 = SeqWithQuality(seq='AAGGAA', qual=[30, 30, 15, 15, 30, 30])
        contig = Contig([seq1, seq1, seq2, seq2, seq3, seq3])

        snp = SeqVariation(location=Location(start=2, end=3), alignment=contig,
                           alleles={'--':[0, 1], 'CC':[2, 3], 'GG':[4, 5]})
        snp2 = remove_bad_quality_alleles(snp)
        assert len(snp2.alleles) == 2

        #It check a complex SeqVar
        seq1 = SeqWithQuality(seq='AA--AA', qual=[30, 30, 0, 0, 30, 30])
        seq2 = SeqWithQuality(seq='AAC-AA', qual=[30, 30, 30, 0, 30, 30])
        seq3 = SeqWithQuality(seq='AAGGAA', qual=[30, 30, 15, 15, 30, 30])
        contig = Contig([seq1, seq1, seq2, seq2, seq3, seq3])

        snp = SeqVariation(location=Location(start=2, end=3), alignment=contig,
                           alleles={'--':[0, 1], 'C-':[2, 3], 'GG':[4, 5]})
        snp2 = remove_bad_quality_alleles(snp)
        assert len(snp2.alleles) == 2

        #It check a complex SeqVar
        seq1 = SeqWithQuality(seq='AA---AA', qual=[30, 30, 0, 0, 0, 30, 30])
        seq2 = SeqWithQuality(seq='AAC--AA', qual=[30, 30, 30, 0, 0, 30, 30])
        seq3 = SeqWithQuality(seq='AAGGGAA', qual=[30, 30, 15, 15, 15, 30, 30])
        seq4 = SeqWithQuality(seq='AAGGGAA', qual=[30, 30, 15, 30, 30, 30, 30])
        contig = Contig([seq1, seq1, seq2, seq2, seq3, seq4])

        snp = SeqVariation(location=Location(start=2, end=4), alignment=contig,
                           alleles={'---':[0, 1], 'C--':[2, 3], 'GGG':[4, 5]})
        snp2 = remove_bad_quality_alleles(snp)
        assert len(snp2.alleles) == 2

         #The same test but with locatable sequences. One of the sequences does
        # not have quality, We provide default quelity
        seq1 = locate_sequence(sequence=SeqWithQuality(seq='AATAA', \
                                                     qual=[30, 30, 30, 30, 30]),
                               location=0)
        seq2 = locate_sequence(sequence=SeqWithQuality(seq='AACAA', \
                                                     qual=[30, 30, 30, 30, 30]),
                               location=0)
        seq3 = SeqWithQuality(seq='AAGAA' )
        contig = Contig([seq1, seq1, seq2, seq2, seq3, seq3])
        snp = SeqVariation(location=2, alignment=contig,
                           alleles={'T':[0, 1], 'C':[2, 3], 'G':[4, 5]})
        remove_bad_quality_alleles = \
        create_bad_quality_allele_remover(qual_threshold=25,default_quality=25)

        snp2 = remove_bad_quality_alleles(snp)
        assert len(snp2.alleles) == 3


    def test_close_to_consensus_end(self):
        'It checks the filtering out of the seqvars close to the end.'
        con = 'ATCTGACTT'
        seq = 'TT' + con + 'TT'
        consensus = locate_sequence(con, location=2)
        cont = Contig([seq, seq], consensus=consensus)
        snp = SeqVariation(alleles={'A':2, 'T':3}, location=1, alignment=cont)
        seqvar_close_to_consensus_limit = create_seqvar_close_to_limit_filter(2)

        assert seqvar_close_to_consensus_limit(snp)
        snp = SeqVariation(alleles={'A':2, 'T':3}, location=3, alignment=cont)
        assert seqvar_close_to_consensus_limit(snp)
        snp = SeqVariation(alleles={'A':2, 'T':3}, location=4, alignment=cont)
        assert not seqvar_close_to_consensus_limit(snp)
        snp = SeqVariation(alleles={'A':2, 'T':3}, location=10, alignment=cont)
        assert seqvar_close_to_consensus_limit(snp)
        snp = SeqVariation(alleles={'A':2, 'T':3}, location=12, alignment=cont)
        assert seqvar_close_to_consensus_limit(snp)

        #without LocatableSequences
        con = 'ATCTGACTT'
        seq = con
        cont = Contig([seq, seq], consensus=con)
        snp = SeqVariation(alleles={'A':2, 'T':3}, location=2, alignment=cont)
        seqvar_close_to_consensus_limit = create_seqvar_close_to_limit_filter(3)
        assert seqvar_close_to_consensus_limit(snp)
        seqvar_close_to_consensus_limit = create_seqvar_close_to_limit_filter(1)
        assert not seqvar_close_to_consensus_limit(snp)

        #some errors
        #no consensus
        #pylint: disable-msg=W0704
        cont = Contig([seq, seq])
        snp = SeqVariation(alleles={'A':2, 'T':3}, location=2, alignment=cont)
        seqvar_close_to_consensus_limit = create_seqvar_close_to_limit_filter(3)
        try:
            seqvar_close_to_consensus_limit(snp)
            self.fail('ValueError expected')
        except ValueError:
            pass
        #no contig
        snp = SeqVariation(alleles={'A':2, 'T':3}, location=2)
        try:
            seqvar_close_to_consensus_limit(snp)
            self.fail('ValueError expected')
        except ValueError:
            pass

    @staticmethod
    def test_pic_filter():
        '''It test if we can filter by pic '''
        snp1 = SeqVariation(alleles={'A':200, 'T':200})
        snp2 = SeqVariation(alleles={'A':300, 'T':100})
        pic_filter = create_pic_filter(0.35)
        assert pic_filter(snp1) == True
        assert pic_filter(snp2) == False

    @staticmethod
    def test_enzyme_filter1():
        'It test the enzyme filter'
        con  = SeqWithQuality(seq='Actgactgactgtca', name='consensus')
        seq  = SeqWithQuality(seq='Actgacttactgtca', name='seq')
        consensus = locate_sequence(con, location=0)
        cont = Contig([seq, seq], consensus=consensus)
        snp = SeqVariation(alleles={'C':2, 'T':3}, location=7, alignment=cont)

        cap_enzime = create_cap_enzyme_filter(True)
        assert cap_enzime(snp) == True










if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_SeqVariation_init']
    unittest.main()