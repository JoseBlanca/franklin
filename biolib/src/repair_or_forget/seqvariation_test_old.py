'''
Created on 2009 mar 25

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
from biolib.seqvariation import (CONFIG, SeqVariation,
                                 seqvars_in_contigs,
                                 calculate_pic,
                                 cap_enzime)
from biolib.locatable_sequence import locate_sequence
from biolib.contig import Contig
from biolib.seqs import SeqRecord
from test.test_utils import SeqRecord

class SeqVariationConfigTest(unittest.TestCase):
    '''Here we will check if the SeqVariation module can be configured.
    '''
    @staticmethod
    def test_config():
        '''It checks that we can configure the module.'''
        assert isinstance(CONFIG.min_num_of_reads, int)
        assert isinstance(CONFIG.only_snp, bool)
        indel_char = CONFIG.indel_char
        assert isinstance(indel_char, str)
        assert len(indel_char) == 1

class SeqVariationTest(unittest.TestCase):
    '''Here we will check if the SeqVariation module works as it should.

    Is the support for snps and indels.
    '''
    @staticmethod
    def test_init():
        '''It tests the init'''
        snp = SeqVariation(name='a snp', location=20,
                           alignment='some_contig_instance',
                           alleles={'A':2, 'T':3})
        assert snp.name == 'a snp'
        assert snp.location == 20
        assert snp.alignment == 'some_contig_instance'
        assert len(snp.alleles) == 2

        #a simple init
        snp = SeqVariation(alleles={'A':2, 'T':3})
        assert len(snp.alleles) == 2

        #some alleles are removed due to the min_number_of_reads
        CONFIG.min_num_of_reads = 3
        snp = SeqVariation(alleles={'A':2, 'T':3})
        assert len(snp.alleles) == 1

        #some alleles are removed due to only_snps
        CONFIG.only_snp = True
        indel = CONFIG.indel_char
        snp = SeqVariation(alleles={'A':3, indel:3})
        assert len(snp.alleles) == 1
        CONFIG.only_snp = False

    @staticmethod
    def test_init_with_seq_lists():
        '''It test that SeqVariation also work with {'A':[seq1], 'T':[seq2]}'''
        CONFIG.min_num_of_reads = 2
        snp = SeqVariation(alleles={'A':[1, 2], 'T':['seq1', 'seq2'], '-':[1]})
        assert len(snp.alleles) == 2

    @staticmethod
    def test_kind():
        '''it checks that we can get the kind of variation.'''
        inchar = CONFIG.indel_char

        seq_var = SeqVariation(alleles={'A':3, inchar:3})
        assert seq_var.is_indel()

        seq_var = SeqVariation(alleles={'A':[1, 2, 3], inchar:[4, 5, 6]})
        assert seq_var.is_indel()

        seq_var = SeqVariation(alleles={'AA':3, inchar * 2:3})
        assert seq_var.is_indel()
        seq_var = SeqVariation(alleles={'AA':[1, 2, 3], inchar * 2:[4, 5, 6]})
        assert seq_var.is_indel()

        seq_var = SeqVariation(alleles={'AA' + inchar :3, inchar * 3:3})
        assert seq_var.is_complex()

        seq_var = SeqVariation(alleles={'A':3, 'T':4})
        assert seq_var.is_snp()

        seq_var = SeqVariation(alleles={'A':3, 'T':4, inchar:3})
        assert seq_var.is_complex()

        seq_var = SeqVariation(alleles={'A' + inchar + 'A':3,
                                        'A' + inchar * 2 :4})
        assert seq_var.is_complex()

    @staticmethod
    def test_sorted_alleles():
        'It checks that we can get the alleles sorted by the reads.'
        snp = SeqVariation(alleles={'A':2, 'T':3})
        sorted_alleles = snp.sorted_alleles()
        assert sorted_alleles[0][0] == 'T'
        assert sorted_alleles[0][1] == 3
        assert sorted_alleles[1][0] == 'A'
        assert sorted_alleles[1][1] == 2

        #with lists
        snp = SeqVariation(alleles={'A':(1, 2), 'T':(2, 3, 4)})
        sorted_alleles = snp.sorted_alleles()
        assert sorted_alleles[0][0] == 'T'
        assert sorted_alleles[0][1] == (2, 3, 4)
        assert sorted_alleles[1][0] == 'A'
        assert sorted_alleles[1][1] == (1, 2)

def check_alleles(expected_alleles, contig):
    '''It checks that all the alleles for all variations are found'''
    detected_vars = 0
#    for seqvar in seqvariations_in_alignment(contig):


    for seqvar in seqvars_in_contigs([contig]):
        expected = expected_alleles[detected_vars]
        found = seqvar.alleles.keys()
        diff = set(expected).difference(found)
        assert not diff
        detected_vars += 1
    assert detected_vars == len(expected_alleles)

class VariationGeneratorTest(unittest.TestCase):
    '''It test if we can get the variations'''
    @staticmethod
    def test_seq_variation_generator():
        '''It checks if we can get the variations.'''

        CONFIG.min_num_of_reads = 2
        #can we find snp?
        allele1 = 'AtA'
        allele2 = 'CtT'
        seqs = [SeqRecord(allele1), SeqRecord(allele2), SeqRecord(allele1),
                SeqRecord(allele2)]
        contig = Contig(sequences = seqs)
        expected_alleles = (('A','C'), ('A','T'))
        check_alleles(expected_alleles, contig)

        #can we find indels with length 1?
        allele1 = '-t-'
        allele2 = 'CtT'
        seqs = [SeqRecord(allele1), SeqRecord(allele2), SeqRecord(allele1),
                SeqRecord(allele2)]
        contig = Contig(sequences = seqs)
        expected_alleles = (('-','C'), ('-','T'))
        check_alleles(expected_alleles, contig)

        #can we find indels with length 1?
        allele1 = '-tn'
        allele2 = 'CtT'
        seqs = [SeqRecord(allele1), SeqRecord(allele2), SeqRecord(allele1),
                SeqRecord(allele2)]
        contig = Contig(sequences = seqs)
        expected_alleles = (('-','C'),)
        check_alleles(expected_alleles, contig)

        allele1 = SeqRecord('--')
        allele2 = SeqRecord('-a')
        allele3 = SeqRecord('aa')
        seqs = [allele1, allele2, allele3, allele3, allele3]
        contig = Contig(sequences = seqs)
        expected_alleles = ()
        check_alleles(expected_alleles, contig)


        allele1 = SeqRecord('AAA')
        allele2 = SeqRecord('--A')
        allele3 = SeqRecord('A--')
        seqs = [allele1, allele1, allele2, allele2, allele3, allele3]
        contig = Contig(sequences = seqs)
        expected_alleles = (('AAA', '--A', 'A--'), )
        check_alleles(expected_alleles, contig)

        allele1 = SeqRecord('AAA')
        allele2 = SeqRecord('A--')
        allele3 = SeqRecord('--A')
        allele4 = SeqRecord('-AA')
        seqs = [allele1, allele1, allele2, allele2, allele3, allele4]
        contig = Contig(sequences = seqs)
        expected_alleles = (('AA','--'), )
        check_alleles(expected_alleles, contig)

        allele1 = SeqRecord('AAA')
        seqs = [allele1, allele1, allele1]
        contig = Contig(sequences = seqs)
        expected_alleles = (() )
        check_alleles(expected_alleles, contig)

    @staticmethod
    def test_seq_variation_with_locseq():
        '''It checks if we can get the variations using LocatableSequence.'''
        # here we check if it works when the seqs in the contigs are
        # locatable sequences
        allele1 = locate_sequence(sequence=SeqRecord('AA'), location=0)
        allele2 = locate_sequence(sequence=SeqRecord('TA'), location=1)
        seqs = [allele1, allele1, allele2, allele2]
        contig = Contig(sequences = seqs)
        expected_alleles = (('AT'), )
        check_alleles(expected_alleles, contig)

class SeqVariationCaracterization(unittest.TestCase):
    '''It tests seqvar caracterization functions  '''
    @staticmethod
    def test_calculate_pic():
        'It checks that we are able to calculate the PIC values'
        snp1 = SeqVariation(alleles={'A':200, 'T':200})
        snp2 = SeqVariation(alleles={'A':300, 'T':100})
        pic_1 = calculate_pic(snp1)
        pic_2 = calculate_pic(snp2)
        assert pic_1 > pic_2

class SeqVariationrEnzime(unittest.TestCase):
    ''' It checks if we have problems with remaps and it functions'''

    @staticmethod
    def test_remap():
        '''It test if the remap external program works '''
        con  = SeqRecord(seq='Actgactgactgtca', name='consensus')
        seq  = SeqRecord(seq='Actgacttactgtca', name='seq')
        consensus = locate_sequence(con, location=0)
        cont = Contig([seq, seq], consensus=consensus)
        snp = SeqVariation(alleles={'C':2, 'T':3}, location=7, alignment=cont)

        enzymes = cap_enzime(snp, True)
        assert ['HinfI', 'TscAI'] == enzymes

        # We need to test it with locations
        con       = SeqRecord(seq='Actgactgactgtca', name='consensus')
        seq       = SeqRecord(seq='Actgacttactgtca', name='seq')
        consensus = locate_sequence(con, location=0)
        contig    = Contig(consensus=consensus)
        contig.append_to_location(seq, 0)
        contig.append_to_location(seq, 0)
        snp = SeqVariation(alleles={'C':2, 'T':3}, location=(7, 8),
                           alignment=contig)
        enzymes = cap_enzime(snp, True)
        assert ['HinfI', 'TscAI'] == enzymes


    @staticmethod
    def test_ecori():
        '''It test a know enzyme reaction ecori '''

        seq_eco   = 'gaattc'
        seq_noeco = 'gaatcc'

        con  ='ATGATGATG' + seq_eco + 'ATGATGATGTGGGAT'
        seq1 = 'ATGATGATG' + seq_eco + 'ATGATGATGTGGGAT'
        seq2 = 'ATGATGATG' + seq_noeco + 'ATGATGATGTGGGAT'
        cont = Contig([seq1, seq2], consensus=con)

        snp = SeqVariation(alleles={'C':2, 'T':3}, location=13, alignment=cont)

        enzymes = cap_enzime(snp)
        assert 'EcoRI' in enzymes



if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_SeqVariation_init']
    unittest.main()
