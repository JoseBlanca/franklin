'''
Created on 2009 mar 25

@author: peio
'''
import unittest
from biolib.SeqVariation import (CONFIG, SeqVariation,
                                 seqvariations_in_alignment,
                                 second_allele_read_times)
from biolib.contig import Contig, locate_sequence
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
    for seqvar in seqvariations_in_alignment(contig):
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
        

class SeqVariationFilteringTest(unittest.TestCase):
    'It checks the filtering methods.'
    @staticmethod
    def test_filter_by_number_reads():
        'It checks that filter the reads with not much reads.'
        #It checks the number of times that the second most abundant allele 
        #has been read
        snp = SeqVariation(alleles={'A':2, 'T':3, 'C':4})
        assert second_allele_read_times(snp, times=3) == True
        assert second_allele_read_times(snp, times=4) == False


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_SeqVariation_init']
    unittest.main()
