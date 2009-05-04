'''
Created on 2009 mar 25

@author: peio
'''
import unittest
from biolib.SeqVariation import (CONFIG, SeqVariation,
                                 seqvariations_in_alignment,
                                 second_allele_read_times,
                                 remove_bad_quality_alleles,
                                 seqvar_close_to_consensus_limit,
                                 calculate_pic,
                                 cap_enzime)
from biolib.contig import Contig, locate_sequence, Location
from biolib.Seqs import SeqWithQuality
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
        snp2 = remove_bad_quality_alleles(snp, qual_threshold=25)
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
        snp2 = remove_bad_quality_alleles(snp, qual_threshold=25)
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
        snp2 = remove_bad_quality_alleles(snp, qual_threshold=25, \
                                          default_quality=25)
        assert len(snp2.alleles) == 3
        
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
        self.failUnlessRaises(ValueError, remove_bad_quality_alleles, \
                              snp, qual_threshold=25)
        
        # The  gaps does not hace quality , so we ignore the quelity but we use
        # them 
        seq1 = SeqWithQuality(seq='AA-AA', qual=[30, 30, 0, 30, 30])
        seq2 = SeqWithQuality(seq='AACAA', qual=[30, 30, 30, 30, 30])
        seq3 = SeqWithQuality(seq='AAGAA', qual=[30, 30, 15, 30, 30])
        contig = Contig([seq1, seq1, seq2, seq2, seq3, seq3])
        snp = SeqVariation(location=2, alignment=contig,
                           alleles={'-':[0, 1], 'C':[2, 3], 'G':[4, 5]})
        snp2 = remove_bad_quality_alleles(snp, qual_threshold=25)
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
        snp2 = remove_bad_quality_alleles(snp, qual_threshold=25)
        assert snp2 is None        
        
        # It checks indels with more than one column
        seq1 = SeqWithQuality(seq='AA--AA', qual=[30, 30, 0, 0, 30, 30])
        seq2 = SeqWithQuality(seq='AACCAA', qual=[30, 30, 30, 30, 30, 30])
        seq3 = SeqWithQuality(seq='AAGGAA', qual=[30, 30, 15, 15, 30, 30])
        contig = Contig([seq1, seq1, seq2, seq2, seq3, seq3])
        
        snp = SeqVariation(location=Location(start=2, end=3), alignment=contig,
                           alleles={'--':[0, 1], 'CC':[2, 3], 'GG':[4, 5]})
        snp2 = remove_bad_quality_alleles(snp, qual_threshold=25)
        assert len(snp2.alleles) == 2
        
        #It check a complex SeqVar
        seq1 = SeqWithQuality(seq='AA--AA', qual=[30, 30, 0, 0, 30, 30])
        seq2 = SeqWithQuality(seq='AAC-AA', qual=[30, 30, 30, 0, 30, 30])
        seq3 = SeqWithQuality(seq='AAGGAA', qual=[30, 30, 15, 15, 30, 30])
        contig = Contig([seq1, seq1, seq2, seq2, seq3, seq3])
        
        snp = SeqVariation(location=Location(start=2, end=3), alignment=contig,
                           alleles={'--':[0, 1], 'C-':[2, 3], 'GG':[4, 5]})
        snp2 = remove_bad_quality_alleles(snp, qual_threshold=25)
        assert len(snp2.alleles) == 2
        
        #It check a complex SeqVar
        seq1 = SeqWithQuality(seq='AA---AA', qual=[30, 30, 0, 0, 0, 30, 30])
        seq2 = SeqWithQuality(seq='AAC--AA', qual=[30, 30, 30, 0, 0, 30, 30])
        seq3 = SeqWithQuality(seq='AAGGGAA', qual=[30, 30, 15, 15, 15, 30, 30])
        seq4 = SeqWithQuality(seq='AAGGGAA', qual=[30, 30, 15, 30, 30, 30, 30])
        contig = Contig([seq1, seq1, seq2, seq2, seq3, seq4])
        
        snp = SeqVariation(location=Location(start=2, end=4), alignment=contig,
                           alleles={'---':[0, 1], 'C--':[2, 3], 'GGG':[4, 5]})
        snp2 = remove_bad_quality_alleles(snp, qual_threshold=25)
        assert len(snp2.alleles) == 2
        
    def test_close_to_consensus_end(self):
        'It checks the filtering out of the seqvars close to the end.'
        con = 'ATCTGACTT'
        seq = 'TT' + con + 'TT'
        consensus = locate_sequence(con, location=2)
        cont = Contig([seq, seq], consensus=consensus)
        snp = SeqVariation(alleles={'A':2, 'T':3}, location=1, alignment=cont)
        assert seqvar_close_to_consensus_limit(snp, max_distance=2)
        snp = SeqVariation(alleles={'A':2, 'T':3}, location=3, alignment=cont)
        assert seqvar_close_to_consensus_limit(snp, max_distance=2)
        snp = SeqVariation(alleles={'A':2, 'T':3}, location=4, alignment=cont)
        assert not seqvar_close_to_consensus_limit(snp, max_distance=2)
        snp = SeqVariation(alleles={'A':2, 'T':3}, location=10, alignment=cont)
        assert seqvar_close_to_consensus_limit(snp, max_distance=2)
        snp = SeqVariation(alleles={'A':2, 'T':3}, location=12, alignment=cont)
        assert seqvar_close_to_consensus_limit(snp, max_distance=2)

        #without LocatableSequences
        con = 'ATCTGACTT'
        seq = con
        cont = Contig([seq, seq], consensus=con)
        snp = SeqVariation(alleles={'A':2, 'T':3}, location=2, alignment=cont)
        assert seqvar_close_to_consensus_limit(snp, max_distance=3)
        assert not seqvar_close_to_consensus_limit(snp, max_distance=1)

        #some errors
        #no consensus
        #pylint: disable-msg=W0704
        cont = Contig([seq, seq])
        snp = SeqVariation(alleles={'A':2, 'T':3}, location=2, alignment=cont)
        try:
            seqvar_close_to_consensus_limit(snp, max_distance=3)
            self.fail('ValueError expected')
        except ValueError:
            pass
        #no contig
        snp = SeqVariation(alleles={'A':2, 'T':3}, location=2)
        try:
            seqvar_close_to_consensus_limit(snp, max_distance=3)
            self.fail('ValueError expected')
        except ValueError:
            pass

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
        con  = 'Actgactgactgtca'
        seq  = 'Actgacttactgtca'
        consensus = locate_sequence(con, location=0)
        cont = Contig([seq, seq], consensus=consensus)
        snp = SeqVariation(alleles={'C':2, 'T':3}, location=7, alignment=cont)
            
        enzimes = cap_enzime(snp, True)
        assert ['HinfI', 'TscAI'] == enzimes
        
    


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_SeqVariation_init']
    unittest.main()
