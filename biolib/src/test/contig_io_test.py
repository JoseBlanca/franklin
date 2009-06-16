'''
Created on 2009 mar 11

@author: peio
'''
import unittest
from biolib.contig_io import get_parser
import biolib
import os.path
from test_utils import SeqWithQuality

DATA_DIR = os.path.join(os.path.split(biolib.__path__[0])[0], 'data')

class CafTest(unittest.TestCase):
    ''' It tests '''
    
    @staticmethod
    def test_reads():
        ''' we check if we can take the reads from the caf file'''
        fhand = open(os.path.join(DATA_DIR, 'example.caf'), 'r')
        caf_parser = get_parser(fhand, format='caf')
        num_reads = 0
        for read in caf_parser.reads():
            num_reads  += 1
            assert read
        # We check if our read method finds all reads
        # in caf example file (21)
        assert num_reads ==   21  

    @staticmethod
    def test_contigs():
        ''' It checks if the contig method returns contigs'''
        fhand = open(os.path.join(DATA_DIR, 'example.caf'), 'r')
        caf_parser = get_parser(fhand, format='caf')
        num_contig = 0
        for contig in caf_parser.contigs():
            num_contig += 1
            assert contig
        # We check if our contig method finds all contigs in
        # caf example file  (1) 
        assert num_contig == 1 
   
    @staticmethod
    def test_contigs2():
        ''' It checks if the contig method returns contigs, '''
        fhand = open(os.path.join(DATA_DIR, 'example2.caf'), 'r')
        caf_parser = get_parser(fhand, format='caf')
        num_contig = 0
        for contig in caf_parser.contigs():
            num_contig += 1
            assert contig
        # We check if our contig method finds all contigs in
        # caf example file  (1) 
        assert num_contig == 1
         
    @staticmethod   
    def test_check_consensus_seq():
        ''' It checks if we get the consensus seq. It only checks the first 10 
        nucleotides TTCAAGCGAT and the quality '''
        
        fhand = open(os.path.join(DATA_DIR, 'example.caf'), 'r')
        caf_parser = get_parser(fhand, format='caf')
        for contig in caf_parser.contigs():
            assert str(contig.consensus.sequence.seq[:10]) == 'TTCAAGCGAT'
 
    @staticmethod
    def xtest_long_index():
        ''' It checks if we can open a long test'''
        long_caf_file_name = '/home/peio/eucalyptus_out.caf'
        long_caf_file      = CafParser(long_caf_file_name)
        assert long_caf_file
    
    @staticmethod
    def test_alignement_seq():
        '''It checks if we locate the reads in good coordinates'''
        fhand = open(os.path.join(DATA_DIR, 'example3.caf'), 'r')
        caf_parser = get_parser(fhand, format='caf')
        for contig  in caf_parser.contigs():
            #consensus start
            assert contig.consensus.location.start == 5
            # These are reverse strands
            assert contig[204].location.start == 97
            assert contig[198].location.start == 2
            #these are forward strands
            assert contig[76].location.start  == 1
            assert contig[19].location.start  == 36
    @staticmethod
    def test_read_seq():
        '''It checks if the reads have the correct secuence. We will check
        randomonly selected columns'''
        fhand = open(os.path.join(DATA_DIR, 'example3.caf'), 'r')
        caf_parser = get_parser(fhand, format='caf')
        for contig in caf_parser.contigs():
            #forward
            assert str(contig[83][46:56]) == 'GGCCGGG-GC'
            assert str(contig[19][41:51]) == 'TTA-CGGCCG'
            #reverses
            read =  contig[204][139:149]
            assert str(read) == 'ATCCACTTTT'
            read = contig[198][8:18]
            assert str(read) == 'CTCCCTGTGN'
 
class AceTest(unittest.TestCase):
    ''' It tests the ace alignment parser'''
    
    def test_contig(self):
        '''It tests that we can get a read by its name.'''
        fhand = open(os.path.join(DATA_DIR, 'example3.ace'), 'r')
        ace_parser = get_parser(fhand, format='ace')
        #if we ask for a wrong contig we get an error
        try:
            ace_parser.contig('not_in_file')
            self.fail('KeyError expected')
            #pylint: disable-msg=W0704
        except ValueError:
            pass
        #now for a real contig
        contig = ace_parser.contig('eucalyptus_lrc1')
        assert contig
 
    @staticmethod   
    def test_check_consensus_seq():
        ''' It checks if we get the consensus seq. It only checks the first 10 
        nucleotides TTCAAGCGAT and the quality '''
 
        fhand = open(os.path.join(DATA_DIR, 'example.caf'), 'r')
        ace_parser = get_parser(fhand, format='ace')
        for contig in ace_parser.contigs():
            assert str(contig.consensus.sequence.seq[:10]) == 'TTCAAGCGAT'

    @staticmethod
    def test_contigs():
        '''It checks if the contig method returns contigs'''
        fhand = open(os.path.join(DATA_DIR, 'example.ace'), 'r')
        ace_parser = get_parser(fhand, format='ace')
        num_contig = 0
        for contig in ace_parser.contigs():
            num_contig += 1
            assert contig
        # We check if our contig method finds all contigs in
        # caf example file  (1) 
        assert num_contig == 8

    @staticmethod
    def test_read_seq():
        '''It checks if the reads have the correct secuence. We will check
        randomonly selected columns'''
        fhand = open(os.path.join(DATA_DIR, 'example3.caf'), 'r')
        ace_parser = get_parser(fhand, format='ace')
        for contig in ace_parser.contigs():
            for read in contig:
                if read.sequence.name == 'E3MFGYR01AWFJG':
                    assert str(read[46:56]) == 'GGCCGGG-GC'
                if read.sequence.name == 'E3MFGYR02JMR2C':
                    assert str(read[41:51]) == 'TTA-CGGCCG'

class CigarTest(unittest.TestCase):
    'Write and read cigar.'
    seq1 = 'AAGCTCANCTTGGACCACCGACTCTCGANTGNNTCGCCGCGGGAGCCGGNTGGANAACCT'
    seq1 = SeqWithQuality(seq1, name='hs989235.cds')
    seq1 = locate_sequence(seq1, location=6)
    seq2 = 'AAGCTCATCTTGG-CCACCGACTCTCGCTTGCGCCGCCGCGGGAGCCGG-TGGA-AACCT'
    seq2 = SeqWithQuality(seq2, name='hsnfg9.embl')
    seq2 = locate_sequence(seq2, location=25690)
    contig = Contig([seq1, seq2])
    print contig.format('cigar')


if __name__ == "__main__":
    unittest.main()
