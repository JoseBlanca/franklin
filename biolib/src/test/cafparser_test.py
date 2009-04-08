'''
Created on 2009 mar 11

@author: peio
'''
from __future__ import with_statement
import unittest
from biolib.cafparser import CafParser
import biolib
import os.path

DATA_DIR = os.path.join(os.path.split(biolib.__path__[0])[0], 'data')

class Test(unittest.TestCase):
    ''' It tests '''
    
    @staticmethod
    def test_reads():
        ''' we check if we can take the reads from the caf file'''
        fname = os.path.join(DATA_DIR, 'example.caf')
        caf_parser  = CafParser(fname)
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
        fname = os.path.join(DATA_DIR, 'example.caf')
        caf_parser  = CafParser(fname)
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
        fname = os.path.join(DATA_DIR, 'example2.caf')
        caf_parser  = CafParser(fname)
        num_contig = 0
        for contig in caf_parser.contigs():
            num_contig += 1
            assert contig
        # We check if our contig method finds all contigs in
        # caf example file  (1) 
        assert num_contig == 1
         
    @staticmethod   
    def xtest_check_consensus_seq():
        ''' It checks if we get the consensus seq. It only checks the first 10 
        nucleotides TTCAAGCGAT and the quality '''
        
        fname = os.path.join(DATA_DIR, 'example.caf')
        caf_parser  = CafParser(fname)
        for contig in caf_parser.contigs():
            assert contig.consensus.sequence.seq[:10] == 'TTCAAGCGAT'
    
    @staticmethod
    def xtest_long_index():
        ''' It checks if we can open a long test'''
        long_caf_file_name = '/home/peio/eucalyptus_out.caf'
        long_caf_file      = CafParser(long_caf_file_name)
        assert long_caf_file
        
        
if __name__ == "__main__":
    unittest.main()
