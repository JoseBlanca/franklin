'''
Created on 2009 mar 11

@author: peio
'''
from __future__ import with_statement
import unittest
from biolib.cafparser import CafFile

from caf_example_file import ExampleCafFile

class Test(unittest.TestCase):
    ''' It tests '''
        
    def test_reads(self):
        ''' we check if we can take the reads from the caf file'''
        with ExampleCafFile() as fname:
            caf_parser  = CafFile(fname)
            num_reads = 0
            for read in caf_parser.reads():
                num_reads  += 1
                assert read
            # We check if our read method finds all reads
            # in caf example file (21)
            assert num_reads ==   21  

    def test_contigs(self):
        ''' It checks if the contig method returns contigs'''
        with ExampleCafFile() as fname:
            caf_parser  = CafFile(fname)

            num_contig = 0
            for contig in caf_parser.contigs():
                num_contig += 1
                assert contig
            # We check if our contig method finds all contigs in
            # caf example file  (1) 
            assert num_contig == 1 
        
    def test_check_consensus_seq(self):
        ''' It checks if we get the consensus seq. It only checks the first 10 
        nucleotides TTCAAGCGAT and the quality '''
        
        with ExampleCafFile() as fname:
            caf_parser  = CafFile(fname)
            for contig in caf_parser.contigs():
                print contig
                assert  contig.consensus[:10] == 'TTCAAGCGAT'

    
    def xtest_long_index(self):
        ''' It checks if we can open a long test'''
        long_caf_file_name         = '/home/peio/eucalyptus_out.caf'
        long_caf_file      = CafFile(long_caf_file_name)
        
        
if __name__ == "__main__":
    unittest.main()
