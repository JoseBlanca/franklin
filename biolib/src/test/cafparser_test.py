'''
Created on 2009 mar 11

@author: peio
'''
import unittest
from biolib.cafparser import CafFile

class Test(unittest.TestCase):
    ''' It tests '''
    def setUp(self):
        ''' It defines the files to use in the tests and creates 
        the index object to use in other tests '''
        self._example_caf_file_name = '../doc/Caf_example_file.caf'
        self._example_caf_file      = CafFile(self._example_caf_file_name)
        
        
    def test_index_creation(self):
        '''It checks if the indexes have been created '''
        caf_file = self._example_caf_file
        assert len(caf_file._qual_index) == 21 
        assert len(caf_file._seq_index)  == 22
        assert len(caf_file._dna_index)  == 22
        assert len(caf_file._type_index) == 22
           
    def test_reads(self):
        ''' we check if we can take the reads from the caf file'''
        caf_file = self._example_caf_file
        
        num_reads = 0
        for read in caf_file.reads():
            num_reads  += 1
            assert read
        # We check if our read method finds all reads
        # in caf example file (21)
        assert num_reads ==   21  

    def test_contigs(self):
        ''' It checks if the contig method returns contigs'''
        caf_file = self._example_caf_file

        num_contig = 0
        for contig in caf_file.contigs():
            num_contig += 1
            assert contig
#            print contig_object
        # We check if our contig method finds all contigs in
        # caf example file  (1) 
        assert num_contig == 1 
        
    def test_check_consensus_seq(self):
        ''' It checks if we get the consensus seq. It only checks the first 10 
        nucleotides TTCAAGCGAT and the quality '''
        
        caf_file = self._example_caf_file
        for contig in caf_file.contigs():
            print contig
            assert  contig.consensus[:10] == 'TTCAAGCGAT'

    
    def xtest_long_index(self):
        ''' It checks if we can open a long test'''
        long_caf_file_name         = '/home/peio/eucalyptus_out.caf'
        long_caf_file      = CafFile(long_caf_file_name)
        
        
if __name__ == "__main__":
    unittest.main()