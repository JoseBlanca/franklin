'''
Created on 2009 mar 11

@author: peio
'''
import unittest
from biolib.cafparser import CafFile

class Test(unittest.TestCase):
    ''' It tests '''
    def setUp(self):
        self._file2test = '/home/peio/eucalyptus_out.caf'
        
    def test_caf_parser(self):
        ''' It tests if we can create and caf file giving a file name'''
#        caf_file_name = '../doc/Caf_example_file.caf'
        caf_file = CafFile(self._file2test)
        assert caf_file

    def test_reads(self):
        ''' we check if we can take the reads from the caf file'''
#        caf_file_name = '../doc/Caf_example_file.caf'
        caf_file = CafFile(self._file2test)
        
        for read in caf_file.reads():
            read_object = caf_file.read_record2read(read['name'])
            assert read_object
    
    def test_contigs(self):
        ''' It checks if the contig method returns contigs'''
#        caf_file_name = '../doc/Caf_example_file.caf'
        caf_file = CafFile(self._file2test)
        for contig in caf_file.contigs():
            contig_object = caf_file.contig_record2contig(contig['name'])
            print contig_object 
            

if __name__ == "__main__":
    unittest.main()