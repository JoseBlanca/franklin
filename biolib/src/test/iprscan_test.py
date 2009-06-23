'''
Created on 2009 mai 21

@author: peio
'''
import unittest, os, biolib
from biolib.iprscan import xml_iprscan_parser_iter

DATA_DIR = os.path.join(os.path.split(biolib.__path__[0])[0], 'data')

class IprsacnTest(unittest.TestCase):
    '''It test all functions in iprscan parser '''
    
    @staticmethod
    def test_iprscan_parse():
        '''Tests iprscan parser '''
        fname = os.path.join(DATA_DIR, 'iprscan.xml')
        fhand = open(fname, 'r')
        protein_ids = ["PRRB_MYCTU", "Q9RHD9_PSEAE", "RS16_ECOLI"]
        cont = 0
        for protein_result in xml_iprscan_parser_iter(fhand):
            assert  protein_result['id'] in protein_ids
            cont += 1
        assert cont == 3
        
        
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testiprscan_parse']
    unittest.main()