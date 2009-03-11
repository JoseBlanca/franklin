'''
Created on 2009 mar 11

@author: peio
'''
import unittest
from cafparser.cafparser import CafFile

class Test(unittest.TestCase):
    
    def test_caf_parser(self):
        # Primero creamos el objeto
        caf_file = '../doc/Caf_example_file.caf'
        caf_index = CafFile(caf_file)
        assert caf_index
    
#     Show contigs in file
#        for contig in caf_index.contigs():
#            print contig 
#     Show reads in file
        for read in caf_index.reads():
            print read['name'],read
#        return a DNA sec of a item:
#        content = caf_index._return_section(read['BaseQuality'])
#        DNA     = caf_index._get_base_quality(content)
#        name=read['name'] 
#    print name,DNA
#    for read in caf_index.reads():
#        caf_index._get_item_info(read['name'])
#        for contig in caf_index.contigs():
#            caf_index._get_item_info(contig['name'])


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()