'''It tests the parser for the ace alignment files'''

import unittest
from biolib.aceparser import AceParser
import biolib
import os.path

DATA_DIR = os.path.join(os.path.split(biolib.__path__[0])[0], 'data')

class AceTest(unittest.TestCase):
    ''' It tests the ace alignment parser'''
    
    def test_contig(self):
        '''It tests that we can get a read by its name.'''
        filen = os.path.join(DATA_DIR, 'example3.ace')
        
        parser = AceParser(filen)
        #if we ask for a wrong contig we get an error
        try:
            parser.contig('not_in_file')
            self.fail('KeyError expected')
            #pylint: disable-msg=W0704
        except ValueError:
            pass
        #now for a real contig
        contig = AceParser(filen).contig('eucalyptus_lrc1')
        print contig
 
if __name__ == "__main__":
    unittest.main()
