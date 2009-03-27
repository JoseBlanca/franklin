'''
Created on 2009 mar 25

@author: peio
'''
import unittest
from biolib.SeqVariation import SeqVariation, SeqVarpConf, seq_var_in_alignment
from biolib.contig import Contig, Location
#from test.test_utils import Seq, SeqWithQuality, SeqRecord, Seqmut
from biolib.SeqRecord import SeqRecord
class SeqVariationTest(unittest.TestCase):
    ''' Here we will check if the SeqVariation module works as it should do
    '''
    @staticmethod
    def test_seq_variation_init():
        '''It checks if the '''
       
        contig = Contig(SeqRecord('aa', name='seq1'),
                        SeqRecord('at',name='seq2'))
        for seq_var in seq_var_in_alignment(contig):
            print type(seq_var)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_SeqVariation_init']
    unittest.main()