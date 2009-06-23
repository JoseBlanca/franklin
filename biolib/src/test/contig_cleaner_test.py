'''
Created on 2009 api 28

@author: peio
'''
import unittest
from biolib.contig import Contig
from biolib.locatable_sequence import locate_sequence
from biolib.contig_cleaner import contig_strip, water_alignment_strip
from Bio.Seq import Seq
class StripTest(unittest.TestCase):
    "It test the contig strip functions, in all methods"
    
    @staticmethod
    def test_contig_strip():
        '''It tests the contig strip function'''
        contig = Contig()
        seq = 'tTCTTAAGGTTATGCGTACGTGCAGTag'
        mask = (1, 26)
        contig.append_to_location(seq, mask=mask)
        seq = 'tTCTTAAGGTTATGCGTACGTGCAGTag'
        mask = (1, 26)
        contig.append_to_location(seq, mask=mask)
        contig1 = contig_strip(contig, 3)
        assert str(contig1[0]).strip() == 'TAAGGTTATGCGTACGTGCA'
        contig2 = contig_strip(contig, 5)
        assert str(contig2[0]).strip() == 'AGGTTATGCGTACGTG'
        
        contig3 = Contig()
        seq = Seq('tTCTTAAGGTTAag')
        mask = (2, 12)
        contig3.append_to_location(seq, mask=mask, forward=False, strand=-1)
        contig = contig_strip(contig3, 2)
        assert str(contig[0]).strip() == 'AACCTTA'
    
    @staticmethod
    def test_contig_water_strip():
        '''It test the contig water strip function '''
        contig = Contig(consensus=locate_sequence(sequence='AATTCCGG', \
                                                  location=1))
        seq = 'tAATTCCGGt'
        mask = (1, 8)
        contig.append_to_location(seq, mask=mask)
        
        contig = water_alignment_strip(contig)
        assert str(contig[0]).strip() == 'AATTCCGG'
        
              
        contig = Contig(consensus=locate_sequence(sequence='AATTCCGG', \
                                                  location=1))
        seq = 'taaTTCCggt'
        mask = (3, 6)
        contig.append_to_location(seq, mask=mask)
        
        contig = water_alignment_strip(contig)
        assert str(contig[0]).strip() == 'TTCC'


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()