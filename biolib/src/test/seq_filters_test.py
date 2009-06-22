'Tests for the filtering module'

from biolib.seq_filters import ifiltering_map, create_blast_filter
from biolib.Seqs import Seq

from itertools import ifilter

import unittest

class FilteringIteratorTest(unittest.TestCase):
    '''It test the ifiltering_map function'''
    @staticmethod
    def test_ifiltering_map():
        'It test both combined'
        assert list(ifiltering_map(lambda x: x%2, [1, 2, 3])) == [1, 1]

class BlastFilteringTest(unittest.TestCase):
    'It tests that we can filter out sequence using a blast result'
    @staticmethod
    def test_ifiltering_map():
        'It test both combined'
        #a random sequence
        seq1 = Seq('AACTACGTAGCTATGCTGATGCTAGTCTAGCTAGTCGTAGTCTGATCGTAGTCAGTT')
        #an arabidopsis cdna sequence
        seq2 = Seq('ATGGTGGGTGGCAAGAAGAAAACCAAGATATGTGACAAAGTGTCACATGAGAAGATAG')
        #we keep the sequences with a good hit in the blast
        blast_filter1 = create_blast_filter(expect=1e-10, keep_better_hits=True,
                                            database='tair7_cdna',
                                            program='blastn')
        assert list(ifilter(blast_filter1, [seq1, seq2])) == [seq2]
        #we can also keep the sequences with a bad blast result
        blast_filter2 = create_blast_filter(expect=1e-10,
                                            keep_better_hits=False,
                                            database='tair7_cdna',
                                            program='blastn')
        assert list(ifilter(blast_filter2, [seq1, seq2])) == [seq1]


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testiprscan_parse']
    unittest.main()
