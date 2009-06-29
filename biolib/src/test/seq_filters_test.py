'Tests for the filtering module'

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of biolib.
# biolib is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# biolib is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with biolib. If not, see <http://www.gnu.org/licenses/>.

from biolib.seq_filters import (ifiltering_map, create_blast_filter,
                                create_ssaha_filter)
from biolib.seqs import Seq

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
    def test_blast_filtering():
        'It test that we can filter using blast'
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

adaptors = '''>adaptor1
AAAAACTAGCTAGTCTACTGATCGTATGTCAAAA
>adaptor2
AAAAATACTCTGATCGATCGGATCTAGCATGCAAA
'''

class TooManyAdaptorsTest(unittest.TestCase):
    'It tests that we can filter out sequences that have too many adaptors'
    @staticmethod
    def test_ssaha2_filtering():
        'It test that we can filter the chimeric sequences using ssaha2'
        ssaha2_filter = create_ssaha_filter(options='adaptors',
                                            subject=adaptors)
        ssaha2_filter = create_ssaha_filter(subject=adaptors)
        ssaha2_filter(Seq('TACTCTGATCGATCGGATCTAGCATGC'))


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testiprscan_parse']
    unittest.main()
