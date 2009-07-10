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

from biolib.seq_filters import (create_aligner_filter, create_length_filter,
                                create_adaptor_matches_filter,
                                create_comtaminant_filter)
from biolib.seqs import Seq, SeqWithQuality

import unittest
from tempfile import NamedTemporaryFile
from itertools import ifilter

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
        parameters = {'expect':1e-10, 'database':'tair7_cdna',
                      'program':'blastn'}
        match_filters = [{'kind'           : 'min_scores',
                          'score_key'      : 'expect',
                          'max_score_value': 0.001}]

        blast_filter = create_aligner_filter(aligner_cmd='blast',
                                     cmd_parameters=parameters,
                                     match_filters=match_filters )
        assert  blast_filter(seq1) == False
        assert  blast_filter(seq2) == True

ADAPTORS = '''>adaptor1
AAAAACTAGCTAGTCTACTGATCGTATGTCAAAA
>adaptor2
AAAAATACTCTGATCGATCGGATCTAGCATGCAAA
'''

class TooManyAdaptorsTest(unittest.TestCase):
    'It tests that we can filter out sequences that have too many adaptors'
    @staticmethod
    def test_too_many_adaptors_filtering():
        'It test that we can filter the chimeric sequences using exonerate'
        adapt1 = 'AAAAACTAGCTAGTCTACTGATCGTATGTCAAAA'
        adapt2 = 'AAAAATACTCTGATCGATCGGATCTAGCATGCAA'
        adaptators = [adapt1, adapt2]
        seq2 = 'ATCGACTACGACATCGACTACGATACGATCAGATCAGATCGATCGACTACTA'

        seq = adapt1 + seq2 + adapt2
        seq1    = SeqWithQuality(seq=seq)
        filter_ = create_adaptor_matches_filter(adaptators)
        assert filter_(seq1)

        seq = adapt1 + seq2 + adapt2 +  seq2 + adapt1
        seq1    = SeqWithQuality(seq=seq)
        filter_ = create_adaptor_matches_filter(adaptators)
        assert not filter_(seq1)

class LengthFilterTest(unittest.TestCase):
    'It checks the length filtering'
    @staticmethod
    def test_length_filter():
        'It test the standard sequence length filter'
        seqs = [SeqWithQuality(length=300), SeqWithQuality(length=100)]
        filter_ = create_length_filter(200)
        filtered_seqs = ifilter(filter_, seqs)
        assert len(list(filtered_seqs)[0]) == 300

        seq1 = 'atataAGATAGATA'
        seq2 = 'atgatgatgAAAAA'

        seqs = [SeqWithQuality(seq=seq1), SeqWithQuality(seq=seq2)]
        filter_ = create_length_filter(6, count_masked=False)
        filtered_seqs = ifilter(filter_, seqs)
        assert list(filtered_seqs)[0].seq == seq1

class ContaminantFilterTest(unittest.TestCase):
    'It test contaminant filter'
    @staticmethod
    def test_contaminant_filter():
        'It tests if the sequence has a contaminant'
        seq1 = 'TTGGCAATCGGTTCCTGGATTGGACTTAGACCCCTACGCATCCTCAAATACCAATACAATTGT'
        seq  = SeqWithQuality(seq=seq1)
        filter_by_contaminant = create_comtaminant_filter('tair7_cdna')
        assert not filter_by_contaminant(seq)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testiprscan_parse']
    unittest.main()
