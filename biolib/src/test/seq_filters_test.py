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

from biolib.seq_filters import  create_filter, strip_seq_by_quality
from biolib.biolib_cmd_utils import create_runner
from biolib.seqs import Seq, SeqWithQuality

import unittest
from tempfile import NamedTemporaryFile

#class FilteringIteratorTest(unittest.TestCase):
#    '''It test the ifiltering_map function'''
#    @staticmethod
#    def test_ifiltering_map():
#        'It test both combined'
#        assert list(ifiltering_map(lambda x: x%2, [1, 2, 3])) == [1, 1]

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

        blast_filter = create_filter(aligner_cmd='blast',
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
        fhand = NamedTemporaryFile()
        fhand.write(ADAPTORS)
        fhand.flush()
        parameters = {'target': fhand.name}
        match_filters = [{'kind'           : 'min_scores',
                          'score_key'      : 'score',
                          'min_score_value': 130}]
        result_filters = [{'kind' :'num_matches','value':2}]
        match_filter = create_filter(aligner_cmd='exonerate',
                                      cmd_parameters=parameters,
                                      match_filters=match_filters,
                                      result_filters=result_filters )
        assert match_filter('TACTCTGATCGATCGGATCTAGCATGC') == False
        result_filters = [{'kind' :'num_matches', 'value':1}]
        match_filter = create_filter(aligner_cmd='exonerate',
                                      cmd_parameters=parameters,
                                      match_filters=match_filters,
                                      result_filters=result_filters )
        assert match_filter('TACTCTGATCGATCGGATCTAGCATGC') == True



class StripSeqByQualitytest(unittest.TestCase):
    'test trim_seq_by_quality '
    @staticmethod
    def test_strip_seq_by_quality():
        'test trim_seq_by_quality '
        qual = [20, 20, 20, 60, 60, 60, 60, 60, 20, 20, 20, 20]
        seq  = 'ataataataata'
        new_seq = strip_seq_by_quality(SeqWithQuality(qual=qual, seq=seq),
                                       quality_treshold=40, min_seq_length=2,
                                       min_quality_bases=3)
        assert new_seq.seq == 'ataat'
        qual = [60, 60, 60, 60, 60, 60, 60]
        seq  = 'ataataa'
        new_seq = strip_seq_by_quality(SeqWithQuality(qual=qual, seq=seq),
                                       quality_treshold=40, min_seq_length=2,
                                       min_quality_bases=3)
        assert  new_seq.seq == 'ataataa'
        qual = [60, 60, 60, 60, 60, 60, 0]
        seq  = 'ataataa'
        new_seq = strip_seq_by_quality(SeqWithQuality(qual=qual, seq=seq),
                                       quality_treshold=40, min_seq_length=2,
                                       min_quality_bases=3)
        assert new_seq.seq == 'ataata'



if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testiprscan_parse']
    unittest.main()
