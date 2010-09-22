'Tests for the filtering module'

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of franklin.
# franklin is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# franklin is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with franklin. If not, see <http://www.gnu.org/licenses/>.

from franklin.seq.seq_filters import (create_aligner_filter, create_length_filter,
                                create_adaptor_matches_filter,
                                create_comtaminant_filter,
                                create_similar_seqs_filter,
                                create_solid_quality_filter)
from franklin.seq.seqs import Seq, SeqWithQuality
from franklin.utils.misc_utils import DATA_DIR
from Bio.Seq import UnknownSeq

import unittest, os
from itertools import ifilter
from franklin.utils.cmd_utils import BLAST_TOOL

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
        
        if BLAST_TOOL == 'blast+':
            arabidopsis_genes = 'arabidopsis_genes+' 
        else:
            arabidopsis_genes = 'arabidopsis_genes'
        parameters = {'expect':1e-10, 'database':arabidopsis_genes,
                      'program':'blastn'}
        match_filters = [{'kind'           : 'min_scores',
                          'score_key'      : 'expect',
                          'max_score_value': 0.001}]

        blastpath = os.path.join(DATA_DIR, 'blast')
        blast_filter = create_aligner_filter(aligner_cmd=BLAST_TOOL,
                                     cmd_parameters=parameters,
                                     match_filters=match_filters,
                                     environment={'BLASTDB':blastpath} )
        assert  blast_filter(seq1) == False
        assert  blast_filter(seq2) == False

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
        adaptors = [adapt1, adapt2]
        seq2 = 'ATCGACTACGACATCGACTACGATACGATCAGATCAGATCGATCGACTACTA'

        seq = adapt1 + seq2 + adapt2
        seq1    = SeqWithQuality(seq=Seq(seq))
        filter_ = create_adaptor_matches_filter(adaptors)
        assert filter_(seq1)

        seq = adapt1 + seq2 + adapt2 +  seq2 + adapt1
        seq1    = SeqWithQuality(seq=Seq(seq))
        filter_ = create_adaptor_matches_filter(adaptors)
        assert not filter_(seq1)

class LengthFilterTest(unittest.TestCase):
    'It checks the length filtering'
    @staticmethod
    def test_length_filter():
        'It test the standard sequence length filter'
        seqs = [SeqWithQuality(UnknownSeq(length=300)),
                SeqWithQuality(UnknownSeq(length=100))]
        filter_ = create_length_filter(200)
        filtered_seqs = ifilter(filter_, seqs)
        assert len(list(filtered_seqs)[0]) == 300

        seq1 = Seq('atataAGATAGATA')
        seq2 = Seq('atgatgatgAAAAA')

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
        seq  = SeqWithQuality(seq=Seq(seq1))
        blastpath = os.path.join(DATA_DIR, 'blast')
        if BLAST_TOOL == 'blast+':
            arabidopsis_genes = 'arabidopsis_genes+' 
        else:
            arabidopsis_genes = 'arabidopsis_genes'
        filter_by_contaminant = create_comtaminant_filter(arabidopsis_genes,
                                              environment={'BLASTDB':blastpath})
        assert not filter_by_contaminant(seq)

class SimilarSeqTest(unittest.TestCase):
    'It test if there are similar seqs in the database'
    @staticmethod
    def test_similar_seqs():
        'It shoul return True if there are similar seqs'
        seq  = 'AATCACCGAGCTCAAGGGTATTCAGGTGAAGAAATTCTTTATTTGGCTCGATGTCGATGAGA'
        seq += 'TCAAGGTCGATCTTCCACCTTCTGATTCAATCTACTTCAAAGTTGGCTTTATCAATAAGAAG'
        seq += 'CTTGATATTGACCAGTTTAAGACTATACATTCTTGTCACGATAATGGTGTCTCTGGCTCTTG'
        seq += 'TGGAGATTCATGGAAGAGTTTTCTCGAGGTAAAGATTAGATCTT'
        seq1 = SeqWithQuality(seq=Seq(seq))
        if BLAST_TOOL == 'blast+':
            arabidopsis_genes = 'arabidopsis_genes+' 
        else:
            arabidopsis_genes = 'arabidopsis_genes'
        db = os.path.join(DATA_DIR, 'blast', arabidopsis_genes)
        filter_ = create_similar_seqs_filter(db=db, blast_program='blastn')
        assert filter_(seq1)

        filter_ = create_similar_seqs_filter(db=db, blast_program='blastn',
                                             inverse=True)
        assert not filter_(seq1)

        filter_ = create_similar_seqs_filter(db=db, blast_program='blastn',
                                             min_sim_seqs=2)
        assert not filter_(seq1)

class SolidFilters(unittest.TestCase):
    'It tests solid filters'
    @staticmethod
    def solid_quality_filter():
        'It test solid quality filters'
        quality = [30, 23, 43, 12, 25, 23, 30, 12, 0]
        seq = 'A' * len(quality)
        sequence =  SeqWithQuality(seq=Seq(seq), qual=quality)
        
        filter_ = create_solid_quality_filter(length=5, threshold=25)
        assert filter_(sequence)

        filter_ = create_solid_quality_filter(length=5, threshold=35)
        assert not filter_(sequence)
        
        quality = [30, 23, 43, 12, 25, 23, 30, 12, 0, 0]
        seq = 'A' * len(quality)
        sequence =  SeqWithQuality(seq=Seq(seq), qual=quality)
        filter_ = create_solid_quality_filter(length=5, threshold=25, 
                                              call_missing=True)
        assert not filter_(sequence)

if __name__ == "__main__":
#    import sys;sys.argv = ['', 'SolidFilters.solid_quality_filter']
    unittest.main()
