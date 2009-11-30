'''
Created on 26/11/2009

@author: jose
'''

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

import unittest, os
from biolib.seq.seqs import SeqWithQuality, Seq
from biolib.seq.seq_analysis import (infer_introns_for_cdna,
                                     _infer_introns_from_matches,
                                     look_for_similar_sequences)
from biolib.utils.misc_utils import DATA_DIR

class IntronTest(unittest.TestCase):
    'It test that we can locate introns'
    @staticmethod
    def test_basic_introns_for_cdna():
        'It tests that we can locate introns'
        seq1 =  'ATGATAATTATGAAAAATAAAATAAAATTTAATTATATAATTCATTTCATCTAATCGTACAA'
        seq1 += 'GCTAGATATTACTATATCAACAACTTTGTGTATAAAAAGGGCAAGAAATTAAGCATTATCGT'
        seq1 += 'GTGAGCCACTTTTTCTATATCTAGAGATAGAAGGTTTAAAATCATGTCTCTAATTGGAAAGC'
        seq1 += 'TTGTGAGTGAATTAGAGATCAATGCAGCTGCTGAGAAATTTTACGAAATATTCAAAGATCAA'
        seq1 += 'TGTTTTCAGGTTCCCAATATAACCCCCAGATGCATTCAACAAGTTGAAATTCATGGTACTAA'
        seq1 += 'TTGGGATGGCCATGGACATGGCTCTATCAAGTCTTGGTATTACACTATTGATGGCAAGGCAG'
        seq1 += 'AAGTTTTTAAGGAACGGGTCGAGTTTCACGATGATAAATTGTTGATAGTCTTGGATGGAGTG'
        seq1 += 'GGAGGAGATGTGTTCAAAAATTATAAAAGCTTTAAACCAGCTTACCAATTTGTACCTAAGGA'
        seq1 += 'TCGTAACCATTGCCAGGCAATTCTGAGTATAGAGTATGAGAAACTTCATCATGGGTCTCCTG'
        seq1 += 'ATCCTCATAAGTATATTGACCTCATGATTGGTATCACTAACGACATTGGATCTCACATTAAA'
        seq1 += 'TAAGTATTTAATGTCTGTCACATTCTCAAGTGTGGCTTGTTAATTTGTTGTGGGAAAGTTAT'
        seq1 += 'ATTTTATTTTGAAGTCATTTTCGTGTGGTTGATTATGTATCTTTGCTATTTTGCTTTTATAT'
        seq1 += 'TTCAATAAGTTATATGGTTTATATAATATTACAAAGTAAATAAAATCCAAGGATCATCCCTT'
        seq1 += 'GTTTATGTTTCGTTATTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
        seq1 += 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
        seq1 += 'AAAAAAAAAAAAAAAAAAAGGCCCCCCCCCCCCCCCAAAAAAATTAAAAAACCCCCCCCCCC'
        seq1 += 'CGGGGGGGGCCC'
        seq = SeqWithQuality(name='seq', seq=Seq(seq1))
        genomic_db = os.path.join(DATA_DIR, 'blast', 'arabidopsis_genes')
        introns = infer_introns_for_cdna(seq, genomic_db=genomic_db)
        assert introns == [359]

        seq2  = 'AAACCAACTTTCTCCCCATTTTCTTCCTCAAACCTCCATCAATGGCTTCCTTCTCCAGAATC'
        seq2 += 'CTCTCCCCATTTTCACTATTTCTTCTGATTCTTGTCATCTCCACTCAAACCCACCTCTCCTT'
        seq2 += 'TTCAGCAAGGGATCTTCTCCTCAAGTCATCTGATATCCACGATCTTCTTCCCCTTTACGGTT'
        seq2 += 'TTCCAGTCGGTCTCTTACCCAGCAATGTCAAGTCCTACACTCTCTCAGACGATGGTAGCTTC'
        seq2 += 'GTAATCGAACTCGATAGCGCTTGCTATGTCCAGTTCGCTGATCTGGTCTATTACGGCAAGAC'
        seq2 += 'GATCAAGGGGAAATTGAGCTATGGGTCATTGAGCGATGTTTCTGGGATTCAAGTCAAGAAGT'
        seq2 += 'TGTTCGCCTGGCTTCCTATTACTGGAATGAGGGTTACTTCAGACTCTAAATCCATCGAGTTT'
        seq2 += 'CAGGTTGGGTTCTTGTCTGAGGCTTTGCCGTTCAGCATGTTTGAGTCCATTCCTACATGCAG'
        seq2 += 'AAAGAAAGCTTGCCTAGAAGGGAAAACAGAGGCAGTGTGAGGTGGAAATAATAGCTTTCCAA'
        seq2 += 'AACGCTTATCCTTTTCATTGGGTGAGAGAAGCATGTTGGTCTTTGCAAGAAGAATAATGTAA'
        seq2 += 'TCTTTGTTTTTATGTCATGAACCTACGGTGTCCATTTTTAATCTTTTTCTTACATGTTCATC'
        seq2 += 'TATATTTATATC-ATATCATAAATATTCTCACATGTTTACCTAATGTTTTCTTTCAATAATA'
        seq2 += 'TTATCTTTTTACGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
        seq = SeqWithQuality(name='seq', seq=Seq(seq2))
        introns = infer_introns_for_cdna(seq, genomic_db=genomic_db)
        assert introns == []

    @staticmethod
    def test_basic_introns_for_cdna_extra():
        'It tests that we can locate introns'
        #no intron, two hsps that follow
        # \
        #  \
        #
        #    \
        #     \
        #      \
        hsp1 = {'query_start': 1,
                'query_end': 10,
                'query_strand': 1,
                'subject_start': 1,
                'subject_end': 10,
                'subject_strand': 1}
        hsp2 = {'query_start': 11,
                'query_end': 20,
                'query_strand': 1,
                'subject_start': 11,
                'subject_end': 20,
                'subject_strand': 1}
        match = {'match_parts':[hsp1, hsp2]}
        alignment = {'matches':[match]}
        alignments = iter([alignment])
        introns = _infer_introns_from_matches(alignments)
        assert introns == []
        # \
        #  \
        #
        #       \
        #        \
        #         \
        hsp1 = {'query_start': 1,
                'query_end': 10,
                'query_strand': 1,
                'subject_start': 1,
                'subject_end': 10,
                'subject_strand': 1}
        hsp2 = {'query_start': 11,
                'query_end': 20,
                'query_strand': 1,
                'subject_start': 21,
                'subject_end': 30,
                'subject_strand': 1}
        match = {'match_parts':[hsp1, hsp2]}
        alignment = {'matches':[match]}
        alignments = iter([alignment])
        introns = _infer_introns_from_matches(alignments)
        assert introns == [10]
        # \
        #  \
        #
        # \
        #  \
        #   \
        hsp1 = {'query_start': 10,
                'query_end': 100,
                'query_strand': 1,
                'subject_start': 10,
                'subject_end': 100,
                'subject_strand': 1}
        hsp2 = {'query_start': 110,
                'query_end': 200,
                'query_strand': 1,
                'subject_start': 10,
                'subject_end': 100,
                'subject_strand': 1}
        match = {'match_parts':[hsp1, hsp2]}
        alignment = {'matches':[match]}
        alignments = iter([alignment])
        introns = _infer_introns_from_matches(alignments)
        assert introns == []

        #reverse
                #no intron, two hsps that follow
        # \
        #  \
        #
        #    \
        #     \
        #      \
        hsp1 = {'query_start': 1,
                'query_end': 10,
                'query_strand': 1,
                'subject_start': 1,
                'subject_end': 10,
                'subject_strand': 1}
        hsp2 = {'query_start': 11,
                'query_end': 20,
                'query_strand': 1,
                'subject_start': 11,
                'subject_end': 20,
                'subject_strand': -1}
        match = {'match_parts':[hsp1, hsp2]}
        alignment = {'matches':[match]}
        alignments = iter([alignment])
        introns = _infer_introns_from_matches(alignments)
        assert introns == []
        # \
        #  \
        #
        #       \
        #        \
        #         \
        hsp1 = {'query_start': 1,
                'query_end': 10,
                'query_strand': -1,
                'subject_start': 1,
                'subject_end': 10,
                'subject_strand': 1}
        hsp2 = {'query_start': 11,
                'query_end': 20,
                'query_strand': -1,
                'subject_start': 21,
                'subject_end': 30,
                'subject_strand': 1}
        match = {'match_parts':[hsp1, hsp2]}
        alignment = {'matches':[match]}
        alignments = iter([alignment])
        introns = _infer_introns_from_matches(alignments)
        assert introns == [10]
        # \
        #  \
        #
        # \
        #  \
        #   \
        hsp1 = {'query_start': 10,
                'query_end': 100,
                'query_strand': -1,
                'subject_start': 10,
                'subject_end': 100,
                'subject_strand': 1}
        hsp2 = {'query_start': 110,
                'query_end': 200,
                'query_strand': -1,
                'subject_start': 10,
                'subject_end': 100,
                'subject_strand': -1}
        match = {'match_parts':[hsp1, hsp2]}
        alignment = {'matches':[match]}
        alignments = iter([alignment])
        introns = _infer_introns_from_matches(alignments)
        assert introns == []

    @staticmethod
    def test_look_for_similar_seqs():
        'It test if we can look for simialr seqs in a database'
        seq  = 'AATCACCGAGCTCAAGGGTATTCAGGTGAAGAAATTCTTTATTTGGCTCGATGTCGATGAGA'
        seq += 'TCAAGGTCGATCTTCCACCTTCTGATTCAATCTACTTCAAAGTTGGCTTTATCAATAAGAAG'
        seq += 'CTTGATATTGACCAGTTTAAGACTATACATTCTTGTCACGATAATGGTGTCTCTGGCTCTTG'
        seq += 'TGGAGATTCATGGAAGAGTTTTCTCGAGGTAAAGATTAGATCTT'
        seq1 = SeqWithQuality(seq=Seq(seq))
        db = os.path.join(DATA_DIR, 'blast', 'arabidopsis_genes')
        similar_seqs = look_for_similar_sequences(seq1, db=db,
                                                  blast_program='blastn')
        assert similar_seqs == ['AT5G19860.1']


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()