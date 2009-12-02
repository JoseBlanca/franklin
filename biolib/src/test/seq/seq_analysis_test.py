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
                                     look_for_similar_sequences,
                                     est2genome_parser)
from biolib.utils.misc_utils import DATA_DIR
from biolib.utils.seqio_utils import FileSequenceIndex

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
        introns = infer_introns_for_cdna(seq, genomic_db=genomic_db,
                                         method='blast')
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
        introns = infer_introns_for_cdna(seq, genomic_db=genomic_db,
                                         method='blast')
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
    def test_infer_introns_est2genome_method():
        'It tests the est2genome method of infering introns'
        seq  = 'GAAAAGATGTGATTGGTGAAATAAGTTTGCCTCAATTCTCTTGTGCCGAAGTTCCAAAGAAGC'
        seq += 'AGTTGGTGAATGAGCAGCCAGTACCCGAAAAATCGAGCAAAGATTTTGTGATGTATGTTGGAG'
        seq += 'GTCTAGCATGGGGGATGGACTGGTGTCCCCAAGCTCATGAAAATAGGGATGCTCCTATGAAAA'
        seq += 'GTGAGTTTGTCGCAATTGCTCCTCATCCTCCTGATTCATCATATCACAAGACTGATGCCTCAC'
        seq += 'TTACAGGCAGAGGTGTAATTCAGATATGGTGCCTGCCAGATCTCATTCAAAAAGATATAATTG'
        seq += 'TGAAAGAAGATTATTTTGCTCAGGTTAACAAAAAACCGTATAGAAATTTGACAAGAAGTGAAG'
        seq += 'CAGGTACGGGAGAAGTATCTGGACCTCAAAAACCAAGAGGAAGACCAAAAAAGAACCCTGGTA'
        seq += 'AAGCAGTCCAGGCAAAAGCATCTAGACCACAAAATCCAAGAGGAAGACCGAGAAAGAAGCCTG'
        seq += 'TTACTGAATCTTTAGGTGATAGAGATAGTGAAGACCACAGTTTACAACCTCTTGCTATAGAGT'
        seq += 'GGTCGCTGCAATCAACAGAACTTTCTGTAGATTTGTCTTGTGGAAATATGAATAAAGCCCAAG'
        seq += 'TAGATATTGCGCTGAGTCAAGAAAGATGTATTAATGCGGCAT'
        seq1 = SeqWithQuality(seq = Seq(seq))
        genomic_db = os.path.join(DATA_DIR, 'blast', 'tomato_genome')
        genomic_seqs_index = FileSequenceIndex(open(genomic_db), 'fasta')
        introns = infer_introns_for_cdna(seq1, genomic_db,
                                         genomic_seqs_index=genomic_seqs_index)
        assert introns == ['478', '572', '613']

    @staticmethod
    def test_parse_est2genome():
        'It test the parser of the est2genome output'
        output = '''Note Best alignment is between reversed est and forward
genome, but splice sites imply REVERSED GENE
Exon       476  99.8 2270227 2270704 scaffold06070     1   478 SGN-U562593
-Intron    -20   0.0 2270705 2271433 scaffold06070
Exon        94 100.0 2271434 2271527 scaffold06070   479   572 SGN-U562593
-Intron    -20   0.0 2271528 2272627 scaffold06070
Exon        41 100.0 2272628 2272668 scaffold06070   573   613 SGN-U562593
-Intron    -20   0.0 2272669 2272767 scaffold06070
Exon        57  98.3 2272768 2272826 scaffold06070   614   672 SGN-U562593

Span       608  99.7 2270227 2272826 scaffold06070     1   672 SGN-U562593

Segment    476  99.8 2270227 2270704 scaffold06070     1   478 SGN-U562593
Segment     94 100.0 2271434 2271527 scaffold06070   479   572 SGN-U562593
Segment     41 100.0 2272628 2272668 scaffold06070   573   613 SGN-U562593
Segment     57  98.3 2272768 2272826 scaffold06070   614   672 SGN-U562593'''
        result = est2genome_parser(output)
        assert len(result['cdna']['introns']) == 3
        assert result['cdna']['exons'][0]['start'] == '1'
        assert result['cdna']['exons'][2]['start'] == '573'
        assert result['genomic']['exons'][0]['start'] == '2270227'
        assert result['genomic']['exons'][2]['start'] == '2272628'
        assert result['cdna']['introns'][0] == '478'
        assert result['cdna']['introns'][2] == '613'
        assert result['genomic']['introns'][0]['start'] == '2270705'
        assert result['genomic']['introns'][2]['start'] == '2272669'



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
        assert similar_seqs[0]['name'] == 'AT5G19860.1'
        assert similar_seqs[0]['query_start'] == 1
        assert similar_seqs[0]['subject_start'] == 323


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()