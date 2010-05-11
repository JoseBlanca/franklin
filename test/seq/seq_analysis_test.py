'''
Created on 26/11/2009

@author: jose
'''

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

import unittest, os
from franklin.seq.seqs import SeqWithQuality, Seq
from franklin.seq.seq_analysis import (infer_introns_for_cdna,
                                     look_for_similar_sequences,
                                     est2genome_parser,
                                     build_sequence_clusters)
from franklin.utils.misc_utils import DATA_DIR
from Bio import SeqIO
class IntronTest(unittest.TestCase):
    'It test that we can locate introns'

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
        genomic_db = os.path.join(DATA_DIR, 'blast', 'tomato_genome2')
        genomic_seqs_index = SeqIO.index(genomic_db, 'fasta')
        introns = infer_introns_for_cdna(seq1, genomic_db,
                                         genomic_seqs_index=genomic_seqs_index)
        assert introns == [478, 572, 613]

    @staticmethod
    def test_parse_est2genome():
        'It test the parser of the est2genome output'
        output = '''Note Best alignment is between reversed est and forward
genome, but splice sites imply REVERSED GENE
Exon       476  99.8 2270227 2270704 scaffold06070     1   478 SGN-U562593
-Intron    -20   0.0 2270705 2271433 scaffold06070
Exon        94 100.0 2271434 2271527 scaffold06070   479   572 SGN-U562593
+Intron    -20   0.0 2271528 2272627 scaffold06070
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
        assert result['cdna']['exons'][0]['start'] == 1
        assert result['cdna']['exons'][2]['start'] == 573
        assert result['genomic']['exons'][0]['start'] == 2270227
        assert result['genomic']['exons'][2]['start'] == 2272628
        assert result['cdna']['introns'][0] == 478
        assert result['cdna']['introns'][2] == 613
        assert result['genomic']['introns'][0]['start'] == 2270705
        assert result['genomic']['introns'][2]['start'] == 2272669

    @staticmethod
    def test_look_for_similar_seqs():
        'It test if we can look for simialr seqs in a database'
        seq  = 'AATCACCGAGCTCAAGGGTATTCAGGTGAAGAAATTCTTTATTTGGCTCGATGTCGATGAGA'
        seq += 'TCAAGGTCGATCTTCCACCTTCTGATTCAATCTACTTCAAAGTTGGCTTTATCAATAAGAAG'
        seq += 'CTTGATATTGACCAGTTTAAGACTATACATTCTTGTCACGATAATGGTGTCTCTGGCTCTTG'
        seq += 'TGGAGATTCATGGAAGAGTTTTCTCGAGGTAAAGATTAGATCTT'
        seq1 = SeqWithQuality(seq=Seq(seq))
        database = os.path.join(DATA_DIR, 'blast', 'arabidopsis_genes')
        similar_seqs = look_for_similar_sequences(seq1, database=database,
                                                  blast_program='blastn')

        assert similar_seqs[0]['name']          == 'AT5G19860.1'
        assert similar_seqs[0]['query_start']   == 1
        assert similar_seqs[0]['subject_start'] == 323

    @staticmethod
    def test_build_sequence_clusters():
        'It test the build sequences_cluster function'
        blast_fhand = open(os.path.join(DATA_DIR, 'cluster_blast.xml'))
        aligner_config = {'results':{'blast':blast_fhand}}
        clusters = build_sequence_clusters(aligner_config)
        assert len(clusters) == 1
        assert clusters[0] == ['seq1', 'seq2']

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
