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
from tempfile import NamedTemporaryFile
from Bio import SeqIO

from franklin.seq.seqs import SeqWithQuality, Seq
from franklin.seq.seq_analysis import (infer_introns_for_cdna,
                                     look_for_similar_sequences,
                                     est2genome_parser,
                                     do_transitive_clustering_on_blast,
                                     select_longer_sequence_from_cluster,
                                     do_transitive_clustering,
                                     do_transitive_clustering_all,
                                     get_seqs_without_match)
from franklin.utils.misc_utils import TEST_DATA_DIR
from franklin.seq.alignment_result import BlastParser

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
        tomato_genome = 'tomato_genome2+'

        genomic_db = os.path.join(TEST_DATA_DIR, 'blast', tomato_genome)
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
        database = os.path.join(TEST_DATA_DIR, 'blast', 'arabidopsis_genes+')
        similar_seqs = look_for_similar_sequences(seq1, database=database,
                                                  blast_program='blastn')


        assert similar_seqs[0]['name']          == 'AT5G19860.1'
        assert similar_seqs[0]['query_start']   == 1
        assert similar_seqs[0]['subject_start'] == 323

    @staticmethod
    def test_get_seq_without_match():
        'It gets the seqs without match'
        blast_fhand = open(os.path.join(TEST_DATA_DIR,
                                        'transitive_cluster.blastout.xml'),
                           'rt')

        seqs = [SeqWithQuality(name='seq1', seq=Seq('aa')),
                SeqWithQuality(name='seq2', seq=Seq('aa')),
                SeqWithQuality(name='seq3', seq=Seq('aa')),
                SeqWithQuality(name='seq4', seq=Seq('aa')),
                SeqWithQuality(name='seq5', seq=Seq('aa')),
                SeqWithQuality(name='seq6', seq=Seq('aa'))]
        alignments = BlastParser(fhand=blast_fhand)
        no_matched_seqs = get_seqs_without_match(alignments, seqs)
        assert 'seq6' in  no_matched_seqs
        assert 'seq5' in  no_matched_seqs

class TestTransitiveClustering(unittest.TestCase):
    'It tests the transitive clustering'

    def test_transitive_clustering(self):
        'We do a transitive clustering'

        blast_fhand = open(os.path.join(TEST_DATA_DIR,
                                        'transitive_cluster.blastout.xml'),
                           'rt')
        filter1 = {'kind': 'score_threshold',
                   'score_key': 'similarity',
                   'min_score': 98,
                  }
        filter2 = {'kind': 'min_length',
                   'min_num_residues': 50,
                   'length_in_query': True
                  }
        filters = [filter1, filter2]

        clusters = do_transitive_clustering_on_blast(blast_fhand, filters)
        assert set([u'seq3', u'seq2', u'seq1']) in clusters
        assert set([u'seq4']) in clusters

        # with the secuences
        blast_fhand = open(os.path.join(TEST_DATA_DIR,
                                        'transitive_cluster.blastout.xml'),
                           'rt')

        seqs = [SeqWithQuality(name='seq1', seq=Seq('aa')),
                SeqWithQuality(name='seq2', seq=Seq('aa')),
                SeqWithQuality(name='seq3', seq=Seq('aa')),
                SeqWithQuality(name='seq4', seq=Seq('aa')),
                SeqWithQuality(name='seq5', seq=Seq('aa')),
                SeqWithQuality(name='seq6', seq=Seq('aa'))]

        clusters, no_matched = do_transitive_clustering_all(blast_fhand, seqs,
                                                            filters)
        assert set([u'seq3', u'seq2', u'seq1']) in clusters
        assert set([u'seq4']) in clusters
        assert 'seq5' in no_matched
        assert 'seq6' in no_matched

    @staticmethod
    def test_do_transitive_clustering():
        similar_pairs = [(1,2), (2,3), (4,5), (6,7), (1,6)]
        clusters = do_transitive_clustering(similar_pairs)
        assert len(clusters) == 2
        assert clusters[1] == set([1, 2, 3, 6, 7])
        assert clusters[0] == set([4, 5])


class TestSelectLongerSeq(unittest.TestCase):
    'It tests the selection of the longer sequences'

    def test_select_longer_seq(self):
        'We select one sequence for each group'
        seqs = '>seq1\nAC\n>seq2\nACTG\n>seq3\nACTT\n>seq4\nA\n'
        seqs_fhand = NamedTemporaryFile()
        seqs_fhand.write(seqs)
        seqs_fhand.flush()
        seqs_fname = seqs_fhand.name
        seqs_index = SeqIO.index(seqs_fname, 'fasta')
        clusters = [('seq1', 'seq2'), ('seq3', 'seq4')]
        long_seqs = select_longer_sequence_from_cluster(seqs_index, clusters)
        assert [seq.name for seq in long_seqs] == ['seq2', 'seq3']


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
