'''
Created on 15/01/2010

@author: peio
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

from franklin.seq.seq_annotation import (create_microsatellite_annotator,
                                         create_ortholog_annotator,
                                         create_description_annotator,
                                         create_orf_annotator,
                                         create_go_annotator,
                                         create_cdna_intron_annotator)
from franklin.seq.seqs import SeqWithQuality, Seq
from franklin.utils.misc_utils import DATA_DIR
import unittest, tempfile, os

class AnnotationTests(unittest.TestCase):
    'Annotations tests'
    @staticmethod
    def test_orthologs_annotator():
        'It test the ortholog annotator'
        blast_fhand = open(os.path.join(DATA_DIR, 'melon_tair.xml'))
        reverse_blast_fhand = open(os.path.join(DATA_DIR, 'tair_melon.xml'))
        blast = {'blast':blast_fhand}
        reverse_blast = {'blast':reverse_blast_fhand}
        ortho_annotator = create_ortholog_annotator(blast, reverse_blast,
                                                    species='arabidopsis')
        sequence = SeqWithQuality(seq=Seq('aaa'), name='melon1')
        sequence = ortho_annotator(sequence)
        assert sequence.annotations['arabidopsis-orthologs'] == ['tair1']

        sequence = SeqWithQuality(seq=Seq('aaa'), name='melon2')
        sequence = ortho_annotator(sequence)
        assert sequence.annotations['arabidopsis-orthologs'] == ['tair2']

    @staticmethod
    def test_get_description_with_funct():
        'It tests if we can get description for seqs in blasts. with mod funct'
        # test with a modifier function
        blast_fhand = open(os.path.join(DATA_DIR, 'blast2.xml'))
        blast = {'blast':blast_fhand,
                 'modifier':lambda(x):x.split('|')[2]}
        descrip_annotator = create_description_annotator([blast])
        sequence = SeqWithQuality(seq=Seq('aaa'), name='CUTC021854')
        sequence = descrip_annotator(sequence)
        assert sequence.description == 'Similar to ankyrin repeat family protein'
        sequence = SeqWithQuality(seq=Seq('aaa'), name='CUTC021853')
        sequence = descrip_annotator(sequence)
        assert sequence.description == 'Similar to DNA-binding protein-related'

    @staticmethod
    def test_microsatelite_annotator():
        'It test the srrs annotator'
        seq = 'atgatgatgatgatgatgatgatgatgatggcgcgcgcgcgcgcgcgcgcgcgcgcg'
        ssr_annot = create_microsatellite_annotator()
        seq1 = SeqWithQuality(seq=Seq(seq))
        ssr_annot(seq1)
        assert seq1.features[0].qualifiers['score'] == 27

    @staticmethod
    def test_orf_annotator():
        'It tests that we can annotate orfs'
        seq = 'CTACTTACTAGCTTTAGTAAATCCTTCTAACCCTCGGTAAAAAAAAAAAAGAGGCATCAAATG'
        seq += 'GCTTCATCCATTCTCTCATCCGCCGNTGTGGCCTTTGNCAACAGGGCTTCCCCTGCTCAAGCT'
        seq += 'AGCATGGGGGCACCATTCACTGGCCTAAAATCCGCCGCTGCTTTCCCNGTNACTCGCANGACC'
        seq += 'AACGACATCACCACTTTGGTTAGCAATGGGGGAAGAGTTCAGGGCNTGAAGGTGTGCCCACCA'
        seq += 'CTTGGATTGAAGAAGTTCGAGACTCTTTCTTACCTTCCTGATATGAGTAACGAGCAATTGGGA'
        seq += 'AAGGAAGTTGACTACCTTCTCAGGAAGGGATGGATTCCCTGCATTGAATTCGACATTCACAGT'
        seq += 'GGATTCGTTTACCGTGAGACCCACAGGTCACCAGGATACTTCGATGGACGCTACTGGACCATG'
        seq += 'TGGAAGCTGCCCATGTTTGGCTGCACCGAT'
        seq1 = SeqWithQuality(seq=Seq(seq))
        matrix_fpath = os.path.join(DATA_DIR, 'At.smat')
        annotator = create_orf_annotator(parameters={'matrix':matrix_fpath})
        annotator(seq1)
        assert len(seq1.features) == 1
        assert seq1.features[0].type == 'orf'


    @staticmethod
    def test_intron_annotator():
        'We can annotate introns in cdnas comparing with genomic'
        seq = 'GAAAAGATGTGATTGGTGAAATAAGTTTGCCTCAATTCTCTTGTGCCGAAGTTCCAAAGAAGC'
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
        seq = SeqWithQuality(name='seq', seq=Seq(seq))
        blast_db_path = os.path.join(DATA_DIR, 'blast')
        genomic_db = os.path.join(blast_db_path, 'tomato_genome2')
        intron_annotator = create_cdna_intron_annotator(genomic_db=genomic_db,
                                            genomic_seqs_fhand=open(genomic_db))
        intron_annotator(seq)
        intron_feat = seq.features[0]
        assert intron_feat.location.start.position == 478
        assert intron_feat.type == 'intron'

    @staticmethod
    def test_go_annotator():
        'It test the go annotator'
        blast = open(os.path.join(DATA_DIR, 'blastResult.xml'))
        fhand, annot_fpath = tempfile.mkstemp()
        os.close(fhand)
        go_annotator = create_go_annotator(blast)
        seq = SeqWithQuality(name='seq1', seq=Seq('aaaa'))

        go_annotator(seq)
        assert 'GO:0009853' in seq.annotations['GOs']

        os.remove(annot_fpath)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testiprscan_parse']
    unittest.main()
