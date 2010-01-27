'''
Created on 15/01/2010

@author: peio
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

import os, unittest
from biolib.utils.misc_utils import DATA_DIR
from biolib.seq.seqs import SeqWithQuality
from biolib.seq.seq_annotation import (create_microsatelite_annotator,
                                       create_ortholog_annotator,
                                       create_description_annotator,
                                       create_orf_annotator)

class AnnotationTests(unittest.TestCase):
    'Annotations tests'
    @staticmethod
    def test_orthologs_annotator():
        'It test the ortholog annotator'
        blast_fhand  = open(os.path.join(DATA_DIR, 'melon_tair.xml'))
        reverse_blast_fhand = open(os.path.join(DATA_DIR, 'tair_melon.xml'))
        blast = {'results':{'blast':blast_fhand}}
        reverse_blast = {'results':{'blast':reverse_blast_fhand}}
        ortho_annotator = create_ortholog_annotator(blast, reverse_blast,
                                                    species='arabidopsis')
        sequence = SeqWithQuality(seq='aaa', name='melon1')
        sequence = ortho_annotator(sequence)
        assert sequence.annotations['arabidopsis-orthologs'] == ['tair1']

        sequence = SeqWithQuality(seq='aaa', name='melon2')
        sequence = ortho_annotator(sequence)
        assert sequence.annotations['arabidopsis-orthologs'] == ['tair2']

    @staticmethod
    def test_get_description_with_funct():
        'It tests if we can get description for seqs in blasts. with mod funct'
        # test with a modifier function
        blast_fhand = open(os.path.join(DATA_DIR, 'blast2.xml'))
        blast = {'results':{'blast':blast_fhand}, 'modifier':lambda(x):x.split('|')[2]}
        descrip_annotator = create_description_annotator([blast])
        sequence = SeqWithQuality(seq='aaa', name='CUTC021854')
        sequence = descrip_annotator(sequence)
        assert sequence.annotations['description'] == \
                                                'ankyrin repeat family protein'
        sequence = SeqWithQuality(seq='aaa', name='CUTC021853')
        sequence = descrip_annotator(sequence)
        assert sequence.annotations['description'] == \
                                                'DNA-binding protein-related'

    @staticmethod
    def test_microsatelite_annotator():
        'It test the srrs annotator'
        seq = 'atgatgatgatgatgatgatgatgatgatggcgcgcgcgcgcgcgcgcgcgcgcgcg'
        ssr_annot = create_microsatelite_annotator()
        seq1 = SeqWithQuality(seq=seq)
        ssr_annot(seq1)
        assert seq1.features[0].qualifiers['score'] == 27

    @staticmethod
    def test_orf_annotator():
        'It tests that we can annotate orfs'
        seq  = 'CTACTTACTAGCTTTAGTAAATCCTTCTAACCCTCGGTAAAAAAAAAAAAGAGGCATCAAATG'
        seq += 'GCTTCATCCATTCTCTCATCCGCCGNTGTGGCCTTTGNCAACAGGGCTTCCCCTGCTCAAGCT'
        seq += 'AGCATGGGGGCACCATTCACTGGCCTAAAATCCGCCGCTGCTTTCCCNGTNACTCGCANGACC'
        seq += 'AACGACATCACCACTTTGGTTAGCAATGGGGGAAGAGTTCAGGGCNTGAAGGTGTGCCCACCA'
        seq += 'CTTGGATTGAAGAAGTTCGAGACTCTTTCTTACCTTCCTGATATGAGTAACGAGCAATTGGGA'
        seq += 'AAGGAAGTTGACTACCTTCTCAGGAAGGGATGGATTCCCTGCATTGAATTCGACATTCACAGT'
        seq += 'GGATTCGTTTACCGTGAGACCCACAGGTCACCAGGATACTTCGATGGACGCTACTGGACCATG'
        seq += 'TGGAAGCTGCCCATGTTTGGCTGCACCGAT'
        seq1 = SeqWithQuality(seq=seq)
        matrix_fpath = os.path.join(DATA_DIR, 'At.smat')
        annotator = create_orf_annotator(parameters={'matrix':matrix_fpath})
        annotator(seq1)
        assert len(seq1.features) == 1
        assert seq1.features[0].type == 'orf'


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testiprscan_parse']
    unittest.main()
