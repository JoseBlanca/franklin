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

import unittest, tempfile, os
from Bio.SeqFeature import  FeatureLocation
from franklin.snv.snv_annotation import INVARIANT, SNP, INSERTION, DELETION
from franklin.seq.seq_annotation import (create_microsatellite_annotator,
                                         create_ortholog_annotator,
                                         create_description_annotator,
                                         create_orf_annotator,
                                         create_go_annotator,
                                         create_cdna_intron_annotator,
                                         create_prot_change_annotator)
#                                         create_polia_annotator)

from franklin.seq.seqs import SeqWithQuality, Seq, SeqFeature
from franklin.utils.misc_utils import DATA_DIR
from franklin.utils.cmd_utils import b2gpipe_runner, BLAST_TOOL

class AnnotationTests(unittest.TestCase):
    'Annotations tests'
    @staticmethod
    def test_orthologs_annotator():
        'It test the ortholog annotator'
        blast_fhand = open(os.path.join(DATA_DIR, 'melon_tair.xml'))
        reverse_blast_fhand = open(os.path.join(DATA_DIR, 'tair_melon.xml'))
        blast = {'blast':blast_fhand,
                 'subj_def_as_acc':True}
        reverse_blast = {'blast':reverse_blast_fhand,
                         'subj_def_as_acc':True}
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
    def xtest_polia_annotator():
        'It test the polia annotator'
        polia_annot = create_polia_annotator()

        seq = 'TCGCATCGATCATCGCAGATCGACTGATCGATCGATCAAAAAAAAAAAAAAAAAAAAAAA'
        seq1 = SeqWithQuality(seq=Seq(seq))
        polia_annot(seq1)
        print seq1.features

        seq = 'TTTTTTTTTTTTTTTTTTTTCGCATCGATCATCGCAGATCGACTGATCGAGTAGTGATGT'
        seq1 = SeqWithQuality(seq=Seq(seq))
        polia_annot(seq1)
        print seq1.features


        seq = 'CGCATCGATCATCAAAAAAAAAAAAAAAAAAAAAAAAAAAGCAGATCGACTGATCGAGTA'
        seq1 = SeqWithQuality(seq=Seq(seq))
        polia_annot(seq1)
        print seq1.features

    @staticmethod
    def test_snv_prot_change_annotator():
        'It test the snv_prot_changepolia annotator'
        # first annotate orf
        seq  = 'ATGGCTTCATCCATTCTCTCATCCGCCGNTGTGGCCTTTGNCAACAGGGCTTCCCCTGCTCA'
        seq += 'AGCTAGCATGGGGGCACCATTCACTGGCCTAAAATCCGCCGCTGCTTTCCCNGTTTTATGTA'
        seq += 'CTGTTTTNACTCGCANGACCAACGACATCACCACTTTGGTTAGCAATGGGGGAAGAGTTCAG'
        seq += 'GGCNTGAAGGTGTGCCCACCACTTGGATTGAAGAAGTTCGAGACTCTTTCTTACCTTCCTGA'
        seq += 'TATGAGTAACGAGCAATTGGGAAAGGAAGTTGACTACCTTCTCAGGAAGGGATGGATTCCCT'
        seq += 'GCATTGAATTCGACATTCACAGTGGATTCGTTTACCGTGAGACCCACAGGTCACCAGGATAC'
        seq += 'TTCGATGGACGCTACTGGACCATGTGGAAGCTGCCCATGTTTGGCTGCACCGAT'

        orf  = 'ATGGCTTCATCCATTCTCTCATCCGCCGNTGTGGCCTTTGNCAACAGGGCTTCCCTGCTCAA'
        orf += 'GCTAGCATGGGGGCACCATTCACTGGCCTAAAATCCGCCGCTGCTTTCCCNGTNACTCGCAN'
        orf += 'GACCAACGACATCACCACTTTGGTTAGCAATGGGGGAAGAGTTCAGGGCNTGAAGGTGTGCC'
        orf += 'CACCACTTGGATTGAAGAAGTTCGAGACTCTTTCTTACCTTCCTGATATGAGTAACGAGCAA'
        orf += 'TTGGGAAAGGAAGTTGACTACCTTCTCAGGAAGGGATGGATTCCCTGCATTGAATTCGACAT'
        orf += 'TCACAGTGGATTCGTTTACCGTGAGACCCACAGGTCACCAGGATACTTCGATGGACGCTAC'
        orf += 'TGGACCATGTGGAAGCTGCCCATGTTTGGCTGCACCGAT'

        qualifiers = {'strand':'forward', 'dna':Seq(orf), 'prot':'prot'}
        feature = SeqFeature(location=FeatureLocation(0, len(seq)), type='orf',
                             qualifiers=qualifiers)
        alleles = {('A', SNP): None, ('G', INVARIANT):None}
        snv = SeqFeature(location=FeatureLocation(24, 24), type='snv',
                          qualifiers = {'alleles':alleles})
        alleles = {('A', SNP): None, ('T', INVARIANT):None}
        snv2 = SeqFeature(location=FeatureLocation(56, 56), type='snv',
                          qualifiers = {'alleles':alleles})
        alleles = {('T', SNP): None, ('A', INVARIANT):None}
        snv3 = SeqFeature(location=FeatureLocation(399, 399), type='snv',
                          qualifiers = {'alleles':alleles})
        alleles = {('T', INVARIANT):None, ('AG', INSERTION):None}
        snv4 = SeqFeature(location=FeatureLocation(250, 250), type='snv',
                          qualifiers = {'alleles':alleles})
        alleles = {('G', INVARIANT):None, ('--', DELETION):None}
        snv5 = SeqFeature(location=FeatureLocation(251, 251), type='snv',
                          qualifiers = {'alleles':alleles})



        sequence = SeqWithQuality(seq=Seq(seq), name='query')
        #print str(sequence[257:265].seq)
        #print sequence[256]
        sequence.features = [feature, snv, snv2, snv3, snv4, snv5]
        annotator = create_prot_change_annotator()
        annotator(sequence)
        [feature, snv, snv2, snv3, snv4, snv5] = sequence.features
        assert snv.qualifiers['protein_change']['kind'] == 'substitution'
        assert snv3.qualifiers['protein_change']['kind'] == 'substitution'
        assert snv4.qualifiers['protein_change']['kind'] == 'breakage'
#
        assert snv.qualifiers['protein_change']['location'] == 'codon_1'


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

    def test_intron_annotator(self):
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
        if BLAST_TOOL == 'blast':
            tomato_genome = 'tomato_genome2'
        else:
            tomato_genome = 'tomato_genome2+'
        genomic_db = os.path.join(blast_db_path, tomato_genome)
        intron_annotator = create_cdna_intron_annotator(genomic_db=genomic_db,
                                            genomic_seqs_fhand=open(genomic_db))
        intron_annotator(seq)
        intron_feat = seq.features[0]
        assert intron_feat.location.start.position == 478
        assert intron_feat.type == 'intron'

        #test an error
        genomic_seqs = tempfile.NamedTemporaryFile(suffix='.fasta')
        genomic_seqs.write('>a_seq\nATCGTG\n')
        intron_annotator = create_cdna_intron_annotator(genomic_db=genomic_db,
                                                genomic_seqs_fhand=genomic_seqs)
        try:
            intron_annotator(seq)
            self.fail('RuntimeError expected')
        except RuntimeError as error:
            assert genomic_seqs.name in str(error)

    @staticmethod
    def test_go_annotator():
        'It test the go annotator'
        blast = open(os.path.join(DATA_DIR, 'blastResult.xml'))
        prop_fpath = os.path.join(DATA_DIR, 'b2gPipe.properties')
        fhand, annot_fpath = tempfile.mkstemp()
        os.close(fhand)
        b2gpipe_runner(blast, annot_fpath, prop_fpath=prop_fpath)
        blast2go = annot_fpath
        go_annotator = create_go_annotator(blast2go)
        seq = SeqWithQuality(name='seq1', seq=Seq('aaaa'))

        go_annotator(seq)
        assert 'GO:0009853' in seq.annotations['GOs']

        os.remove(annot_fpath)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'AnnotationTests.test_get_description_with_funct']
    unittest.main()
