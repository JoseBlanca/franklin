'''
Created on 16/03/2010

@author: jose
'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of project.
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

import unittest, os, shutil

from os.path import join

from franklin.utils.misc_utils import NamedTemporaryDir, DATA_DIR
from franklin.seq.seqs import SeqWithQuality, Seq
from franklin.seq.writers import write_seqs_in_file
from franklin.backbone.create_project import create_project
from franklin.backbone.backbone_runner import do_analysis
from franklin.backbone.analysis import BACKBONE_BASENAMES, BACKBONE_DIRECTORIES

THREADS = False

class AnnotationTest(unittest.TestCase):
    'It test the ortholog analysis'

    @staticmethod
    def test_ortholog_annotation_analysis():
        'We can annotate orthologs'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'

        config = {'blast':{'arabidopsis': {'path':'/path/to/tair',
                                           'species':'arabidopsis'},
                          'arabidopsis2':{'path':'/path/to/tair2',
                                           'species':'arabidopsis2'}},

                  'Annotation':{'ortholog_annotation':{'ortholog_databases':
                                            ['arabidopsis', 'arabidopsis2']}},
                'General_settings':{'threads':THREADS}}


        settings_path = create_project(directory=test_dir.name,
                                       name=project_name,
                                       configuration=config)
        project_dir = join(test_dir.name, project_name)

        # create blast results
        melon_tair_blastdir = join(project_dir, 'annotations', 'blast',
                                   'melon', 'tair')
        melon_tair2_blastdir = join(project_dir, 'annotations', 'blast',
                                   'melon', 'tair2')
        os.makedirs(melon_tair_blastdir)
        os.makedirs(melon_tair2_blastdir)
        tair_melon_blastdir = join(project_dir, 'annotations', 'blast',
                                   'tair', 'melon')
        tair2_melon_blastdir = join(project_dir, 'annotations', 'blast',
                                   'tair2', 'melon')
        os.makedirs(tair_melon_blastdir)
        os.makedirs(tair2_melon_blastdir)
        blast_fname = BACKBONE_BASENAMES['blast_basename'] + '.tblastx.xml'
        shutil.copy(join(DATA_DIR, 'melon_tair.xml'),
                   join(melon_tair_blastdir, blast_fname))
        shutil.copy(join(DATA_DIR, 'melon_tair.xml'),
                   join(melon_tair2_blastdir, blast_fname))
        shutil.copy(join(DATA_DIR, 'tair_melon.xml'),
                   join(tair_melon_blastdir, blast_fname))
        shutil.copy(join(DATA_DIR, 'tair_melon.xml'),
                   join(tair2_melon_blastdir, blast_fname))

        #some melon file to annotate
        input_dir = join(project_dir, BACKBONE_DIRECTORIES['annotation_input'])
        os.makedirs(input_dir)
        seq1 = SeqWithQuality(Seq('A'), id='melon1')
        seq2 = SeqWithQuality(Seq('A'), id='melon2')
        write_seqs_in_file([seq1, seq2],
                           open(join(input_dir, 'melon.fasta'), 'a'))

        do_analysis(project_settings=settings_path, kind='annotate_orthologs',
                    silent=True)
        repr_fpath = join(project_dir, 'annotations', 'repr', 'melon.0.repr')
        repr_ = open(repr_fpath).read()
        assert 'arabidopsis-orthologs' in repr_
        assert 'arabidopsis2-orthologs' in repr_

        do_analysis(project_settings=settings_path, kind='write_annotations',
                    silent=True)

        ort_fpath = join(project_dir, 'annotations', 'result', 'melon.orthologs')
        assert os.path.exists(ort_fpath)
        assert "tair1" in open(ort_fpath).read()

        orf_fpath = join(project_dir, 'annotations', 'result', 'melon.orf')
        assert not os.path.exists(orf_fpath)


        os.chdir('/tmp')
        test_dir.close()

    @staticmethod
    def test_description_annotation_analysis():
        'We can annotate with description'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'
        arab_blastdb = join(DATA_DIR, 'blast', 'arabidopsis_genes+')
        config = {'blast':{'arabidopsis': {'path': arab_blastdb,
                                           'species':'arabidopsis'}},
                  'Annotation':{'description_annotation':{
                                                    'description_databases':
                                                              ['arabidopsis']}},
                'General_settings':{'threads':THREADS}
                 }

        settings_path = create_project(directory=test_dir.name,
                                       name=project_name,
                                       configuration=config)
        project_dir = join(test_dir.name, project_name)

        #some melon file to annotate
        input_dir = join(project_dir, BACKBONE_DIRECTORIES['annotation_input'])
        os.makedirs(input_dir)
        seq1 = SeqWithQuality(Seq('GGGTCAGCGAAGAACCACTACAAGAGGTGGAAGAGCGAAA'),
                              id='CUTC021854')
        seq2 = SeqWithQuality(Seq('Atagtagcatcagatgagcatcgacttctagctagctagct'),
                               id='CUTC021853')
        write_seqs_in_file([seq1, seq2],
                           open(join(input_dir, 'melon.fasta'), 'a'))

        do_analysis(project_settings=settings_path,
                    kind='annotate_descriptions', silent=True)
        repr_fpath = join(project_dir, 'annotations', 'repr', 'melon.0.repr')
        result = open(repr_fpath).read()
        assert 'yet another one' in result

        os.chdir('/tmp')
        test_dir.close()

    @staticmethod
    def test_cdna_intron_annoation_analysis():
        'We can annotate introns'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'
        blast_db_path = os.path.join(DATA_DIR, 'blast')
        genomic_db = os.path.join(blast_db_path, 'tomato_genome2')
        config = {'Annotation':
                        {'Cdna_intron_annotation':{'genomic_db': genomic_db,
                                                   'genomic_seqs':genomic_db}},
                'General_settings':{'threads':THREADS}}
        settings_path = create_project(directory=test_dir.name,
                                       name=project_name,
                                      configuration=config)
        project_dir = join(test_dir.name, project_name)
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
        annot_input_dir = join(project_dir, 'annotations', 'input')
        os.makedirs(annot_input_dir)

        #create some seqs to annotate
        fasta = '>seq\n%s\n' % seq
        fhand = open(os.path.join(annot_input_dir, 'seqs.fasta'), 'w')
        fhand.write(fasta)
        fhand.close()
        do_analysis(project_settings=settings_path, kind='annotate_introns',
                    silent=True)
        repr_fpath = join(project_dir, 'annotations', 'repr', 'seqs.0.repr')

        assert "type='intron'" in  open(repr_fpath).read()
        os.chdir('/tmp')
        test_dir.close()

    @staticmethod
    def test_microsatellite_annoation_analysis():
        'We can annotate introns'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'
        settings_path = create_project(directory=test_dir.name,
                                       name=project_name,
                         configuration={'General_settings':{'threads':THREADS}})
        project_dir = join(test_dir.name, project_name)
        seq = 'GAAAAGATGTGATTGGTGAAATAAGTTTGCCTCAATTCTCTTGTGCCGAAGTTCCAAAGAAGC'
        seq += 'AGTTGGTGAATGAGCAGCCAGTACCCGAAAAATCGAGCAAAGATTTTGTGATGTATGTTGGAG'
        seq += 'GTCTAGCATGGGGGATGGACTGGTGTCCCCAAGCTCATGAAAATAGGGATGCTCCTATGAAAA'
        seq += 'GAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA'
        seq += 'GTGAGTTTGTCGCAATTGCTCCTCATCCTCCTGATTCATCATATCACAAGACTGATGCCTCAC'
        seq += 'TTACAGGCAGAGGTGTAATTCAGATATGGTGCCTGCCAGATCTCATTCAAAAAGATATAATTG'
        seq += 'TGAAAGAAGATTATTTTGCTCAGGTTAACAAAAAACCGTATAGAAATTTGACAAGAAGTGAAG'
        seq += 'CAGGTACGGGAGAAGTATCTGGACCTCAAAAACCAAGAGGAAGACCAAAAAAGAACCCTGGTA'
        seq += 'AAGCAGTCCAGGCAAAAGCATCTAGACCACAAAATCCAAGAGGAAGACCGAGAAAGAAGCCTG'
        seq += 'TTACTGAATCTTTAGGTGATAGAGATAGTGAAGACCACAGTTTACAACCTCTTGCTATAGAGT'
        seq += 'GGTCGCTGCAATCAACAGAACTTTCTGTAGATTTGTCTTGTGGAAATATGAATAAAGCCCAAG'
        seq += 'TAGATATTGCGCTGAGTCAAGAAAGATGTATTAATGCGGCAT'
        annot_input_dir = join(project_dir, 'annotations', 'input')
        os.makedirs(annot_input_dir)

        #create some seqs to annotate
        fasta = '>seq\n%s\n' % seq
        fhand = open(os.path.join(annot_input_dir, 'seqs.fasta'), 'w')
        fhand.write(fasta)
        fhand.close()
        do_analysis(project_settings=settings_path,
                    kind='annotate_microsatellites', silent=True)
        repr_fpath = join(project_dir, 'annotations', 'repr', 'seqs.0.repr')
        result = open(repr_fpath).read()
        assert "type='microsatellite'" in  result

        do_analysis(project_settings=settings_path, kind='write_annotations',
                    silent=True)
        ssr_fpath = join(project_dir, 'annotations', 'result', 'seqs.ssr')
        assert os.path.exists(ssr_fpath)
        assert "Seqname"  in open(ssr_fpath).read()


        os.chdir('/tmp')
        test_dir.close()

    @staticmethod
    def test_orf_annoation_analysis():
        'We can annotate orfs'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'
        matrix = os.path.join(DATA_DIR, 'At.smat')
        config = {'Annotation':{'orf_annotation': {'estscan_matrix':matrix}},
                  'General_settings':{'threads':THREADS}}

        settings_path = create_project(directory=test_dir.name,
                                       name=project_name,
                                       configuration=config)
        project_dir = join(test_dir.name, project_name)
        seq = 'CTACTTACTAGCTTTAGTAAATCCTTCTAACCCTCGGTAAAAAAAAAAAAGAGGCATCAAATG'
        seq += 'GCTTCATCCATTCTCTCATCCGCCGNTGTGGCCTTTGNCAACAGGGCTTCCCCTGCTCAAGCT'
        seq += 'AGCATGGGGGCACCATTCACTGGCCTAAAATCCGCCGCTGCTTTCCCNGTNACTCGCANGACC'
        seq += 'AACGACATCACCACTTTGGTTAGCAATGGGGGAAGAGTTCAGGGCNTGAAGGTGTGCCCACCA'
        seq += 'CTTGGATTGAAGAAGTTCGAGACTCTTTCTTACCTTCCTGATATGAGTAACGAGCAATTGGGA'
        seq += 'AAGGAAGTTGACTACCTTCTCAGGAAGGGATGGATTCCCTGCATTGAATTCGACATTCACAGT'
        seq += 'GGATTCGTTTACCGTGAGACCCACAGGTCACCAGGATACTTCGATGGACGCTACTGGACCATG'
        seq += 'TGGAAGCTGCCCATGTTTGGCTGCACCGAT'

        annot_input_dir = join(project_dir, 'annotations', 'input')
        os.makedirs(annot_input_dir)

        #create some seqs to annotate
        fasta = '>seq\n%s\n' % seq
        fhand = open(os.path.join(annot_input_dir, 'seqs.fasta'), 'w')
        fhand.write(fasta)
        fhand.close()
        do_analysis(project_settings=settings_path,
                    kind='annotate_orfs', silent=True)
        repr_fpath = join(project_dir, 'annotations', 'repr', 'seqs.0.repr')
        result = open(repr_fpath).read()
        assert "type='orf'" in  result
        do_analysis(project_settings=settings_path, kind='write_annotations',
                    silent=True)
        seq_fpath = join(project_dir, 'annotations', 'result', 'seqs.orf_seq.fasta')
        pep_fpath = join(project_dir, 'annotations', 'result', 'seqs.orf_pep.fasta')

        assert 'ATCCGCCGNTGTGGCCTTTGNCAACAGGGCTTCCCCT' in open(seq_fpath).read()
        assert 'QASMGAPFTGLKSAAAFPVTRXTNDITTLVSNG' in open(pep_fpath).read()

        os.chdir('/tmp')
        test_dir.close()

    @staticmethod
    def test_go_annotation_analysis():
        'We can annotate gos'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'
        nr_path = os.path.join(DATA_DIR, 'blast', 'arabidopsis_genes+')
        config = {'blast':{'nr': {'path': nr_path,
                                           'species':'nr'}},
                  'Annotation':{'go_annotation':{'blast_database':'nr',
                                                 'create_dat_file':True,
                                                 'java_memory':2048}
                 }, 'General_settings':{'threads':THREADS}}

        settings_path = create_project(directory=test_dir.name,
                                       name=project_name,
                                       configuration=config)
        project_dir = join(test_dir.name, project_name)
        seq = 'CTTCATCCATTCTCTCATCCGCCGNTGTGGCCTTTGNCAACAGGGCTTCCCCTGCTCAAGCT'
        seq += 'AGCATGGGGGCACCATTCACTGGCCTAAAATCCGCCGCTGCTTTCCCNGTNACTCGCANGACC'
        seq += 'AACGACATCACCACTTTGGTTAGCAATGGGGGAAGAGTTCAGGGCNTGAAGGTGTGCCCACCA'
        seq += 'CTTGGATTGAAGAAGTTCGAGACTCTTTCTTACCTTCCTGATATGAGTAACGAGCAATTGGGA'
        seq += 'AAGGAAGTTGACTACCTTCTCAGGAAGGGATGGATTCCCTGCATTGAATTCGACATTCACAGT'
        seq += 'GGATTCGTTTACCGTGAGACCCACAGGTCACCAGG'

        annot_input_dir = join(project_dir, 'annotations', 'input')
        os.makedirs(annot_input_dir)

        #create some seqs to annotate
        fasta = '>seq1\n%s\n' % seq
        fhand = open(os.path.join(annot_input_dir, 'seqs.fasta'), 'w')
        fhand.write(fasta)
        fhand.close()
        bdir = join(project_dir, 'annotations', 'blast', 'seqs',
                    'arabidopsis_genes+')
        os.makedirs(bdir)
        shutil.copy(join(DATA_DIR, 'blastResult.xml'),
                    join(bdir, 'blast.tblastx.xml'))

        do_analysis(project_settings=settings_path, kind='annotate_gos',
                    silent=True)
        repr_fpath = join(project_dir, 'annotations', 'repr', 'seqs.0.repr')
        result = open(repr_fpath).read()
        assert 'GO:0019253' in result
        assert os.path.exists(os.path.join(project_dir, 'annotations',
                                           'result', 'seqs.b2g.dat'))
        assert os.path.exists(os.path.join(project_dir, 'annotations',
                                           'result', 'seqs.b2g.annot'))
        do_analysis(project_settings=settings_path, kind='annotate_gos',
                    silent=True)

if    __name__ == "__main__":
#    import sys;sys.argv = ['', 'AnnotationTest.test_go_annotation_analysis']
    unittest.main()
