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

from os.path import join, exists

from franklin.utils.misc_utils import NamedTemporaryDir, DATA_DIR
from franklin.seq.seqs import SeqWithQuality, Seq
from franklin.seq.writers import write_seqs_in_file
from franklin.backbone.create_project import create_project
from franklin.backbone.backbone_runner import do_analysis
from franklin.backbone.analysis import BACKBONE_BASENAMES, BACKBONE_DIRECTORIES
from franklin.utils.cmd_utils import BLAST_TOOL
from franklin.seq.readers import seqs_in_file

THREADS = False

class AnnotationTest(unittest.TestCase):
    'It test the ortholog analysis'

    @staticmethod
    def test_ortholog_annotation_analysis():
        'We can annotate orthologs'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'

        config = {'blast':{'arabidopsis': {'path':'/path/to/tair',
                                           'species':'arabidopsis',
                                           'kind':'nucl'},
                          'arabidopsis2':{'path':'/path/to/tair2',
                                           'species':'arabidopsis2',
                                           'kind': 'nucl'}},

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
        pickle_fpath = join(project_dir, BACKBONE_DIRECTORIES['annotation_dbs'],
                          'melon.0.pickle')
        pickle = open(pickle_fpath).read()
        assert 'arabidopsis-orthologs' in pickle
        assert 'arabidopsis2-orthologs' in pickle

        do_analysis(project_settings=settings_path, kind='write_annotations',
                    silent=True)

        orf_fpath = join(project_dir, 'annotations', 'features',
                         'melon.orthologs')
        assert os.path.exists(orf_fpath)
        assert "tair1" in open(orf_fpath).read()

        orf_fpath = join(project_dir, 'annotations', 'features', 'melon.orf')
        assert not os.path.exists(orf_fpath)

        do_analysis(project_settings=settings_path, kind='annotation_stats',
                    silent=True)
        stats_fpath = join(project_dir, 'annotations', 'features', 'stats',
                           'melon.txt')
        result = open(stats_fpath).read()
        expected = '''Orthologs
_________
Sequences with arabidopsis orthologs: 2
Number of arabidopsis orthologs: 2
Sequences with arabidopsis2 orthologs: 2
Number of arabidopsis2 orthologs: 2'''
        #print result
        assert expected in result

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
        repr_fpath = join(project_dir, BACKBONE_DIRECTORIES['annotation_dbs'],
                          'melon.0.pickle')
        result = open(repr_fpath).read()
        assert 'yet another one' in result

        do_analysis(project_settings=settings_path, kind='annotation_stats',
                    silent=True)
        stats_fpath = join(project_dir, 'annotations', 'features', 'stats',
                           'melon.txt')
        result = open(stats_fpath).read()
        expected = '''Annotation statistics
---------------------
Number of sequences: 2
Sequences with description: 1'''
        assert expected in result

        os.chdir('/tmp')
        test_dir.close()

    @staticmethod
    def test_cdna_intron_annoation_analysis():
        'We can annotate introns'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'
        blast_db_path = os.path.join(DATA_DIR, 'blast')
        if BLAST_TOOL == 'blast':
            genomic_db = os.path.join(blast_db_path, 'tomato_genome2')
        else:
            genomic_db = os.path.join(blast_db_path, 'tomato_genome2+')
        config = {'Annotation':
                        {'Cdna_intron_annotation':{'genomic_db': genomic_db,
                                                'genomic_seq_file':genomic_db}},
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
        pickle_fpath = join(project_dir, BACKBONE_DIRECTORIES['annotation_dbs'],
                          'seqs.0.pickle')
        assert 'intron' in open(pickle_fpath).read()

        do_analysis(project_settings=settings_path, kind='annotation_stats',
                    silent=True)
        stats_fpath = join(project_dir, 'annotations', 'features', 'stats',
                           'seqs.txt')
        result = open(stats_fpath).read()
        expected = '''Sequences with intron: 1
Number of introns: 3'''
        assert expected in result

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
        pickle_fpath = join(project_dir, BACKBONE_DIRECTORIES['annotation_dbs'],
                          'seqs.0.pickle')
        result = open(pickle_fpath).read()
        assert 'microsatellite' in  result

        do_analysis(project_settings=settings_path, kind='write_annotations',
                    silent=True)
        ssr_fpath = join(project_dir, 'annotations', 'features', 'seqs.ssr')
        assert os.path.exists(ssr_fpath)
        assert "Seqname"  in open(ssr_fpath).read()

        do_analysis(project_settings=settings_path, kind='annotation_stats',
                    silent=True)
        stats_fpath = join(project_dir, 'annotations', 'features', 'stats',
                           'seqs.txt')
        result = open(stats_fpath).read()
        expected = '''Sequences with microsatellites: 1
Microsatellite types:
\tdinucleotide: 1
Microsatellite locations:
\tunknown: 1'''
        assert expected in result

        os.chdir('/tmp')
        test_dir.close()

    @staticmethod
    def test_orf_annotation_analysis():
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
        repr_fpath = join(project_dir, BACKBONE_DIRECTORIES['annotation_dbs'],
                          'seqs.0.pickle')
        result = open(repr_fpath).read()
        assert 'orf' in  result
        do_analysis(project_settings=settings_path, kind='write_annotations',
                    silent=True)

        seq_fpath = join(project_dir, 'annotations', 'features',
                         'seqs.orf_seq.fasta')
        pep_fpath = join(project_dir, 'annotations', 'features',
                         'seqs.orf_pep.fasta')

        assert 'ATCCGCCGNTGTGGCCTTTGNCAACAGGGCTTCCCCT' in open(seq_fpath).read()
        assert 'QASMGAPFTGLKSAAAFPVTRXTNDITTLVSNG' in open(pep_fpath).read()

        do_analysis(project_settings=settings_path, kind='annotation_stats',
                    silent=True)
        stats_fpath = join(project_dir, 'annotations', 'features', 'stats',
                           'seqs.txt')
        result = open(stats_fpath).read()
        expected = '''Sequences with ORF: 1
Number of ORFs: 1'''
        assert expected in result

        os.chdir('/tmp')
        test_dir.close()

    @staticmethod
    def test_go_annotation_analysis():
        'We can annotate gos'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'
        nr_path = os.path.join(DATA_DIR, 'blast', 'arabidopsis_genes+')
        b2g = os.path.join(DATA_DIR, 'b2gPipe.properties')
        config = {'blast':{'nr': {'path': nr_path,
                                           'species':'nr'}},
                  'Annotation':{'go_annotation':{'blast_database':'nr',
                                                 'create_dat_file':True,
                                                 'java_memory':2048,
                                                 'b2g_properties_file':b2g,
                                                 'blast2go_path':None}
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
        repr_fpath = join(project_dir, BACKBONE_DIRECTORIES['annotation_dbs'],
                          'seqs.0.pickle')
        result = open(repr_fpath).read()
        assert 'GO:0019253' in result
        assert os.path.exists(os.path.join(project_dir, 'annotations',
                                           'features', 'seqs.b2g.dat'))
        assert os.path.exists(os.path.join(project_dir, 'annotations',
                                           'features', 'seqs.b2g.annot'))
        do_analysis(project_settings=settings_path, kind='annotate_gos',
                    silent=True)

        do_analysis(project_settings=settings_path, kind='annotation_stats',
                    silent=True)
        stats_fpath = join(project_dir, 'annotations', 'features', 'stats',
                           'seqs.txt')
        result = open(stats_fpath).read()
        expected = '''Sequences with GOs: 1
Number of GOs: 12'''
        assert expected in result

    @staticmethod
    def test_protein_change_annotation_analysis():
        'We can annotate protein changes'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'
        matrix = os.path.join(DATA_DIR, 'At.smat')
        configuration = {'Snvs':{'min_quality':20},
                         'Sam_processing':{'add_default_qualities':True},
                         'Annotation':{'orf_annotation':
                                                     {'estscan_matrix':matrix}},
                         'General_settings':{'threads':THREADS}}

        settings_path = create_project(directory=test_dir.name,
                                       name=project_name,
                                       configuration=configuration)
        project_dir = join(test_dir.name, project_name)
        #setup the original reads
        reads_dir = join(project_dir, 'reads')
        clean_reads_dir = join(reads_dir, 'cleaned')
        os.mkdir(reads_dir)
        os.mkdir(clean_reads_dir)

        solexa = '@seq1\n'
        solexa += 'TCATTGAAAGTTGAAACTGATAGTAGCAGAGTTTTTTCCTCTGTTTGG\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIIIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa = '@seq2\n'
        solexa += 'ATATGATTGAAGATATTTCTGGGCTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa = '@seq14\n'
        solexa += 'ATATGATTGAAGATATTTCTGGGCTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa = '@seq15\n'
        solexa += 'ATATGATTGAAGATATTTCTGGGCTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa = '@seq12\n'
        solexa += 'ATATGATTGAAGATATTTCTGGACTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa = '@seq13\n'
        solexa += 'ATATGATTGAAGATATTTCTGGACTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa = '@seq16\n'
        solexa += 'ATATGATTGAAGATATTTCTGGACTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'

        sanger = '>seq3\n'
        sanger += 'GATATGATTGAAGATATTTCTGGGCTTTAAGGGTTCTTGAGGATTTATAGGAGATACTGA'
        sanger += 'GATTCTGGAATCTCTGAGTTTCTGGGTTCAAGTTGCACTGACCATTGTTGGATTTGTAGA'
        sanger += 'TTGTTTCTTCTTTCATTAGGCATTGATTATGGGTAAATGCGTGGGTACATATAATATATA'
        sanger += 'TCTGTTGAATGCAATTTACACATTGACTGAGGAACAACATGAACATGGCAGCTTTCTCAA'
        sanger += 'AATTGAACCACAGAAGGCTTAAAAGCAAAGTCTTTGGAGAATCAGACTAAGCTTGAGA\n'
        sanger += '>seq4\n'
        sanger += 'TCCATACTTTACTCTATCTCTTTCTGTGTTTGGTAACACGCGAGGATTGGATGATAT'
        sanger += 'CTATATATTCTAATGTGGACTAAAAATGTGTGTGTGTGTATGAAGATGGGAAGCCGGAAG'
        sanger += 'TCATCAAGGAGAAAAGAGAGATAGACGACGAGAAGATGGCGTTGACGTTCAGAGGACTAG'
        sanger += 'AGGGTCATGTGATGGAGAAGTACAAGAAGTATGAGGTTATCTTACAGTTCATTCCCAAGT'
        sanger += 'CGAACGAAGGCTGCGTCTGCAAAGTCACTCTGATATGGGAGAATCGCAACGAAGACTCCC'
        sanger += '>seq5\n'
        sanger += 'TCCATACTTTACTCTATCTCTTTCTGTGTTTGGTAACACGCGAGGATTGGATGATAT'
        sanger += 'CTATATATTCTAAAGTGGACTAAAAATGTGTGTGTGTGTATGAAGATGGGAAGCCGGAAG'
        sanger += 'TCATCAAGGAGAAAAGAGAGATAGACGACGAGAAGATGGCGTTGACGTTCAGAGGACTAG'
        sanger += 'AGGGTCATGTGATGGAGAAGTACAAGAAGTATGAGGTTATCTTACAGTTCATTCCCAAGT'
        sanger += 'CGAACGAAGGCTGCGTCTGCAAAGTCACTCTGATATGGGAGAATCGCAACGAAGACTCCC'

        fpath_sanger = join(clean_reads_dir, 'lb_hola1.pl_sanger.sm_hola.fasta')
        fpath_solexa = join(clean_reads_dir,
                                    'lb_hola2.pl_illumina.sm_hola.sfastq')
        open(fpath_sanger, 'w').write(sanger)
        open(fpath_solexa, 'w').write(solexa)

        fpath_sanger2 = join(clean_reads_dir, 'lb_adios.pl_sanger.fasta')
        fpath_solexa2 = join(clean_reads_dir,
                                    'lb_adios.pl_illumina.sfastq')
        open(fpath_sanger2, 'w').write(sanger)
        open(fpath_solexa2, 'w').write(solexa)

        #the reference
        reference_dir = join(project_dir, 'mapping/reference')
        os.makedirs(reference_dir)
        reference_fpath = join(reference_dir, 'reference.fasta')
        out = open(reference_fpath, 'w')
        for line in open(join(DATA_DIR, 'blast/arabidopsis_genes')):
            out.write(line)

        do_analysis(project_settings=settings_path, kind='mapping', silent=True)
        mapping_dir = join(project_dir, 'mapping')
        singular_mapping_dir = sorted(os.listdir(mapping_dir))[0]
        singular_mapping_dir = join(mapping_dir, singular_mapping_dir)
        assert exists(join(singular_mapping_dir, 'bams',
                            'by_readgroup', 'lb_hola2.pl_illumina.sm_hola.bam'))
        result_dir = join(mapping_dir, 'bams')
        assert exists(result_dir)
        result_dir_by_lib = join(result_dir, 'by_readgroup')
        assert exists(result_dir_by_lib)

        do_analysis(project_settings=settings_path, kind='merge_bams',
                    silent=True)
        assert exists(join(result_dir, 'merged.0.bam'))
        assert exists(join(result_dir, 'merged.0.bam.bai'))

        #we realign the mapping using GATK
        do_analysis(project_settings=settings_path, kind='realign_bam',
                    silent=True)
        assert exists(join(result_dir, 'merged.1.bam'))

        annot_input_dir = join(project_dir, 'annotations', 'input')
        os.makedirs(annot_input_dir)
        os.symlink(reference_fpath, join(annot_input_dir, 'reference.fasta'))
        do_analysis(project_settings=settings_path, kind='annotate_snvs',
                    silent=True)

        do_analysis(project_settings=settings_path, kind='annotate_orfs',
                    silent=True)

        do_analysis(project_settings=settings_path, kind='annotate_prot_change',
                    silent=True)

        result_file = join(project_dir, 'annotations', 'db',
                           'reference.2.pickle')

        seqs = list(seqs_in_file(open(result_file)))
        snv = seqs[2].features[0]
        assert snv.qualifiers['protein_change']['kind'] == 'substitution'
        assert snv.qualifiers['protein_change']['location'] == 'codon_1'

        os.chdir('/tmp')
        test_dir.close()



if    __name__ == "__main__":
    #import sys;sys.argv = ['', 'AnnotationTest.test_protein_change_annotation_analysis']
    unittest.main()
