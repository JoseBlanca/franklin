'''
Created on 16/03/2010

@author: jose
'''

import unittest, os, shutil

from os.path import join

from franklin.utils.misc_utils import NamedTemporaryDir, DATA_DIR
from franklin.seq.seqs import SeqWithQuality, Seq
from franklin.seq.writers import write_seqs_in_file
from franklin.backbone.create_project import create_project
from franklin.backbone.backbone_runner import do_analysis
from franklin.backbone.analysis import BACKBONE_BASENAMES, BACKBONE_DIRECTORIES

class OrthologTest(unittest.TestCase):
    'It test the ortholog analysis'

    @staticmethod
    def test_ortholog_annotation_analysis():
        'We can annotate orthologs'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'

        config = {'blast':{'arabidopsis': {'path':'/path/to/tair',
                                    'species':'arabidopsis',
                                    'kind': 'nucl'}},
                  'Annotation':{'ortholog_annotation':{'ortholog_databases':
                                                       ['arabidopsis']}}
                 }

        settings_path = create_project(directory=test_dir.name,
                                       name=project_name,
                                       configuration=config)
        project_dir = join(test_dir.name, project_name)

        # create blast results
        melon_tair_blastdir = join(project_dir, 'annotations', 'blast',
                                   'melon', 'tair')
        os.makedirs(melon_tair_blastdir)
        tair_melon_blastdir = join(project_dir, 'annotations', 'blast',
                                   'tair', 'melon')
        os.makedirs(tair_melon_blastdir)
        blast_fname = BACKBONE_BASENAMES['blast_basename'] + '.tblastx.xml'
        shutil.copy(join(DATA_DIR, 'melon_tair.xml'),
                   join(melon_tair_blastdir, blast_fname))
        shutil.copy(join(DATA_DIR, 'tair_melon.xml'),
                   join(tair_melon_blastdir, blast_fname))

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
        assert 'arabidopsis-orthologs' in open(repr_fpath).read()

        os.chdir('/tmp')
        test_dir.close()

    @staticmethod
    def test_description_annotation_analysis():
        'We can annotate with description'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'
        arab_blastdb = join(DATA_DIR, 'blast', 'arabidopsis_genes+')
        config = {'blast':{'arabidopsis': {'path': arab_blastdb,
                                           'species':'arabidopsis',
                                           'kind': 'nucl'}},
                  'Annotation':{'description_annotation':{
                                                    'description_databases':
                                                              ['arabidopsis']}}
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

        do_analysis(project_settings=settings_path, kind='annotate_description',
                    silent=True)
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
                                                   'genomic_seqs':genomic_db}}}
        settings_path = create_project(directory=test_dir.name,
                                       name=project_name,
                                      configuration=config)
        project_dir = join(test_dir.name, project_name)
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
                                       name=project_name)
        project_dir = join(test_dir.name, project_name)
        seq  = 'GAAAAGATGTGATTGGTGAAATAAGTTTGCCTCAATTCTCTTGTGCCGAAGTTCCAAAGAAGC'
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
                    kind='annotate_microsatellite', silent=True)
        repr_fpath = join(project_dir, 'annotations', 'repr', 'seqs.0.repr')
        result = open(repr_fpath).read()
        assert "type='microsatellite'" in  result

        do_analysis(project_settings=settings_path, kind='write_annotation',
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
        config =  {'Annotation':{'orf_annotation': {'estscan_matrix':matrix}}}

        settings_path = create_project(directory=test_dir.name,
                                       name=project_name,
                                       configuration=config)
        project_dir = join(test_dir.name, project_name)
        seq  = 'CTACTTACTAGCTTTAGTAAATCCTTCTAACCCTCGGTAAAAAAAAAAAAGAGGCATCAAATG'
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
                    kind='annotate_orf', silent=True)
        repr_fpath = join(project_dir, 'annotations', 'repr', 'seqs.0.repr')
        result = open(repr_fpath).read()
        assert "type='orf'" in  result
        os.chdir('/tmp')
        test_dir.close()

    @staticmethod
    def test_go_annotation_analysis():
        'We can annotate gos'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'
        nr_path = os.path.join(DATA_DIR, 'blast', 'arabidopsis_genes+')
        config = {'blast':{'nr': {'path': nr_path,
                                           'species':'nr',
                                           'kind': 'nucl'}},
                  'Annotation':{'go_annotation':{'blast_database':'nr',
                                                 'create_dat_file':True}
                 }}

        settings_path = create_project(directory=test_dir.name,
                                       name=project_name,
                                       configuration=config)
        project_dir = join(test_dir.name, project_name)
        seq  = 'CTTCATCCATTCTCTCATCCGCCGNTGTGGCCTTTGNCAACAGGGCTTCCCCTGCTCAAGCT'
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

        do_analysis(project_settings=settings_path, kind='annotate_go',
                    silent=True)
        repr_fpath = join(project_dir, 'annotations', 'repr', 'seqs.0.repr')
        result = open(repr_fpath).read()
        print result
        assert 'GO:0019253' in result
        assert os.path.exists(os.path.join(project_dir, 'annotations',
                                           'go_files', 'seqs.dat'))
        print "remove this part of the test after you test it"
        return

        test_dir = NamedTemporaryDir()
        project_name = 'backbone'
        nr_path = '/srv/databases/blast/nr+'
        config = {'blast':{'nr': {'path': nr_path,
                                           'species':'nr',
                                           'kind': 'nucl'}},
                  'Annotation':{'go_annotation':{'blast_database':'nr'}}
                 }

        settings_path = create_project(directory=test_dir.name,
                                       name=project_name,
                                       configuration=config)
        project_dir = join(test_dir.name, project_name)
        seq  = 'CTTCATCCATTCTCTCATCCGCCGNTGTGGCCTTTGNCAACAGGGCTTCCCCTGCTCAAGCT'
        seq += 'AGCATGGGGGCACCATTCACTGGCCTAAAATCCGCCGCTGCTTTCCCNGTNACTCGCANGACC'
        seq += 'AACGACATCACCACTTTGGTTAGCAATGGGGGAAGAGTTCAGGGCNTGAAGGTGTGCCCACCA'
        seq += 'CTTGGATTGAAGAAGTTCGAGACTCTTTCTTACCTTCCTGATATGAGTAACGAGCAATTGGGA'
        seq += 'AAGGAAGTTGACTACCTTCTCAGGAAGGGATGGATTCCCTGCATTGAATTCGACATTCACAGT'
        seq += 'GGATTCGTTTACCGTGAGACCCACAGGTCACCAGG'

        annot_input_dir = join(project_dir, 'annotations', 'input')
        os.makedirs(annot_input_dir)

        #create some seqs to annotate
        fasta = '>seq\n%s\n' % seq
        fhand = open(os.path.join(annot_input_dir, 'seqs.fasta'), 'w')
        fhand.write(fasta)
        fhand.close()
        bdir = join(project_dir, 'annotations', 'blast', 'seqs',
                    'arabidopsis_genes+')

        do_analysis(project_settings=settings_path, kind='annotate_go')
        repr_fpath = join(project_dir, 'annotations', 'repr', 'seqs.repr')
        result = open(repr_fpath).read()
        print result




if    __name__    ==    "__main__":
    #import    sys;sys.argv    =    ['',    'SamTest.test_realignbam']
    unittest.main()
