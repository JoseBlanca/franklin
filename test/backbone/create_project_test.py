'''
Created on 26/01/2010

@author: jose
'''
import unittest, os.path
from os.path import join, exists

from configobj import ConfigObj

from franklin.utils.misc_utils import NamedTemporaryDir, DATA_DIR
from franklin.backbone.create_project import create_project
from franklin.backbone.analysis import BACKBONE_DIRECTORIES
from franklin.backbone.backbone_runner import do_analysis

READS_454 = '''@FKU4KFK07H6D2L
GGTTCAAGGTTTGAGAAAGGATGGGAAGAAGCCAAATGCCTACATTGCTGATACCACTACGGCAAATGCTCAAGTTCGGACGCTTGCTGAGACGGTGAGACTGGATGCAAGAACTAAGTTATTGAAT
+
CCFFFFFCC;;;FFF99;HHIHECCHHIIIFHHEEEHHHHHHIIIFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFGGGGFFGGBCCGFFFGG
@FKU4KFK07IJEB9
AATTCCCTTTCCGGTGATCGCCTCTCCTAGCAATGATGCAGCAAAGCCGATCATAGCAACACGGCCGACGAAGAGTTCATTCTGCTTGGTAAAACCAATGCCACCAGAAGTCCCAAAAATTCCATCTTCAACCTTCTGCTTCGGCTTCTCAACCTTCTTGGGCGGGGCTTTAGTTCTGGATTTGAACACAGCAAGAGGGTGAAACCA
+
FFFFFFFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFGGGGFFFFGGGGGGFFFFFFGAA@5557>55577G
@FKU4KFK07IMD9G
ACCCTTCTAGCAGTAGGAATGACTTGACCACCTCCTCTGTGGATGGCATCAGCGTGAAGCACCACATCACAAACTTCAAAACAGATGCCT
+
B>>>GGFFFFFFFFFGGGGHHHHHHHHHIIIIIIIIIIIIIIIIIFFFFFFFFFFFFFFFFFFFFFFFFFFFFFGCB;;;;BGFFFFFFF
@FKU4KFK07IE66H
ATGAACCGCAAAAGCTTGTATGCTGTATTGCCTTGATTTGGTTTCCAAGATTCTTCCCACATATTTAGGAGAGAGTGTAGCTGGAAACCGAGACTTTGGCTGGAGAAGCA
+
FFFFFFFFG7777CGFFFFHHHHHHIHHHHHHHHHHHHHHH==?CCFFFFFFFGGCCCGFFFFFFFFFFFFFFFFFFFFFFFFFGGGFFFFFFFFFFFFFGGGFFFFFFF
'''
READS_ILL= '''@HWI-EAS59:3:1:5:1186#0/1
TGATACCACTGCTTANTCTGCGTTGNTACCA
+
BBBBCBAAABA=BB<%<C@?BBA@9%<B;A@
@HWI-EAS59:3:1:5:169#0/1
GTTGATACCACTGCTNACTCTGCG
+
BBA@ABBBBBBBAA9%;BBBBCAA
@HWI-EAS59:3:1:5:1467#0/1
AAGCAGTGGTATCAANGCAGAGTTCNGCTTC
+
B?@CBCBCBAACBB:%>BBBAC@B<%9>BBB
@HWI-EAS59:3:1:5:651#0/1
CATACCCAATTGCTTNAGATGGAT
+
AAAABBA@@B?7=@<%=<ABBA:A
@HWI-EAS59:3:1:5:1609#0/1
AGAGCATTGCTTTGANTTTAAAAACNCGGTT
+
?BA@B>CCCBCBB:<%>BCA@;A?<%?BBAB
@HWI-EAS59:3:1:5:247#0/1
GATGGATCCCAAGTTNTTGAGGAACNAGAGG
+
BBA?;BBBBBA3AB=%=BBB@A=A=%<><@?
'''

class TestBackbone(unittest.TestCase):
    'It tests the backbone'

    @staticmethod
    def test_create_project():
        'We can create a project'
        test_dir = NamedTemporaryDir()
        settings_path = create_project(directory=test_dir.name,
                                       name='backbone')

        assert settings_path == join(test_dir.name,
                                'backbone', BACKBONE_DIRECTORIES['config_file'])
        settings = ConfigObj(settings_path, unrepr=True)
        assert settings['General_settings']['project_name'] == 'backbone'
        project_path = join(test_dir.name, 'backbone')
        assert settings['General_settings']['project_path'] == project_path
        os.chdir('/tmp')
        test_dir.close()

    @staticmethod
    def xtest_cleaning_analysis_lucy():
        'We can clean the reads'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'
        settings_path = create_project(directory=test_dir.name,
                                       name=project_name)
        project_dir = join(test_dir.name, project_name)
        #setup the original reads
        reads_dir = join(project_dir, 'reads')
        original_reads_dir = join(reads_dir, 'original')
        os.mkdir(reads_dir)
        os.mkdir(original_reads_dir)

        os.makedirs(join(project_dir, 'config_data', 'lucy'))
        lucy_settings = join(project_dir, 'config_data', 'lucy', 'lucy.conf')
        luc_c = open(lucy_settings, 'w')
        luc_c.write(repr({'ps':{'vector_file':'tmp' , 'splice_file':'tmp'}}))
        luc_c.flush()

        #print original_reads_dir
        fpath_454 = join(original_reads_dir, 'pl_454.lb_ps.sfastq')
        fpath_ill = join(original_reads_dir, 'pl_illumina.lb_psi.sfastq')
        open(fpath_454, 'w').write(READS_454)
        open(fpath_ill, 'w').write(READS_ILL)
        do_analysis(project_settings=settings_path, kind='clean_reads')
        cleaned_dir = join(project_dir, 'reads', 'cleaned')
        assert exists(cleaned_dir)

        cleaned_454 = join(cleaned_dir, os.path.basename(fpath_454))
        assert exists(cleaned_454)

    @staticmethod
    def test_cleaning_analysis():
        'We can clean the reads'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'
        project_dir = join(test_dir.name, project_name)
        adaptors_dir = join(project_dir, 'config_data', 'adaptors')
        adaptors_path_454 = join(adaptors_dir, '454_adaptors')
        configuration = {'Cleaning':{'adaptors_file_454':adaptors_path_454}}
        settings_path = create_project(directory=test_dir.name,
                                       name=project_name,
                                       configuration=configuration)

        #setup the original reads
        reads_dir = join(project_dir, 'reads')
        original_reads_dir = join(reads_dir, 'original')
        os.mkdir(reads_dir)
        os.mkdir(original_reads_dir)

        os.makedirs(adaptors_dir)
        adap_fhand = open(adaptors_path_454, 'w')
        adap_fhand.write('''>smart_5_cds_primer_1
GGTTCAAGGTTTGAGAAAGGATGGGAAG''')
        adap_fhand.close()

        #print original_reads_dir
        fpath_454 = join(original_reads_dir, 'pl_454.lb_a.sfastq')
        fpath_ill = join(original_reads_dir, 'pl_illumina.lb_b.sfastq')
        open(fpath_454, 'w').write(READS_454)
        open(fpath_ill, 'w').write(READS_ILL)

        do_analysis(project_settings=settings_path, kind='clean_reads')
        cleaned_dir = join(project_dir, 'reads', 'cleaned')
        assert exists(cleaned_dir)
        cleaned_454 = join(cleaned_dir, os.path.basename(fpath_454))
        assert exists(cleaned_454)
        assert 'GGTTCAAGGTTTGAGAAAGGATGGGAAG' not in open(cleaned_454).read()

        do_analysis(project_settings=settings_path, kind='clean_read_stats')
        clean_stats_dir = join(cleaned_dir, 'stats')
        original_stats_dir = join(original_reads_dir, 'stats')
        assert exists(clean_stats_dir)
        assert exists(original_stats_dir)
        assert exists(join(clean_stats_dir, 'global.diff_qual_distrib.png'))
        assert exists(join(clean_stats_dir, 'pl_454.lb_a.general_stats.dat'))


        do_analysis(project_settings=settings_path,
                    kind='prepare_mira_assembly')
        assembly_input = join(project_dir, 'assembly', 'input')
        assert exists(assembly_input)
        mira_in_454 = join(assembly_input, 'backbone_in.454.fasta')
        mira_in_qul = join(assembly_input, 'backbone_in.454.fasta.qual')
        assert exists(mira_in_454)
        assert exists(mira_in_qul)

        do_analysis(project_settings=settings_path, kind='mira_assembly')
        assembly_dir = join(project_dir, 'assembly')
        singular_assembly_dir = sorted(os.listdir(assembly_dir))[0]
        assert exists(join(assembly_dir, singular_assembly_dir,
                                           'stderr.txt'))

        do_analysis(project_settings=settings_path, kind='select_last_assembly')
        assem_result_dir = join(assembly_dir, 'result')
        assert exists(assem_result_dir)
        assert exists(join(assem_result_dir, 'info'))

        #mira does not do any contig
        result_dir = join(assembly_dir, 'result')
        open(join(result_dir, 'contigs.fasta'), 'w')
        open(join(result_dir, 'contigs.qual'), 'w')
        do_analysis(project_settings=settings_path,
                    kind='set_assembly_as_reference')
        os.chdir('/tmp')
        test_dir.close()


    @staticmethod
    def test_mapping_analysis():
        'We can map the reads'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'
        settings_path = create_project(directory=test_dir.name,
                                       name=project_name,
                                      configuration={
                                                    'Snvs':{'min_quality':20},
                                'Sam_processing':{'add_default_qualities':True},
                                })
        project_dir = join(test_dir.name, project_name)
        #setup the original reads
        reads_dir = join(project_dir, 'reads')
        clean_reads_dir = join(reads_dir, 'cleaned')
        os.mkdir(reads_dir)
        os.mkdir(clean_reads_dir)

        solexa =  '@seq1\n'
        solexa += 'TCATTGAAAGTTGAAACTGATAGTAGCAGAGTTTTTTCCTCTGTTTGG\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIIIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa =  '@seq2\n'
        solexa += 'ATATGATTGAAGATATTTCTGGGCTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa =  '@seq14\n'
        solexa += 'ATATGATTGAAGATATTTCTGGGCTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa =  '@seq15\n'
        solexa += 'ATATGATTGAAGATATTTCTGGGCTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa =  '@seq12\n'
        solexa += 'ATATGATTGAAGATATTTCTGGACTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa =  '@seq13\n'
        solexa += 'ATATGATTGAAGATATTTCTGGACTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa =  '@seq16\n'
        solexa += 'ATATGATTGAAGATATTTCTGGACTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'

        sanger  = '>seq3\n'
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

        fpath_sanger = join(clean_reads_dir, 'lb_hola.pl_sanger.fasta')
        fpath_solexa = join(clean_reads_dir,
                                    'lb_hola.pl_illumina.sfastq')
        open(fpath_sanger, 'w').write(sanger)
        open(fpath_solexa, 'w').write(solexa)

        #the reference
        reference_dir = join(project_dir, 'mapping/reference')
        os.makedirs(reference_dir)
        reference_fpath = join(reference_dir, 'reference.fasta')
        out = open(reference_fpath, 'w')
        for line in open(join(DATA_DIR, 'blast/arabidopsis_genes')):
            out.write(line)

        do_analysis(project_settings=settings_path, kind='mapping')
        mapping_dir = join(project_dir, 'mapping')
        singular_mapping_dir = sorted(os.listdir(mapping_dir))[0]
        singular_mapping_dir = join(mapping_dir, singular_mapping_dir)
        assert exists(join(singular_mapping_dir, 'result',
                                    'by_readgroup', 'lb_hola.pl_illumina.bam'))
        do_analysis(project_settings=settings_path, kind='select_last_mapping')
        result_dir = join(mapping_dir, 'result')
        assert exists(result_dir)
        result_dir_by_lib = join(result_dir, 'by_readgroup')
        assert exists(result_dir_by_lib)

        do_analysis(project_settings=settings_path, kind='merge_bam')
        assert exists(join(result_dir, 'merged.bam'))

        annot_input_dir = join(project_dir, 'annotations', 'input')
        os.makedirs(annot_input_dir)
        os.symlink(reference_fpath, join(annot_input_dir, 'reference.fasta'))
        do_analysis(project_settings=settings_path, kind='annotate_snv')
        repr_fpath = join(project_dir, 'annotations', 'repr', 'reference.repr')
        assert "type='snv'" in  open(repr_fpath).read()

        do_analysis(project_settings=settings_path, kind='filter_snvs')
        do_analysis(project_settings=settings_path, kind='write_annotation')
        vcf_fpath = join(project_dir, 'annotations', 'result', 'reference.vcf')
        vcf = open(vcf_fpath).read()
        assert 'vks' in vcf
        assert 'AT5G19860.1' in vcf

        os.chdir('/tmp')
        test_dir.close()

    @staticmethod
    def test_wsg_asembly_analysis():
        'We can assembly with wsg'
        test_dir = NamedTemporaryDir()
        proj_name = 'backbone'
        project_dir = join(test_dir.name, proj_name)
        settings_path = create_project(directory=test_dir.name, name=proj_name)
        #setup the original reads
        reads_dir = join(project_dir, 'reads')
        clean_reads_dir = join(reads_dir, 'cleaned')
        os.mkdir(reads_dir)
        os.mkdir(clean_reads_dir)

        fpath_sanger = join(clean_reads_dir, 'lb_hola.pl_sanger.sfastq')
        open(fpath_sanger, 'w').write(READS_454)
        fpath_sanger3 = join(clean_reads_dir, 'lb_hola2.pl_sanger.fasta')
        open(fpath_sanger3, 'w').write('>seq\nAAAAAAA\n>seq2\nTTTT\n')
        fpath_sanger2 = join(clean_reads_dir, 'lb_caracola.pl_sanger.sfastq')
        open(fpath_sanger2, 'w').write(READS_454)
        do_analysis(project_settings=settings_path, kind='prepare_wsg_assembly')
        result_dir = join(project_dir, 'assembly', 'input')
        assert exists(result_dir)
        frg_fpath = join(result_dir, 'all_seq.frg')
        assert exists(frg_fpath)
        # This dos not work because of the  input data
        #do_analysis(project_settings=settings_path, kind='wsg_assembly')
        os.chdir('/tmp')
        test_dir.close()

    @staticmethod
    def test_cdna_intron_annoation_analysis():
        'We can annotate introns'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'
        blast_db_path = os.path.join(DATA_DIR, 'blast')
        genomic_db = os.path.join(blast_db_path, 'tomato_genome2')
        config = {'Cdna_intron_annotation':{'genomic_db': genomic_db,
                                            'genomic_seqs':genomic_db}}
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
        do_analysis(project_settings=settings_path, kind='annotate_introns')
        repr_fpath = join(project_dir, 'annotations', 'repr', 'seqs.repr')

        assert "type='intron'" in  open(repr_fpath).read()
        os.chdir('/tmp')
        test_dir.close()

if __name__ == "__main__":
    #import sys;sys.argv = ['TestBackbone.test_mapping_analysis']#, 'Test.testName']
    unittest.main()
