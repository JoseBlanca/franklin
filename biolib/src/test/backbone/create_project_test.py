'''
Created on 26/01/2010

@author: jose
'''
import unittest, os.path, shutil

from biolib.utils.misc_utils import NamedTemporaryDir, DATA_DIR
from biolib.backbone.create_project import create_project
from biolib.backbone.analysis import do_analysis, BACKBONE_DIRECTORIES
from configobj import ConfigObj

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

        assert settings_path == os.path.join(test_dir.name,
                                'backbone', BACKBONE_DIRECTORIES['config_file'])
        settings = ConfigObj(settings_path)
        assert settings['General_settings']['project_name'] == 'backbone'
        project_path = os.path.join(test_dir.name, 'backbone')
        assert settings['General_settings']['project_path'] == project_path
        test_dir.close()

    @staticmethod
    def xtest_cleaning_analysis():
        'We can clean the reads'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'
        settings_path = create_project(directory=test_dir.name,
                                       name=project_name)
        project_dir = os.path.join(test_dir.name, project_name)
        #setup the original reads
        reads_dir = os.path.join(project_dir, 'reads')
        original_reads_dir = os.path.join(reads_dir, 'original')
        os.mkdir(reads_dir)
        os.mkdir(original_reads_dir)

        #print original_reads_dir
        fpath_454 = os.path.join(original_reads_dir, 'pt_454.sfastq')
        fpath_ill = os.path.join(original_reads_dir, 'pt_illumina.sfastq')
        open(fpath_454, 'w').write(READS_454)
        open(fpath_ill, 'w').write(READS_ILL)

        do_analysis(project_settings=settings_path, kind='clean_reads')
        cleaned_dir = os.path.join(project_dir, 'reads', 'cleaned')
        assert os.path.exists(cleaned_dir)
        cleaned_454 = os.path.join(cleaned_dir, os.path.basename(fpath_454))
        assert os.path.exists(cleaned_454)

        do_analysis(project_settings=settings_path,
                    kind='prepare_mira_assembly')
        assembly_input = os.path.join(project_dir, 'assembly', 'input')
        assert os.path.exists(assembly_input)
        mira_in_454 = os.path.join(assembly_input, 'backbone_in.454.fasta')
        mira_in_qul = os.path.join(assembly_input, 'backbone_in.454.fasta.qual')
        assert os.path.exists(mira_in_454)
        assert os.path.exists(mira_in_qul)


        do_analysis(project_settings=settings_path, kind='mira_assembly')
        assembly_dir = os.path.join(project_dir, 'assembly')
        singular_assembly_dir = sorted(os.listdir(assembly_dir))[0]
        assert os.path.exists(os.path.join(assembly_dir, singular_assembly_dir,
                                           'stderr.txt'))

        do_analysis(project_settings=settings_path, kind='select_last_assembly')
        assem_result_dir = os.path.join(assembly_dir, 'result')
        assert os.path.exists(assem_result_dir)
        assert os.path.exists(os.path.join(assem_result_dir, 'info'))

        #mira does not do any contig
        result_dir = os.path.join(assembly_dir, 'result')
        open(os.path.join(result_dir, 'contigs.fasta'), 'w')
        open(os.path.join(result_dir, 'contigs.qual'), 'w')
        do_analysis(project_settings=settings_path,
                    kind='set_assembly_as_reference')
        test_dir.close()


    @staticmethod
    def test_mapping_analysis():
        'We can clean the reads'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'
        settings_path = create_project(directory=test_dir.name,
                                       name=project_name)
        project_dir = os.path.join(test_dir.name, project_name)
        #setup the original reads
        reads_dir = os.path.join(project_dir, 'reads')
        clean_reads_dir = os.path.join(reads_dir, 'cleaned')
        os.mkdir(reads_dir)
        os.mkdir(clean_reads_dir)

        solexa =  '@seq1\n'
        solexa += 'TCATTGAAAGTTGAAACTGATAGTAGCAGAGTTTTTTCCTCTGTTTGG\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIIIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa =  '@seq2\n'
        solexa += 'ATATGATTGAAGATATTTCTGGGCTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIIIIUJUAUGJUUJUDFAOUDJOFSUD\n'

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

        fpath_sanger = os.path.join(clean_reads_dir, 'lb_hola.pt_sanger.fastq')
        fpath_solexa = os.path.join(clean_reads_dir,
                                    'lb_hola.pt_illumina.sfastq')
        open(fpath_sanger, 'w').write(sanger)
        open(fpath_solexa, 'w').write(solexa)

        #the reference
        reference_dir = os.path.join(project_dir, 'mapping/reference')
        os.makedirs(reference_dir)
        reference_fpath = os.path.join(reference_dir, 'reference.fasta')
        out = open(reference_fpath, 'w')
        for line in open(os.path.join(DATA_DIR, 'blast/arabidopsis_genes')):
            out.write(line)

        do_analysis(project_settings=settings_path, kind='mapping')

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
