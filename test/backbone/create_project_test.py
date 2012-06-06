'''
Created on 26/01/2010

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

import unittest, os.path, shutil
from os.path import join, exists
from tempfile import NamedTemporaryFile

from franklin.utils.misc_utils import NamedTemporaryDir, TEST_DATA_DIR
from franklin.backbone.create_project import (create_project,
                                              create_configuration)
from franklin.backbone.analysis import (BACKBONE_DIRECTORIES,
                                        BACKBONE_BASENAMES,
                                        scrape_info_from_fname)
from franklin.backbone.backbone_runner import do_analysis
from franklin.seq.readers import seqs_in_file

READS_NOQUAL = '''>FM195262.1
TACGGCCGGGGTNNCNNANNNNGCATTCTCGCAGGGTCTTTCTACACTATTAGATAAGAT
GGATCCTTCTCAGAGAGTGAAGTTTGTTCAGGAAGTCAAGAAGGTTCTTGGATGATGATA
TGATACCAACACATCCAACACAATATGCGCATGCTACATGTTATTTTTCAAGTACATACA
TAGAAGGATATTGCTTGGCCTTGATTGATCATGTCTGATCTAAGTCGATCATTATTTTCT
TGAAACTTCCTTTCGGACGTGGTGCTATGGTTGATGAATTTGGATGTGTGCGTTCTGCCA
GGTGTAAGCCCAAAGGTTTATACAGACCGAGTTAAGGTTAGGAAGAGCACGAGTGAACTT
>FM187038.1
GCACGAGGGCCGTCGCAGAGGCGGTGCTGAGGTCGAGGGCCGGTCTTGGCAGGCCACAAC
AGCCCACTGGCTCGTTCCTCTTCCTGGGTCCGACTGGCGTGGGGAAAACTGAGCTGGCCA
AGGCCCTAGCCGAACAGCTGTTCGACGACGAGAACCTTCTTGTCCGCATCGACATGTCGG
AGGTCGCTCGCGATCTGGTCATGCAGGAGGTGAGGAGGCACTTCCGCCCTGAGCTGCTGA
ACCGTCTCGACGAGATCGTGATCTTCGATCCTCTGTCCCACGAGCAGCTGAGGAAGGTCG
CTCGCCTTCAGATGAAGGATGTGGCCGTCCGTCTTGCCGAANNNNNCATCGCTCTGGCTG
TGACCGANNNNNCATTGGACATCATCTTGTCTCTCTCTNNNNNNTCNNNNT
'''
READS_SOLID = '''@10_1824_570_F3
NAGATACGTTACGGTCATCCGGGCCAGCTCATCGCNGAAANGNAAGTCA
+
!@A=BB+@:)@9@.<B=67/B=>;@B><@<B)75?!:@=2!7!::>/?;
@10_1993_178_F3
NGCACATACTGGATAGGCCATTCAGAGGAGACAGACAAGGTCAAGGACA
+
!9,=:===5;-;??<A=;4>@?=@>(*(?@?8;6=<8?<;7%??<989@
@10_212_402_F3
NCTGTTAGGACGTTGTCCGACTAGCGTTCTGTCAGACTAACATAACCAA
+
!;?;<9B:?>,:AB;B;>7B=BA??ABB:BA?7?4?<9':>'@8<B%2:
@10_51_1274_F3
NGATCCCGACTCGGAATCTAGAACTCGCCTAGGCTCTCACGACGCTAGG
+
!@??7:9?<<:??*=@B=<;@?@<B?=>?@<>=?AAA9:6<;:54B;@4
@10_577_337_F3
NCACGGACCTCAGTTAACCTAGCTAGGCGTCTGCGATGGTTACTTTCAT
+
!=8%@5?ABB&BBBB@BB7B@<BB<:BA-A:B?;@<?>@B4@?+,6@B8
@10_796_926_F3
NAAGAAGTCGGTCAGCTTGTAATTACACGAGCACCGTGCGAACTGCAAG
+
!9@B8?5?A67B<@>9'B@/=4B;5+@77>;15:/A;>=%>:2):1@><


'''
READS_454 = '''@FKU4KFK07H6D2L
GGTTCAAGGTTTGAGAAAGGATGGGAAGAAGCCAAATGCCTACATTGCTGATACCACTACGGCAAATGCTCAAGTTCGGACGCTTGCTGAGACGGTGAGACTGGATGCAAGAACTAAGTTATTGAATAGTCAGCATGCATGATTAGGCTAAGCCGTAAGCATAGCATGACCCCATTGGCAAAGCTAGCATGATACGACATCATTATAGCGAGAGACGCATATCGAGAATGAGCGATCAGCACATGTCAGCGAGCTACTGACTATCATATATAGCGCAGAGACGACTAGCATCGAT
+
CCFFFFFCC;;;FFF99;HHIHECCHHIIIFHHEEEHHHHHHIIIFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFGGGGFFGGBCCGFFFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
@FKU4KFK07IJEB9
AATTCCCTTTCCGGTGATCGCCTCTCCTAGCAATGATGCAGCAAAGCCGATCATAGCAACACGGCCGACGAAGAGTTCATTCTGCTTGGTAAAACCAATGCCACCAGAAGTCCCAAAAATTCCATCTTCAACCTTCTGCTTCGGCTTCTCAACCTTCTTGGGCGGGGCTTTAGTTCTGGATTTGAACACAGCAAGAGGGTGAAACCA
+
FFFFFFFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFGGGGFFFFGGGGGGFFFFFFGAA@5557>55577G
@FKU4KFK07IMD9G
ACCCTTCTAGCAGTAGGAATGACTTGACCACCTCCTCTGTGGATTGTAGTAGTAGATGATGGCATCAGCGTGAAGCACCACATCACAAACTTCAAAACAGATGCCT
+
B>>>GGFFFFFFFFFGGGGHHHHHHHHHIIIIIIIIIIIIIIIIIFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFGCB;;;;BGFFFFFFF
@FKU4KFK07IE66H
ATGAACCGCAAAAGCTTGTATGCTGTATTGCCTTGATTTGGTTTCCAAGATTCTTCCCACATATATGATGATGATGATGATGTTAGGAGAGAGTGTAGCTGGAAACCGAGACTTTGGCTGGAGAAGCATGCGTATGATCATGACAGATAGC
+
FFFFFFFFG7777CGFFFFHHHHHHIHHHHHHHHHHHHHHH==?CCFFFFFFFGGCCCGFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFGGGFFFFFFFFFFFFFGGGFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
'''

SANGER_QUAL = """@SEX_N_int_8782-8882_agen_plate-001_A01_F.b.ab1 CHROMAT_FILE: SEX_N_int_8782-8882_agen_plate-001_A01_F.b.ab1 PHD_FILE: SEX_N_int_8782-8882_agen_plate-001_A01_F.b.ab1.phd.1 CHEM: term DYE: big TIME: Mon May 24 13:19:19 2010
gacctatagaacaagtcggtccggaattcccgggatcatcaacatgccttccaccgctgcttccaccaaggtcccccagaccaccatgaacttcaacggctactgcgttgttatgtaaataacctcacagtgccttgaagggccaacgagagcattgagaggcagatgccaaccagctttcctccagcagatcttctggcaaaccacccatctcacctgagtgatcctgttctcgcaacatcactcgatgagaaccactgaagatcggtatctcgcctctcaacgtcagaccagcaatctgcttttccctctatcatgttcgagaagcccacctctcgttccgccacaacaacaacctcgcctaggaagcttcttgtcactgcttcatggtggactctgcacaatgggagcacccactgacaacaaccggagcggtgaggatctttccttgatgggaatggaaaaaaacctttgtcaaccggggatatacatccgcacacatcacacacaccatcctttcacataaacaactttcgcttctttttgttggggataacttcacacgggactacttttcgtccacacggagctttccctctttcaattttcgacactcaggcaactaaccgtatcttttggaccctttctctctttcgttgttcgttcgacaccaaacaccctttgcgggctgaggacatcgcctgcgggaaacgatgcaacaacttctcttcggttcacactcttctgcaaacgtcacctgggagtcggtggactcacacataccaacctcaacctcttccctcggctgggagggaggctatttcgggaatgcttctcacttttggttgggggggacttgaccgggcttggggattgggcctttatcccccccttcctttttttttgcgaacctttccttggaaaatttctttgggaaagggaattttaaaagggggggggaaatccttttccttttttacttgttctttttgagaga
+
)))))..)))))*.16FFFFFH412-/4\HJK@@?B=BBBF??:A98==1+++13DJJJJ5+++)0.5==AHHHH@9555<G=74--0;@JA?B==89888AEFMB>>BADDDDDKKBBDDDLGGFA:200335A=AA566=====A=M>=AMMDEJ>>>LLL<<<LAC555G7000C<?//-6--++/-228LLDDBL77/35558<DL@@@DHLLL>;..5ALLLLDDD10.3.//>9D6000,,0HH8//HHDDHLLMH;-0.<<<<<<<=</<<>>><D>>>D888H<...<<MLHHLL>>>>>>HLLLL<=//..,++++0088LLH>>>>>@@DHH//.88>464844<LHH228HL11/H;<<88688HH222??DDD<;;LDH@@@DD@9822388<DDD89345CC223995963-./33<8;;D;;;LHHLKDDDLJDDDHLLLDBBHHHHKHL><<885131--00544AADFLLA=9=.-.,33A<<555<<<A=33223373330599<==J<<833..-9050C9888EAA?F?AAFFFFA@834303722376;8:88=9A88,-,-69:---4.79FCHFFCB//1355C44..36-.0007=3/-0,,,01.3.66?AAADDDEBFA==BBDFBF==552.1.,+-/26:=7-2++)*.1113139<@320;/-.44640+))))****+*****,,0.,**+*/58888855620/----70)))))+-44<25++*****+02++*+,++,,+++-))))),,**))+')'')'*000-,))))))+++('+*+*))))))++++))))),,,++)))))+.2-,,,+***,***+,))))))*-,-+)))))))))*+,,)))++*,10*))))),-.1?88//-++++++++/+****))))))+-**+*.)**)+*+*+*++)'''++-+010-*)))))())))+,00/++632++*+***+**,..-,,,,,,
@SEX_N_int_8782-8882_agen_plate-001_A01_R.g.ab1 CHROMAT_FILE: SEX_N_int_8782-8882_agen_plate-001_A01_R.g.ab1 PHD_FILE: SEX_N_int_8782-8882_agen_plate-001_A01_R.g.ab1.phd.1 CHEM: term DYE: big TIME: Mon May 24 13:19:19 2010
gggggggcggtattctacatcaagtccccactcgtggcacctcgcctcaggaacactttgaccgaccgaactctccgcccgactgcctctgtggccggtgcttcgaatactcatcgacgaacatctccctatttgccactccgggttctaccttgccgagcgttcgcacaccatgtccacgggataacccctatgactggatcgcagagaacatgctggctaggaggtagtccgcgcccgtgaccatggtctctgcaccataaagcgagctctcaatatcttgcagacttcttcccccttcggttactactttacccgtcccaggtcgtctggtgtccagttcctttcaatagacctccagcgcccatcacatatggcgtctgtgcgtataccgcactaagttcccagaagttgtctctacgtggatgattgtgagtgatacgatgtccccggcacccactaccgaactcatggttgtatcagatcctgacgaggagtacgtggaggataggatcgtggcgctactgaactcaaggaaatctccatctgctgatccaagaattctaaaagatcacgacgccggaatcgcaccctagtccctgcagcttcccccgaacaatgcttagagagtattcgatgtttccaacagttatgctttgggccggatcaagggggcccatattccagaacaagctttggtcgcgctccatgcaccggctggacagcctctttcccataagctcggggtggttccacttcctcccgcaccttctttccgccttgggccccggtccggaaaaaaaaggaaaaccggaaagggacaaaagggtgcgaggggaaataccccccattcttctcttttcccggaataagaagaccaaaagggggggggggggccctttgccc
+
47771/***).//.4--,--03/..000...0...-3,--9.,,-,,,.,,.9......-1-,,0---04.,,,,,-;<AA94--....33,,,--../??=90...A<9//.8-,,-13---.0?9<00...044;900////668=555555333000006-,,,)))004,,0879-0-..-,,--++,.,103336630//--<-.-0,,--/22.3A------.-<000--------00000-,,+-**++-,---,,----,--633//8;02//-<...33888A33...0-3.--,--,,,-,--7=DC85540.-84-/+++++--,66--,--399:73251,,0000,,,,,6<A?==??88,.-,))))-766--+))-+++,,-++))-127-.)))))*030--+.+,.;000-,,000,,,44000-,,,-00032--,,,0-01=00+++0122///0001--*)))++423:,,+++//++**))++1153-3,,,0,,,,,+++++0-1-,--410**))))++*)))+)+***+-*++++++++..0004-.++++1,,,,,)+++,-***))))))*,'''')),+,+++-****+))))++)))))'*+++,,)*)))))+)+,)))))))))))))+)**)))))))))))-)))))))+,,*+*****))****---+'''))/,,+++))))))))++))))))+,++-1,))'((''')'(''''..+)''))))'')',,,,,+++.,++*'''''(205/0-,,,-00119-.--+(()))((')))*0..,,**+((*(''++****))+)'('))..900,/,))'')((.0+),,,*))''''''''))))))*+.40.))))))*****0***++
"""

READS_ILL = '''@HWI-EAS59:3:1:5:1186#0/1
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
THREADS = False

class TestBackboneUtils(unittest.TestCase):
    'It tests the backbone utils'


    def test_scrape_info_from_fname(self):
        'it tests scrape_info_from_fname'
        fhand = NamedTemporaryFile(prefix='pl_illumina.sm_test.lb_lib1.',
                                   suffix='sfastq')
        fhand.write(READS_ILL)
        fhand.flush()
        file_info = scrape_info_from_fname(fhand.name)
        assert file_info['lb'] == 'lib1'
        assert file_info['format'] == 'fastq'

        # this should fail
        fhand = NamedTemporaryFile(prefix='pl__illumina.sm_test.lb_lib1.',
                                   suffix='sfastq')
        fhand.write(READS_ILL)
        fhand.flush()
        try:
            file_info = scrape_info_from_fname(fhand.name)
            self.fail()
        except RuntimeError:
            pass


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
        settings = create_configuration(settings_path)
        assert settings['General_settings']['project_name'] == 'backbone'
        project_path = join(test_dir.name, 'backbone')
        assert settings['General_settings']['project_path'] == project_path
        assert settings['Cleaning']['strip_n_percent'] == 2.0
        content = open(settings_path).read()
        assert 'strip_n_percent' in content
        test_dir.close()

    @staticmethod
    def test_cleaning_analysis_lucy():
        'We can clean the reads'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'
        univec = os.path.join(TEST_DATA_DIR, 'blast', 'univec')
        configuration = {'Cleaning':{'vector_database':None}}
        settings_path = create_project(directory=test_dir.name,
                                       name=project_name,
                                       configuration=configuration)

        project_dir = join(test_dir.name, project_name)
        #setup the original reads
        reads_dir = join(project_dir, 'reads')
        original_reads_dir = join(reads_dir, 'raw')
        os.mkdir(reads_dir)
        os.mkdir(original_reads_dir)

        os.makedirs(join(project_dir, 'config_data', 'lucy'))
        lucy_settings = join(project_dir, 'config_data', 'lucy', 'lucy.conf')
        luc_c = open(lucy_settings, 'w')
        luc_c.write(repr({'ps':{'vector_file':'tmp' , 'splice_file':'tmp'}}))
        luc_c.flush()

        #print original_reads_dir
        fpath_noqual = join(original_reads_dir, 'pl_sanger.lb_ps.fasta')
        fpath_qual = join(original_reads_dir, 'pl_sanger.lb_andreas.sfastq')
        fpath_454 = join(original_reads_dir, 'pl_454.lb_ps.sfastq')
        fpath_ill = join(original_reads_dir, 'pl_illumina.lb_psi.sfastq')
        open(fpath_noqual, 'w').write(READS_NOQUAL)
        open(fpath_qual, 'w').write(SANGER_QUAL)
        open(fpath_454, 'w').write(READS_454)
        open(fpath_ill, 'w').write(READS_ILL)
        do_analysis(project_settings=settings_path, kind='clean_reads',
                    silent=True)
        cleaned_dir = join(project_dir, 'reads', 'cleaned')
        assert exists(cleaned_dir)

        cleaned_qual = join(cleaned_dir, os.path.basename(fpath_qual))
        assert 'SEX' in open(cleaned_qual).read()

        cleaned_454 = join(cleaned_dir, os.path.basename(fpath_454))
        assert exists(cleaned_454)

        cleaned_noqual = join(cleaned_dir, os.path.basename(fpath_noqual))
        clean_seqs = open(cleaned_noqual).read()
        assert clean_seqs.startswith('>FM195262.1\nGCATTCTCG')

    @staticmethod
    def test_cleaning_analysis():
        'We can clean the reads'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'
        project_dir = join(test_dir.name, project_name)
        adaptors_dir = join(project_dir, 'config_data', 'adaptors')
        adaptors_path_454 = join(adaptors_dir, '454_adaptors')
        words = ['^ATGAAC', 'TTGATTTGGT']
        univec = os.path.join(TEST_DATA_DIR, 'blast', 'univec+')
        configuration = {'Cleaning':{'vector_database':univec,
                                     'adaptors_file_454':adaptors_path_454,
                                     'short_adaptors_454':words,
                                     'edge_removal':{'454_left':3,
                                                     '454_right':3}},
                         'General_settings':{'threads':THREADS}}

        settings_path = create_project(directory=test_dir.name,
                                       name=project_name,
                                       configuration=configuration)

        #setup the original reads
        reads_dir = join(project_dir, 'reads')
        original_reads_dir = join(reads_dir, 'raw')
        os.mkdir(reads_dir)
        os.mkdir(original_reads_dir)

        os.makedirs(adaptors_dir)
        adap_fhand = open(adaptors_path_454, 'w')
        adap_fhand.write('''>smart_5_cds_primer_1
GGTTCAAGGTTTGAGAAAGGATGGGAAG\n''')
        adap_fhand.close()

        #print original_reads_dir
        fpath_454 = join(original_reads_dir, 'pl_454.lb_a.sfastq')
        fpath_ill = join(original_reads_dir, 'pl_illumina.lb_b.sfastq')
        fpath_solid = join(original_reads_dir, 'pl_solid.lb_prueba.sfastq')

        open(fpath_solid, 'w').write(READS_SOLID)
        open(fpath_454, 'w').write(READS_454)
        open(fpath_ill, 'w').write(READS_ILL)

        do_analysis(project_settings=settings_path, kind='clean_reads',
                    silent=True)
        cleaned_dir = join(project_dir, 'reads', 'cleaned')
        assert exists(cleaned_dir)
        cleaned_454 = join(cleaned_dir, os.path.basename(fpath_454))
        assert exists(cleaned_454)
        seqs = list(seqs_in_file(open(cleaned_454)))
        # It means thar the adaptor has been removed
        seq = seqs[0].seq
        assert 'GGTTCAAGGTTTGAGAAAGGATGGGAAG' not in seq

        seq = seqs[2].seq
        # It means that the starting word has been removed
        assert  seq.startswith('TTCCAAGATTCTTCCCACAT')

        # solid
        cleaned_solid = join(cleaned_dir, os.path.basename(fpath_solid))
        clean_seqs = open(cleaned_solid).read()
        assert '10_1824_570_F3' not in clean_seqs


        do_analysis(project_settings=settings_path,
                    kind='prepare_mira_assembly', silent=True)
        assembly_input = join(project_dir, 'assembly', 'input')
        assert exists(assembly_input)
        mira_in_454 = join(assembly_input, 'backbone_in.454.fasta')
        mira_in_qul = join(assembly_input, 'backbone_in.454.fasta.qual')
        assert exists(mira_in_454)
        assert exists(mira_in_qul)

        do_analysis(project_settings=settings_path, kind='mira_assembly',
                    silent=True)
        assembly_dir = join(project_dir, 'assembly')
        singular_assembly_dir = sorted(os.listdir(assembly_dir))[0]
        test_dir.close()

    @staticmethod
    def test_read_stats_analysis():
        'It test the read statistics'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'
        project_dir = join(test_dir.name, project_name)

        settings_path = create_project(directory=test_dir.name,
                                       name=project_name)

        #setup the original reads
        reads_dir = join(project_dir, 'reads')
        original_reads_dir = join(reads_dir, 'raw')
        os.mkdir(reads_dir)
        os.mkdir(original_reads_dir)
        fpath_454 = join(original_reads_dir, 'pl_454.lb_a.sfastq')
        fpath_ill = join(original_reads_dir, 'pl_illumina.lb_b.sfastq')
        open(fpath_454, 'w').write(READS_454)
        open(fpath_ill, 'w').write(READS_ILL)

        #the cleaned reads
        cleaned_reads_dir = join(reads_dir, 'cleaned')
        os.mkdir(cleaned_reads_dir)
        fpath_454 = join(cleaned_reads_dir, 'pl_454.lb_a.sfastq')
        fpath_ill = join(cleaned_reads_dir, 'pl_illumina.lb_no_raw.sfastq')
        open(fpath_454, 'w').write(READS_454)
        open(fpath_ill, 'w').write(READS_ILL)

        do_analysis(project_settings=settings_path, kind='read_stats',
                    silent=True)

        clean_stats_dir = join(cleaned_reads_dir, 'stats')
        clean_fnames = os.listdir(clean_stats_dir)
        expected_fnames = ['pl_454.lb_a.qual.diff',
                           'pl_illumina.lb_no_raw.qual',
                           'pl_454.lb_a.qual',
                           'pl_illumina.lb_no_raw.length',
                           'pl_454.lb_a.length',
                           'pl_454.lb_a.length.diff']
        for fname in expected_fnames:
            assert fname + '.dat' in clean_fnames
            assert fname + '.svg' in clean_fnames

        statistics_fpath = join(clean_stats_dir,
                                BACKBONE_BASENAMES['statistics_file'])
        content = open(statistics_fpath).read()
        assert content == '''statistics for pl_illumina.lb_no_raw.sfastq
-------------------------------------------
Num sequences: 6
Total sequence length: 172
Sequence length minimum: 24
Sequence length maximum: 31
Sequence length average: 28.67
Sequence length variance: 10.89
Sequence qualities minimum: 4
Sequence qualities maximum: 34
Sequence qualities average: 29.63
Sequence qualities variance: 47.80

statistics for pl_454.lb_a.sfastq
---------------------------------
Num sequences: 4
Total sequence length: 759
Sequence length minimum: 106
Sequence length maximum: 295
Sequence length average: 189.75
Sequence length variance: 4972.69
Sequence qualities minimum: 20
Sequence qualities maximum: 40
Sequence qualities average: 36.99
Sequence qualities variance: 8.19

'''

        boxplot_fpath = join(clean_stats_dir,
                             'pl_illumina.lb_no_raw' + '.qual.boxplot.dat')
        exp = 'distrib\tmean\tstd_deviation\t1st_quartile\tmedian\t3rd_qualtile'
        assert exp in open(boxplot_fpath).read()
        freq_nucl_fpath = join(clean_stats_dir, 'pl_454.lb_a.freq_position.svg')
        nucl_freq = open(freq_nucl_fpath).read()
        assert 'style="fill:#0000ff;stroke:#000000;"/>' in nucl_freq

    @staticmethod
    def test_read_stats_analysis2():

        # another read stats with real data

        clean_fpath = os.path.join(TEST_DATA_DIR, 'clean_stats', 'cleaned',
                                   'lb_sflp2.pl_sanger.sm_t111.sfastq')
        raw_fpath = os.path.join(TEST_DATA_DIR, 'clean_stats', 'raw',
                                   'lb_sflp2.pl_sanger.sm_t111.sfastq')

        test_dir = NamedTemporaryDir()
        project_name = 'backbone'
        project_dir = join(test_dir.name, project_name)

        settings_path = create_project(directory=test_dir.name,
                                       name=project_name)

        #setup the original reads
        reads_dir = join(project_dir, 'reads')
        original_reads_dir = join(reads_dir, 'raw')
        cleaned_reads_dir = join(reads_dir, 'cleaned')
        os.mkdir(reads_dir)
        os.mkdir(original_reads_dir)
        os.mkdir(cleaned_reads_dir)

        shutil.copy(clean_fpath, cleaned_reads_dir)
        shutil.copy(raw_fpath, original_reads_dir)
        do_analysis(project_settings=settings_path, kind='read_stats',
                    silent=True)


    @staticmethod
    def test_remove_output_on_error():
        'We remove files when we have an error on cleaning'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'
        project_dir = join(test_dir.name, project_name)

        configuration = {'Cleaning':{'adaptors_file_454':'AKHSGDASD'},
                         'General_settings':{'threads':THREADS}}


        settings_path = create_project(directory=test_dir.name,
                                       name=project_name,
                                       configuration=configuration)

        #setup the original reads
        reads_dir = join(project_dir, 'reads')
        original_reads_dir = join(reads_dir, 'raw')
        os.mkdir(reads_dir)
        os.mkdir(original_reads_dir)

        #fake solid reads
        fpath_solid = join(original_reads_dir, 'pl_454.lbx_b.sfastq')

        try:
            do_analysis(project_settings=settings_path, kind='clean_reads',
                        silent=True)
        except KeyError:
            pass
        output__fpath = join(reads_dir, 'cleaned', 'pl_454.lb_b.sfastq')
        assert not exists(output__fpath)

class ConfigurationTest(unittest.TestCase):
    'Tests for configirations'
    @staticmethod
    def test_create_configuration():
        'Test default config'
        fhand = NamedTemporaryFile()
        fhand.write('[Snvs]\nmin_num_alleles = 2\n[[edge_removal]]\n454_left=None')
        fhand.flush()
        config_fpath = fhand.name
        config = create_configuration(config_fpath)
        assert config['Snvs']['min_quality'] == 45

        #malformed
        fhand = NamedTemporaryFile()
        fhand.write('[blast]\n[[nr]]\nkind="nucl"\npath="path"\nspecies="all"')
        fhand.flush()
        config_fpath = fhand.name
        config = create_configuration(config_fpath)
        assert config['Snvs']['min_quality'] == 45

class UtilTest(unittest.TestCase):
    'test utils from backbone'

    @staticmethod
    def test_scrape_info_from_fname():
        'scrape info from fpath'
        fhand = NamedTemporaryFile(prefix='st_prot.A.', suffix='.fasta')
        fhand.write('>seq\nTGATGC')
        fhand.flush()
        info = scrape_info_from_fname(fhand.name)
        assert info['st'] == 'prot'

if __name__ == "__main__":
#    import sys;sys.argv = ['', 'TestBackboneUtils.test_scrape_info_from_fname']
    unittest.main()
