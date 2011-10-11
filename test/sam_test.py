'''
It test sam module functions
Created on 05/01/2010

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

from tempfile import NamedTemporaryFile
import unittest, os

from franklin.utils.misc_utils import TEST_DATA_DIR
from StringIO import StringIO

from franklin.sam import (bam2sam, sam2bam, merge_sam, bamsam_converter,
                          add_header_and_tags_to_sam, sort_bam_sam,
                          standardize_sam, realign_bam,
                          bam_distribs, sample_bam, bam_general_stats,
                          remove_unmapped_reads, get_read_group_info)
import pysam

THREADS=4


class SamTest(unittest.TestCase):
    'It test sam tools related functions'

    @staticmethod
    def testbam2sam():
        'It test bam2sam function'
        bampath = os.path.join(TEST_DATA_DIR, 'seq.bam')
        sampath = NamedTemporaryFile(suffix='.sam').name
        bam2sam(bampath, sampath, header=True)
        assert 'SN:SGN-U572743' in open(sampath).readline()

    @staticmethod
    def testsam2bam():
        'It test sam2bam function'
        bampath = os.path.join(TEST_DATA_DIR, 'seq.bam')
        sampath = NamedTemporaryFile(suffix='.sam').name
        bam2sam(bampath, sampath, header=True)
        assert 'SN:SGN-U572743' in open(sampath).readline()

        newbam = NamedTemporaryFile(suffix='.bam')
        sam2bam(sampath, newbam.name)
        newsam = NamedTemporaryFile(suffix='.sam')
        bam2sam(newbam.name, newsam.name, header=True)
        newsam_content = open(newsam.name).read()
        oldsam_content = open(sampath).read()

        assert newsam_content == oldsam_content

    @staticmethod
    def test_format_converter():
        'It test BAM SAM converter'
        #bam to sam
        bampath = os.path.join(TEST_DATA_DIR, 'seq.bam')
        sam = NamedTemporaryFile(suffix='.sam')
        bamsam_converter(bampath, sam.name)
        sam.flush()
        assert 'SN:SGN-U572743' in sam.read()

        #SAM TO bam
        bam = NamedTemporaryFile(suffix='.bam')
        bamsam_converter(sam.name, bam.name)
        bam.flush()

    @staticmethod
    def test_merge_sam():
        'It merges two sams'
        reference = NamedTemporaryFile(suffix='.sam')
        reference.write('''>SGN-U572743
atatata
>SGN-U576692
gcgc''')
        sam1 = NamedTemporaryFile(suffix='.sam')
        sam1.write('''@SQ	SN:SGN-U576692	LN:1714
@SQ	SN:SGN-U572743	LN:833
@RG	ID:g1	LB:g1	SM:g1
@RG	ID:g2	LB:g2	SM:g2
SGN-E221403	0	SGN-U576692	1416	207	168M	*	0	0	AGCCTGATAAAGGTCTGCCTACGTGTTTTAAGTGGAATCCGTTTCCCCATGTCCAAACCTTCTAAATAGTTTTTTGTGTTAGTTCTTGTATGCCACATACAAAAATTAACAAACTCTTTTGCCACATATGTTCCAGCACGTCAAAGCAACATGTATTTGAGCTACTTT	558<///035EB@;550300094>>FBF>>88>BBB200>@FFMMMJJ@@755225889>0..14444::FMF@@764444448@;;84444<//,4,.,<<QFBB;::/,,,.69FBB>9:2/.409;@@>88.7,//55;BDK@11,,093777777884241<:7	AS:i:160	XS:i:0	XF:i:3	XE:i:4	XN:i:0	RG:Z:g2
SGN-E221664	0	SGN-U572743	317	226	254M24S	*	0	0	GGATGATCTTAGAGCTGCCATTCAAAAGATGTTAGACACTCCTGGGCCATACTTGTTGGATGTGATTGTACCTCATCAGGAGCATGTTCTACCGATGATTCCCAGTGGCGGTGCTTTCAAAAATGTGATTACGGAGGGTGATGGGAGACGTTCCTATTGACTTTGAGAAGCTACATAACTAGTTCAAGGCATTGTATTATCTAAAATAAACTTAATATTTATGTTTACTTAAAAGTTTTTCATTGTGTGAAGGAAAAAAAAAAAAAAAAAAAAAAAAA	999@7<22-2***-,206433>:?9<,,,66:>00066=??EEAAA?B200002<<@@@=DB99777864:..0::@833099???<@488>></...<:B<<.,,8881288@BBDDBD885@@;;9:/9.,,,99B99233885558=?DKKKDDAA??DKBB=440/0<8?DEDFBB??6@152@@FBMFIIDDDDDDKKKOK@@@@DD:N688BBDDDBBBKKDEDDBN977?<9<111:<??==BKMPKKBB==99>QQYYYYYYYYYYYYQQ	AS:i:250	XS:i:0	XF:i:0	XE:i:7	XN:i:0	RG:Z:g1
''')
        sam1.flush()
        sam2 = NamedTemporaryFile(suffix='.sam')
        sam2.write('''@SQ	SN:SGN-U576692	LN:1714
@SQ	SN:SGN-U572743	LN:833
@RG	ID:g1	LB:g1	SM:g1
@RG	ID:g3	LB:g3	SM:g3
SGN-E200000	0	SGN-U572743	317	226	254M24S	*	0	0	GGATGATCTTAGAGKTGCCATTCAAAAGATGTTAGACACTCCTGGGCCATACTTGTTGGATGTGATTGTACCTCATCAGGAGCATGTTCTACCGATGATTCCCAGTGGCGGTGCTTTCAAAAATGTGATTACGGAGGGTGATGGGAGACGTTCCTATTGACTTTGAGAAGCTACATAACTAGTTCAAGGCATTGTATTATCTAAAATAAACTTAATATTTATGTTTACTTAAAAGTTTTTCATTGTGTGAAGGAAAAAAAAAAAAAAAAAAAAAAAAA	999@7<22-2***-,206433>:?9<,,,66:>00066=??EEAAA?B200002<<@@@=DB99777864:..0::@833099???<@488>></...<:B<<.,,8881288@BBDDBD885@@;;9:/9.,,,99B99233885558=?DKKKDDAA??DKBB=440/0<8?DEDFBB??6@152@@FBMFIIDDDDDDKKKOK@@@@DD:N688BBDDDBBBKKDEDDBN977?<9<111:<??==BKMPKKBB==99>QQYYYYYYYYYYYYQQ	AS:i:250	XS:i:0	XF:i:0	XE:i:7	XN:i:0	RG:Z:g1
SGN-E40000	0	SGN-U576692	1416	207	168M	*	0	0	AGCCTGATAAAGGTCTGCCTACGTGTTTTAAGTGGAATCCGTTTCCCCATGTCCAAACCTTCTAAATAGTTTTTTGTGTTAGTTCTTGTATGCCACATACAAAAATTAACAAACTCTTTTGCCACATATGTTCCAGCACGTCAAAGCAACATGTATTTGAGCTACTTT	558<///035EB@;550300094>>FBF>>88>BBB200>@FFMMMJJ@@755225889>0..14444::FMF@@764444448@;;84444<//,4,.,<<QFBB;::/,,,.69FBB>9:2/.409;@@>88.7,//55;BDK@11,,093777777884241<:7	AS:i:160	XS:i:0	XF:i:3	XE:i:4	XN:i:0	RG:Z:g3
''')
        sam2.flush()
        sam3	= NamedTemporaryFile(suffix='.sam')
        merge_sam(infiles=[sam1,	sam2],	outfile=sam3,	reference=reference)
        sam3.seek(0)
        sam3_content = sam3.read()

        assert	'SN:SGN-U572743'	in	sam3_content
        assert	'SGN-E200000'	in	sam3_content
        assert	'SGN-E221664'	in	sam3_content

        #open('/home/peio/tmp/samtools/merged.sam',	'w').write(sam3.getvalue())

    @staticmethod
    def	test_readgroup_to_sam():
        'It test that we can add the readgroup the and the header to a bam'
        sam_value	=	'''@SQ	SN:SGN-U576692	LN:1714
@SQ	SN:SGN-U572743	LN:833
SGN-E221403	0	SGN-U576692	1416	207	168M	*	0	0	AGCCTGATAAAGGTCTGCCTACGTGTTTTAAGTGGAATCCGTTTCCCCATGTCCAAACCTTCTAAATAGTTTTTTGTGTTAGTTCTTGTATGCCACATACAAAAATTAACAAACTCTTTTGCCACATATGTTCCAGCACGTCAAAGCAACATGTATTTGAGCTACTTT	558<///035EB@;550300094>>FBF>>88>BBB200>@FFMMMJJ@@755225889>0..14444::FMF@@764444448@;;84444<//,4,.,<<QFBB;::/,,,.69FBB>9:2/.409;@@>88.7,//55;BDK@11,,093777777884241<:7	AS:i:160	XS:i:0	XF:i:3	XE:i:4	XN:i:0
SGN-E221664	0	SGN-U572743	317	226	254M24S	*	0	0	GGATGATCTTAGAGCTGCCATTCAAAAGATGTTAGACACTCCTGGGCCATACTTGTTGGATGTGATTGTACCTCATCAGGAGCATGTTCTACCGATGATTCCCAGTGGCGGTGCTTTCAAAAATGTGATTACGGAGGGTGATGGGAGACGTTCCTATTGACTTTGAGAAGCTACATAACTAGTTCAAGGCATTGTATTATCTAAAATAAACTTAATATTTATGTTTACTTAAAAGTTTTTCATTGTGTGAAGGAAAAAAAAAAAAAAAAAAAAAAAAA	999@7<22-2***-,206433>:?9<,,,66:>00066=??EEAAA?B200002<<@@@=DB99777864:..0::@833099???<@488>></...<:B<<.,,8881288@BBDDBD885@@;;9:/9.,,,99B99233885558=?DKKKDDAA??DKBB=440/0<8?DEDFBB??6@152@@FBMFIIDDDDDDKKKOK@@@@DD:N688BBDDDBBBKKDEDDBN977?<9<111:<??==BKMPKKBB==99>QQYYYYYYYYYYYYQQ	AS:i:250	XS:i:0	XF:i:0	XE:i:7	XN:i:0
'''

        insam = NamedTemporaryFile(prefix='lb_group1.', suffix='.sam')
        insam.write(sam_value)

        outsam = NamedTemporaryFile()
        add_header_and_tags_to_sam(insam,	outsam)

        out_content	=	open(outsam.name).read()
        assert 'RG:Z:group1' in out_content
        assert 'SM:group1' in out_content

        insam = NamedTemporaryFile(prefix='sm_sample1.lb_group1.', suffix='.sam')
        insam.write(sam_value)

        outsam = NamedTemporaryFile()
        add_header_and_tags_to_sam(insam,    outsam)

        out_content = open(outsam.name).read()
        #print out_content
        assert 'RG:Z:sample1_group1' in out_content
        assert 'SM:sample1' in out_content

    @staticmethod
    def test_sort_bam():
        'It test that we can sort bams using picard'
        #sort sam
        sam_fpath = os.path.join(TEST_DATA_DIR, 'samtools', 'seqs.sam')
        sorted_samfhand = NamedTemporaryFile(suffix='.sam')
        sort_bam_sam(sam_fpath, sorted_samfhand.name)

        #sort bam
        bam_fpath = os.path.join(TEST_DATA_DIR, 'samtools', 'seqs.bam')
        sorted_bamfhand = NamedTemporaryFile(suffix='.bam')
        sort_bam_sam(bam_fpath, sorted_bamfhand.name)

        #sort bam to sam
        bam_fpath = os.path.join(TEST_DATA_DIR, 'samtools', 'seqs.bam')
        sorted_samfhand = NamedTemporaryFile(suffix='.sam')
        sort_bam_sam(bam_fpath, sorted_samfhand.name)

        #sort sam
        sam_fpath = os.path.join(TEST_DATA_DIR, 'samtools', 'seqs.sam')
        sorted_bamfhand = NamedTemporaryFile(suffix='.bam')
        sort_bam_sam(sam_fpath, sorted_bamfhand.name)

    @staticmethod
    def test_standarize_sam():
        'It test that we can add default qualities to the sanger reads'
        sam_fhand = StringIO('''@SQ\tSN:SGN-U576692\tLN:1714
@SQ\tSN:SGN-U572743\tLN:833
@RG\tID:g1\tLB:g1\tSM:g1\tPL:sanger
@RG\tID:g3\tLB:g3\tSM:g3\tPL:sanger
SGN-E200000\t64\tSGN-U572743\t317\t226\t254M24S\t*\t0\t0\tGGATGATKTTAGAG\t*\tAS:i:250\tXS:i:0\tXF:i:0\tXE:i:7\tXN:i:0\tRG:Z:g1
SGN-E40000\t0\tSGN-U576692\t1416\t207\t168M\t*\t0\t0\tAGCCTGATAA\t,,09377777\tAS:i:160\tXS:i:0\tXF:i:3\tXE:i:4\tXN:i:0\tRG:Z:g3
SGN-E40000\t20\tSGN-U576692\t1416\t207\t168M\t*\t0\t0\tAGCCTGATAA\t,,09377777\tAS:i:160\tXS:i:0\tXF:i:3\tXE:i:4\tXN:i:0\tRG:Z:g3
''')
        out_fhand = StringIO()
        standardize_sam(sam_fhand, out_fhand, 20, add_def_qual=True)
        lines = out_fhand.getvalue().splitlines()
        assert 'GGATGATNTTAGAG\t55555555555555\t' in lines[4]
        assert lines[6].startswith('SGN-E40000\t20\t*\t0\t0\t*\t*\t0\t0\t')
        assert lines[4].startswith('SGN-E200000\t0')
    @staticmethod
    def test_realignbam():
        'It test the GATK realigner'
        sam_test_dir = os.path.join(TEST_DATA_DIR, 'samtools')
        bam_path = os.path.join(sam_test_dir, 'seqs.bam')
        reference_path = os.path.join(sam_test_dir, 'reference.fasta')
        out_bam = NamedTemporaryFile(suffix='.bam')
        realign_bam(bam_path, reference_path, out_bam.name, threads=THREADS)
        out_bam.close()

    @staticmethod
    def test_remove_unmapped_reads():
        'Tests remove_unmapped_reads'
        sam = NamedTemporaryFile(suffix='.sam')
        sam.write(SAM)
        sam.flush()
        bam_fhand = NamedTemporaryFile()
        sam2bam(sam.name, bam_fhand.name)

        out_bam_fhand = NamedTemporaryFile()
        out_removed_reads_fhand = NamedTemporaryFile()
        remove_unmapped_reads(bam_fhand, out_bam_fhand, out_removed_reads_fhand)
        reads = open(out_removed_reads_fhand.name).read()
        assert '@SGN-E221406' in reads
        assert 'FFMMMJJ@@755225889>0.' in reads

        out_sam = NamedTemporaryFile(suffix='.sam')
        bam2sam(out_bam_fhand.name, out_sam.name, header=True)
        sam_out = open(out_sam.name).read()
        assert 'SGN-U572743' in sam_out
        assert 'SGN-E221403' in sam_out

    @staticmethod
    def test_get_read_group_info():
        'Tests get_read_group_info'
        sam_sample = '''@SQ\tSN:SGN-U576692\tLN:1714
@SQ\tSN:SGN-U572743\tLN:833
@RG\tID:g1\tLB:g1\tSM:g1\tPL:sanger
@RG\tID:g3\tLB:g3\tSM:g3\tPL:sanger
SGN-E200000\t0\tSGN-U572743\t317\t226\t14M\t*\t0\t0\tGGATGATKTTAGAG\t*\tAS:i:250\tXS:i:0\tXF:i:0\tXE:i:7\tXN:i:0\tRG:Z:g1
SGN-E40000\t0\tSGN-U576692\t1416\t207\t10M\t*\t0\t0\tAGCCTGATAA\t,,09377777\tAS:i:160\tXS:i:0\tXF:i:3\tXE:i:4\tXN:i:0\tRG:Z:g3
SGN-E40000\t20\tSGN-U576692\t1416\t207\t10M\t*\t0\t0\tAGCCTGATAA\t,,09377777\tAS:i:160\tXS:i:0\tXF:i:3\tXE:i:4\tXN:i:0\tRG:Z:g3
'''
        sam_fhand = NamedTemporaryFile(suffix='.sam')
        sam_fhand.write(sam_sample)
        sam_fhand.flush()
        bam_fhand = NamedTemporaryFile(suffix='.bam')
        sam2bam(sam_fhand.name, bam_fhand.name)
        bam_fhand.flush()
        bam = pysam.Samfile(bam_fhand.name, 'rb')
        read_gro_i = get_read_group_info(bam)
        assert read_gro_i == {'g3': {'LB': 'g3', 'SM': 'g3', 'PL': 'sanger'},
                              'g1': {'LB': 'g1', 'SM': 'g1', 'PL': 'sanger'}}


SAM = '''@SQ\tSN:SGN-U576692\tLN:1714
@SQ\tSN:SGN-U572743\tLN:833
@RG\tID:g1\tLB:g1\tSM:g1\tPL:454
@RG\tID:g2\tLB:g2\tSM:g2\tPL:illumina
SGN-U572743\t0\tSGN-U572743\t1416\t207\t168M\t*\t0\t0\tAGCCTGATAAAGGTCTGCCTACGTGTTTTAAGTGGAATCCGTTTCCCCATGTCCAAACCTTCTAAATAGTTTTTTGTGTTAGTTCTTGTATGCCACATACAAAAATTAACAAACTCTTTTGCCACATATGTTCCAGCACGTCAAAGCAACATGTATTTGAGCTACTTT\t558<///035EB@;550300094>>FBF>>88>BBB200>@FFMMMJJ@@755225889>0..14444::FMF@@764444448@;;84444<//,4,.,<<QFBB;::/,,,.69FBB>9:2/.409;@@>88.7,//55;BDK@11,,093777777884241<:7\tAS:i:160\tXS:i:0\tXF:i:3\tXE:i:4\tXN:i:0\tRG:Z:g2\tNM:i:1
SGN-E221402\t0\tSGN-U572743\t1416\t207\t168M\t*\t0\t0\tAGCCTGATAAAGGTCTGCCTACGTGTTTTAAGTGGAATCCGTTTCCCCATGTCCAAACCTTCTAAATAGTTTTTTGTGTTAGTTCTTGTATGCCACATACAAAAATTAACAAACTCTTTTGCCACATATGTTCCAGCACGTCAAAGCAACATGTATTTGAGCTACTTT\t558<///035EB@;550300094>>FBF>>88>BBB200>@FFMMMJJ@@755225889>0..14444::FMF@@764444448@;;84444<//,4,.,<<QFBB;::/,,,.69FBB>9:2/.409;@@>88.7,//55;BDK@11,,093777777884241<:7\tAS:i:160\tXS:i:0\tXF:i:3\tXE:i:4\tXN:i:0\tRG:Z:g1\tNM:i:1
SGN-E221403\t0\tSGN-U572743\t1416\t208\t168M\t*\t0\t0\tAGCCTGATAAAGGTCTGCCTACGTGTTTTAAGTGGAATCCGTTTCCCCATGTCCAAACCTTCTAAATAGTTTTTTGTGTTAGTTCTTGTATGCCACATACAAAAATTAACAAACTCTTTTGCCACATATGTTCCAGCACGTCAAAGCAACATGTATTTGAGCTACTTT\t558<///035EB@;550300094>>FBF>>88>BBB200>@FFMMMJJ@@755225889>0..14444::FMF@@764444448@;;84444<//,4,.,<<QFBB;::/,,,.69FBB>9:2/.409;@@>88.7,//55;BDK@11,,093777777884241<:7\tAS:i:160\tXS:i:0\tXF:i:3\tXE:i:4\tXN:i:0\tRG:Z:g2
SGN-E221404\t0\tSGN-U572743\t1416\t210\t168M\t*\t0\t0\tAGCCTGATAAAGGTCTGCCTACGTGTTTTAAGTGGAATCCGTTTCCCCATGTCCAAACCTTCTAAATAGTTTTTTGTGTTAGTTCTTGTATGCCACATACAAAAATTAACAAACTCTTTTGCCACATATGTTCCAGCACGTCAAAGCAACATGTATTTGAGCTACTTT\t558<///035EB@;550300094>>FBF>>88>BBB200>@FFMMMJJ@@755225889>0..14444::FMF@@764444448@;;84444<//,4,.,<<QFBB;::/,,,.69FBB>9:2/.409;@@>88.7,//55;BDK@11,,093777777884241<:7\tAS:i:160\tXS:i:0\tXF:i:3\tXE:i:4\tXN:i:0\tNM:i:2
SGN-E221404\t256\tSGN-U576692\t1416\t167\t168M\t*\t0\t0\tAGCCTGATAAAGGTCTGCCTACGTGTTTTAAGTGGAATCCGTTTCCCCATGTCCAAACCTTCTAAATAGTTTTTTGTGTTAGTTCTTGTATGCCACATACAAAAATTAACAAACTCTTTTGCCACATATGTTCCAGCACGTCAAAGCAACATGTATTTGAGCTACTTT\t558<///035EB@;550300094>>FBF>>88>BBB200>@FFMMMJJ@@755225889>0..14444::FMF@@764444448@;;84444<//,4,.,<<QFBB;::/,,,.69FBB>9:2/.409;@@>88.7,//55;BDK@11,,093777777884241<:7\tAS:i:160\tXS:i:0\tXF:i:3\tXE:i:4\tXN:i:0\tRG:Z:g2\tNM:i:1
SGN-E221405\t0\tSGN-U576692\t1416\t304\t168M\t*\t0\t0\tAGCCTGATAAAGGTCTGCCTACGTGTTTTAAGTGGAATCCGTTTCCCCATGTCCAAACCTTCTAAATAGTTTTTTGTGTTAGTTCTTGTATGCCACATACAAAAATTAACAAACTCTTTTGCCACATATGTTCCAGCACGTCAAAGCAACATGTATTTGAGCTACTTT\t558<///035EB@;550300094>>FBF>>88>BBB200>@FFMMMJJ@@755225889>0..14444::FMF@@764444448@;;84444<//,4,.,<<QFBB;::/,,,.69FBB>9:2/.409;@@>88.7,//55;BDK@11,,093777777884241<:7\tAS:i:160\tXS:i:0\tXF:i:3\tXE:i:4\tXN:i:0\tRG:Z:g1\tNM:i:0\tX0:i:2
SGN-E221406\t4\tSGN-U576692\t1416\t100\t168M\t*\t0\t0\tAGCCTGATAAAGGTCTGCCTACGTGTTTTAAGTGGAATCCGTTTCCCCATGTCCAAACCTTCTAAATAGTTTTTTGTGTTAGTTCTTGTATGCCACATACAAAAATTAACAAACTCTTTTGCCACATATGTTCCAGCACGTCAAAGCAACATGTATTTGAGCTACTTT\t558<///035EB@;550300094>>FBF>>88>BBB200>@FFMMMJJ@@755225889>0..14444::FMF@@764444448@;;84444<//,4,.,<<QFBB;::/,,,.69FBB>9:2/.409;@@>88.7,//55;BDK@11,,093777777884241<:7\tAS:i:160\tXS:i:0\tXF:i:3\tXE:i:4\tXN:i:0\tRG:Z:g1\tNM:i:3\tX0:i:1
'''

class SamStatsTest(unittest.TestCase):
    'Tests sam stat calculators'

    @staticmethod
    def test_bam_distribs():
        'test bam coverage distrib'
        sam = NamedTemporaryFile(suffix='.sam')
        sam.write(SAM)
        sam.flush()
        bam_fhand = NamedTemporaryFile()
        sam2bam(sam.name, bam_fhand.name)

        summary_fhand = StringIO()

        distribs = bam_distribs(bam_fhand, 'coverage',
                                summary_fhand=summary_fhand)
        expected = [2547]
        assert distribs[('platform', '454')]['distrib'] == expected
        assert 'average: 0.13' in  summary_fhand.getvalue()

        distribs = bam_distribs(bam_fhand, 'mapq')
        assert distribs[('platform', '454')]['distrib'][0] == 1

        distribs = bam_distribs(bam_fhand, 'mapq', sample_size=100)
        assert distribs[('platform', '454')]['distrib'][0] == 1

        distribs = bam_distribs(bam_fhand, 'edit_distance')
        assert distribs[('platform', '454')]['distrib'][0] == 1

    @staticmethod
    def test_sample_bam():
        'it tests sample bam function'
        sam = NamedTemporaryFile(suffix='.sam')
        sam.write(SAM)
        sam.flush()
        bam_fhand = NamedTemporaryFile()
        sam2bam(sam.name, bam_fhand.name)
        bam_fhand.flush()
        out_bam = NamedTemporaryFile(suffix='.bam')
        sample_bam(bam_fhand, out_bam, 2)
        out_sam = NamedTemporaryFile(suffix='.sam')
        bam2sam(out_bam.name, out_sam.name, header=True)

        sam = open(out_sam.name).read().splitlines()
        assert len(sam) == 6

    @staticmethod
    def test_general_mapping_stats():
        'General mapping statistics'
        sam = NamedTemporaryFile(suffix='.sam')
        sam.write(SAM)
        sam.flush()
        bam_fhand = NamedTemporaryFile()
        sam2bam(sam.name, bam_fhand.name)

        out_fhand = StringIO()

        bam_general_stats(bam_fhand, out_fhand)
        result = out_fhand.getvalue()
        assert 'illumina\t3\t100.0' in result
        assert 'Secondary alignments: 1' in result
        assert 'Reads with one X0 best alignment: 1' in result

if	__name__	==	"__main__":
    #import sys;sys.argv = ['', 'SamStatsTest.test_general_mapping_stats']
    unittest.main()
