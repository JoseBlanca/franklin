'''
It test sam module functions
Created on 05/01/2010

@author: peio
'''
from tempfile import NamedTemporaryFile
from test.utils.misc_utils_test import NamedTemporariDirTest
import unittest, os

from franklin.utils.misc_utils import DATA_DIR
from StringIO import StringIO

from franklin.sam import (bam2sam, sam2bam, merge_sam, bamsam_converter,
                        add_header_and_tags_to_sam, sort_bam_sam)

class SamTest(unittest.TestCase):
    'It test sam tools related functions'

    @staticmethod
    def testbam2sam():
        'It test bam2sam function'
        bampath = os.path.join(DATA_DIR, 'seq.bam')
        sampath = NamedTemporaryFile(suffix='.sam').name
        bam2sam(bampath, sampath)
        assert 'SN:SGN-U572743' in open(sampath).readline()

    @staticmethod
    def testsam2bam():
        'It test sam2bam function'
        bampath = os.path.join(DATA_DIR, 'seq.bam')
        sampath = NamedTemporaryFile(suffix='.sam').name
        bam2sam(bampath, sampath)
        assert 'SN:SGN-U572743' in open(sampath).readline()

        newbam = NamedTemporaryFile(suffix='.bam')
        sam2bam(sampath, newbam.name)
        assert newbam.read() == open(bampath).read()

    @staticmethod
    def test_format_converter():
        'It test BAM SAM converter'
        #bam to sam
        bampath = os.path.join(DATA_DIR, 'seq.bam')
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
SGN-E200000	0	SGN-U572743	317	226	254M24S	*	0	0	GGATGATCTTAGAGCTGCCATTCAAAAGATGTTAGACACTCCTGGGCCATACTTGTTGGATGTGATTGTACCTCATCAGGAGCATGTTCTACCGATGATTCCCAGTGGCGGTGCTTTCAAAAATGTGATTACGGAGGGTGATGGGAGACGTTCCTATTGACTTTGAGAAGCTACATAACTAGTTCAAGGCATTGTATTATCTAAAATAAACTTAATATTTATGTTTACTTAAAAGTTTTTCATTGTGTGAAGGAAAAAAAAAAAAAAAAAAAAAAAAA	999@7<22-2***-,206433>:?9<,,,66:>00066=??EEAAA?B200002<<@@@=DB99777864:..0::@833099???<@488>></...<:B<<.,,8881288@BBDDBD885@@;;9:/9.,,,99B99233885558=?DKKKDDAA??DKBB=440/0<8?DEDFBB??6@152@@FBMFIIDDDDDDKKKOK@@@@DD:N688BBDDDBBBKKDEDDBN977?<9<111:<??==BKMPKKBB==99>QQYYYYYYYYYYYYQQ	AS:i:250	XS:i:0	XF:i:0	XE:i:7	XN:i:0	RG:Z:g1
SGN-E40000	0	SGN-U576692	1416	207	168M	*	0	0	AGCCTGATAAAGGTCTGCCTACGTGTTTTAAGTGGAATCCGTTTCCCCATGTCCAAACCTTCTAAATAGTTTTTTGTGTTAGTTCTTGTATGCCACATACAAAAATTAACAAACTCTTTTGCCACATATGTTCCAGCACGTCAAAGCAACATGTATTTGAGCTACTTT	558<///035EB@;550300094>>FBF>>88>BBB200>@FFMMMJJ@@755225889>0..14444::FMF@@764444448@;;84444<//,4,.,<<QFBB;::/,,,.69FBB>9:2/.409;@@>88.7,//55;BDK@11,,093777777884241<:7	AS:i:160	XS:i:0	XF:i:3	XE:i:4	XN:i:0	RG:Z:g3
''')
        sam2.flush()
        sam3	= NamedTemporaryFile(suffix='.sam')
        merge_sam(infiles=[sam1,	sam2],	outfile=sam3,	reference=reference)
        sam3.seek(0)
        sam3_content = sam3.read()
        assert	'SN:SGN-U572743'	in	sam3_content.split('\n')[0]
        assert	'SGN-E200000'	in	sam3_content
        assert	'SGN-E221664'	in	sam3_content

        #open('/home/peio/tmp/samtools/merged.sam',	'w').write(sam3.getvalue())

    @staticmethod
    def	test_readgroup_to_sam():
        'It	test	that	we	can	add	readgroup	tah	and	header	to	a	bam'
        sam_value	=	'''@SQ	SN:SGN-U576692	LN:1714
@SQ	SN:SGN-U572743	LN:833
SGN-E221403	0	SGN-U576692	1416	207	168M	*	0	0	AGCCTGATAAAGGTCTGCCTACGTGTTTTAAGTGGAATCCGTTTCCCCATGTCCAAACCTTCTAAATAGTTTTTTGTGTTAGTTCTTGTATGCCACATACAAAAATTAACAAACTCTTTTGCCACATATGTTCCAGCACGTCAAAGCAACATGTATTTGAGCTACTTT	558<///035EB@;550300094>>FBF>>88>BBB200>@FFMMMJJ@@755225889>0..14444::FMF@@764444448@;;84444<//,4,.,<<QFBB;::/,,,.69FBB>9:2/.409;@@>88.7,//55;BDK@11,,093777777884241<:7	AS:i:160	XS:i:0	XF:i:3	XE:i:4	XN:i:0
SGN-E221664	0	SGN-U572743	317	226	254M24S	*	0	0	GGATGATCTTAGAGCTGCCATTCAAAAGATGTTAGACACTCCTGGGCCATACTTGTTGGATGTGATTGTACCTCATCAGGAGCATGTTCTACCGATGATTCCCAGTGGCGGTGCTTTCAAAAATGTGATTACGGAGGGTGATGGGAGACGTTCCTATTGACTTTGAGAAGCTACATAACTAGTTCAAGGCATTGTATTATCTAAAATAAACTTAATATTTATGTTTACTTAAAAGTTTTTCATTGTGTGAAGGAAAAAAAAAAAAAAAAAAAAAAAAA	999@7<22-2***-,206433>:?9<,,,66:>00066=??EEAAA?B200002<<@@@=DB99777864:..0::@833099???<@488>></...<:B<<.,,8881288@BBDDBD885@@;;9:/9.,,,99B99233885558=?DKKKDDAA??DKBB=440/0<8?DEDFBB??6@152@@FBMFIIDDDDDDKKKOK@@@@DD:N688BBDDDBBBKKDEDDBN977?<9<111:<??==BKMPKKBB==99>QQYYYYYYYYYYYYQQ	AS:i:250	XS:i:0	XF:i:0	XE:i:7	XN:i:0
'''

        insam	=	NamedTemporaryFile(prefix='lb_group1.',	suffix='.sam')
        insam.write(sam_value)

        outsam	=	NamedTemporaryFile()
        add_header_and_tags_to_sam(insam,	outsam)

        out_content	=	open(outsam.name).read()
        assert	'RG:Z:group1'	in	out_content
        assert	'SM:group1'in	out_content

    @staticmethod
    def	xtest_split_bam_by_readgroup():
        'This	test	is	still	a	work	in	progress'
        'It	test	that	we	can	split	the	bam	by	RG'
        reference	=	StringIO('''>SGN-U572743
atatata
>SGN-U576692
gcgc''')
        sam1	=	StringIO('''@SQ	SN:SGN-U576692	LN:1714
@SQ	SN:SGN-U572743	LN:833
@RG	ID:g1	LB:g1	SM:g1
@RG	ID:g2	LB:g2	SM:g2
SGN-E221403	0	SGN-U576692	1416	207	168M	*	0	0	AGCCTGATAAAGGTCTGCCTACGTGTTTTAAGTGGAATCCGTTTCCCCATGTCCAAACCTTCTAAATAGTTTTTTGTGTTAGTTCTTGTATGCCACATACAAAAATTAACAAACTCTTTTGCCACATATGTTCCAGCACGTCAAAGCAACATGTATTTGAGCTACTTT	558<///035EB@;550300094>>FBF>>88>BBB200>@FFMMMJJ@@755225889>0..14444::FMF@@764444448@;;84444<//,4,.,<<QFBB;::/,,,.69FBB>9:2/.409;@@>88.7,//55;BDK@11,,093777777884241<:7	AS:i:160	XS:i:0	XF:i:3	XE:i:4	XN:i:0	RG:Z:g2
SGN-E221664	0	SGN-U572743	317	226	254M24S	*	0	0	GGATGATCTTAGAGCTGCCATTCAAAAGATGTTAGACACTCCTGGGCCATACTTGTTGGATGTGATTGTACCTCATCAGGAGCATGTTCTACCGATGATTCCCAGTGGCGGTGCTTTCAAAAATGTGATTACGGAGGGTGATGGGAGACGTTCCTATTGACTTTGAGAAGCTACATAACTAGTTCAAGGCATTGTATTATCTAAAATAAACTTAATATTTATGTTTACTTAAAAGTTTTTCATTGTGTGAAGGAAAAAAAAAAAAAAAAAAAAAAAAA	999@7<22-2***-,206433>:?9<,,,66:>00066=??EEAAA?B200002<<@@@=DB99777864:..0::@833099???<@488>></...<:B<<.,,8881288@BBDDBD885@@;;9:/9.,,,99B99233885558=?DKKKDDAA??DKBB=440/0<8?DEDFBB??6@152@@FBMFIIDDDDDDKKKOK@@@@DD:N688BBDDDBBBKKDEDDBN977?<9<111:<??==BKMPKKBB==99>QQYYYYYYYYYYYYQQ	AS:i:250	XS:i:0	XF:i:0	XE:i:7	XN:i:0	RG:Z:g1
''')
        sam2	=	StringIO('''@SQ	SN:SGN-U576692	LN:1714
@SQ	SN:SGN-U572743	LN:833
@RG	ID:g1	LB:g1	SM:g1
@RG	ID:g3	LB:g3	SM:g3
SGN-E200000	0	SGN-U572743	317	226	254M24S	*	0	0	GGATGATCTTAGAGCTGCCATTCAAAAGATGTTAGACACTCCTGGGCCATACTTGTTGGATGTGATTGTACCTCATCAGGAGCATGTTCTACCGATGATTCCCAGTGGCGGTGCTTTCAAAAATGTGATTACGGAGGGTGATGGGAGACGTTCCTATTGACTTTGAGAAGCTACATAACTAGTTCAAGGCATTGTATTATCTAAAATAAACTTAATATTTATGTTTACTTAAAAGTTTTTCATTGTGTGAAGGAAAAAAAAAAAAAAAAAAAAAAAAA	999@7<22-2***-,206433>:?9<,,,66:>00066=??EEAAA?B200002<<@@@=DB99777864:..0::@833099???<@488>></...<:B<<.,,8881288@BBDDBD885@@;;9:/9.,,,99B99233885558=?DKKKDDAA??DKBB=440/0<8?DEDFBB??6@152@@FBMFIIDDDDDDKKKOK@@@@DD:N688BBDDDBBBKKDEDDBN977?<9<111:<??==BKMPKKBB==99>QQYYYYYYYYYYYYQQ	AS:i:250	XS:i:0	XF:i:0	XE:i:7	XN:i:0	RG:Z:g1
SGN-E40000	0	SGN-U576692	1416	207	168M	*	0	0	AGCCTGATAAAGGTCTGCCTACGTGTTTTAAGTGGAATCCGTTTCCCCATGTCCAAACCTTCTAAATAGTTTTTTGTGTTAGTTCTTGTATGCCACATACAAAAATTAACAAACTCTTTTGCCACATATGTTCCAGCACGTCAAAGCAACATGTATTTGAGCTACTTT	558<///035EB@;550300094>>FBF>>88>BBB200>@FFMMMJJ@@755225889>0..14444::FMF@@764444448@;;84444<//,4,.,<<QFBB;::/,,,.69FBB>9:2/.409;@@>88.7,//55;BDK@11,,093777777884241<:7	AS:i:160	XS:i:0	XF:i:3	XE:i:4	XN:i:0	RG:Z:g3
''')
        sam3	=	NamedTemporaryFile()
        merge_sam(infiles=[sam1,	sam2],	outfile=sam3,	reference=reference)
        print	open(sam3.name).read()

    @staticmethod
    def test_sort_bam():
        'It test that we can sort bams using picard'
        #sort sam
        sam_fpath = os.path.join(DATA_DIR, 'samtools', 'seqs.sam')
        sorted_samfhand = NamedTemporaryFile(suffix='.sam')
        sort_bam_sam(sam_fpath, sorted_samfhand.name)

        #sort bam
        bam_fpath = os.path.join(DATA_DIR, 'samtools', 'seqs.bam')
        sorted_bamfhand = NamedTemporaryFile(suffix='.bam')
        sort_bam_sam(bam_fpath, sorted_bamfhand.name)

        #sort bam to sam
        bam_fpath = os.path.join(DATA_DIR, 'samtools', 'seqs.bam')
        sorted_samfhand = NamedTemporaryFile(suffix='.sam')
        sort_bam_sam(bam_fpath, sorted_samfhand.name)

        #sort sam
        sam_fpath = os.path.join(DATA_DIR, 'samtools', 'seqs.sam')
        sorted_bamfhand = NamedTemporaryFile(suffix='.bam')
        sort_bam_sam(sam_fpath, sorted_bamfhand.name)




#    temp_file = NamedTemporaryDir()
#    split_by_read_group(sam3,	temp_file.name)
#    new_sams  = os.listdir(temp_file.name)










if	__name__	==	"__main__":
	#import	sys;sys.argv	=	['',	'Test.testName']
	unittest.main()
