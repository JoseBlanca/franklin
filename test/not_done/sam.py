'''
Created on 2011 api 19

@author: peio
'''
import unittest
from StringIO import StringIO
from tempfile import NamedTemporaryFile
from franklin.sam import bam2sam
from franklin.not_done.sam import sam_creator

class SamCreatorTest(unittest.TestCase):
    'It tests all function related to sam creation'

    @staticmethod
    def test_sam_creator():
        'It test sam creator'
        alignement = '''ref\taggttttataaaacAAAAaattaagtctacagag-caacta
sample\taggttttataaaacAAA-aattaagtctacagagtcaacta
read\taggttttataaaacAA-Aaattaagtctacagagtcaacta
read\taggttttataaaacA-AAaattaagtctacagagtcaacta
read\taggttttataaaac-AAAaattaagtctacagagtcaacta'''

        fhand = StringIO(alignement)
        out_bam_fhand = NamedTemporaryFile(suffix = '.bam')
        out_ref_fhand  = NamedTemporaryFile(suffix = '.fasta')


        sam_creator(fhand, out_bam_fhand.name, out_ref_fhand.name)

        ref = open(out_ref_fhand.name).read()

        assert  '>ref\naggttttataaaacAAAAaattaagtctacagagcaacta' in  ref
        out_sam_fhand = NamedTemporaryFile(suffix = '.sam')
        bam2sam(out_bam_fhand.name, out_sam_fhand.name)
        bam = open(out_sam_fhand.name).read()

        assert  '16M1D17M1I6M' in bam
        assert  '15M1D18M1I6M' in bam
        assert  '14M1D19M1I6M' in bam

        alignement = '''ref\taggttttataaaac----aattaagtctacagagcaacta
sample\taggttttataaaacAAATaattaagtctacagagcaacta
read\taggttttataaaac****aaAtaa
read1\t ggttttataaaac****aaAtaaTt
read2\t     ttataaaacAAATaattaagtctaca
read3\t        CaaaT****aattaagtctacagagcaac
read4\t          aaT****aattaagtctacagagcaact
read5\t            T****aattaagtctacagagcaacta'''

        fhand = StringIO(alignement)
        out_bam_fhand = NamedTemporaryFile(suffix = '.bam')
        out_ref_fhand  = NamedTemporaryFile(suffix = '.fasta')


        sam_creator(fhand, out_bam_fhand.name, out_ref_fhand.name)
        out_sam_fhand = NamedTemporaryFile(suffix = '.sam')
        bam2sam(out_bam_fhand.name, out_sam_fhand.name)
        bam = open(out_sam_fhand.name).read()
        assert 'AGGTTTTATAAAACAAATA' in bam
        assert '20M' in bam

if __name__    ==    "__main__":
    #import sys;sys.argv = ['', 'SamCreatorTest.test_sam_creator']
    unittest.main()

