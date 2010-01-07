'''
It test sam module functions
Created on 05/01/2010

@author: peio
'''
import unittest, os

from biolib.utils.misc_utils import DATA_DIR
from biolib.utils.misc_utils import NamedTemporaryDir

from biolib.sam import (bam2sam, sam2bam, add_tag_to_bam, merge_bam)
class SamTest(unittest.TestCase):
    'It test sam tools related functions'

    @staticmethod
    def testbam2sam():
        'It test bam2sam function'
        tempdir = NamedTemporaryDir()
        bampath = os.path.join(DATA_DIR, 'seq.bam')
        sampath = os.path.join(tempdir.get_name(), 'seq.sam')
        bam2sam(bampath, sampath)
        assert 'SN:SGN-U572743' in open(sampath).readline()
        tempdir.close()

    @staticmethod
    def testsam2bam():
        'It test sam2bam function'
        tempdir = NamedTemporaryDir()
        bampath = os.path.join(DATA_DIR, 'seq.bam')
        sampath = os.path.join(tempdir.get_name(), 'seq.sam')
        bam2sam(bampath, sampath)
        assert 'SN:SGN-U572743' in open(sampath).readline()

        newbam = os.path.join(tempdir.get_name(), 'new.bam')
        sam2bam(sampath, newbam)
        assert open(newbam).read() == open(bampath).read()
        tempdir.close()

    @staticmethod
    def test_add_tag_to_bam():
        'It tests the add tag to bam'
        tags = {'RG':'group1'}
        tempdir = NamedTemporaryDir()
        bampath = os.path.join(DATA_DIR, 'seq.bam')
        newbam = os.path.join(tempdir.get_name(), 'new.bam')
        add_tag_to_bam(bampath, tags, newbam)
        sam = bam2sam(newbam)
        assert  "RG:Z:group1" in sam
        tempdir.close()

    @staticmethod
    def test_merge_bam():
        'It test that we can merge bams'
        bampath = os.path.join(DATA_DIR, 'seq.bam')
        bampath_list = [bampath, bampath]
        tempdir = NamedTemporaryDir()
        newbam = os.path.join(tempdir.get_name(), 'new.bam')
        header = '@RG    ID:group1 LB:group1'
        merge_bam(bampath_list, newbam, header=header)
        merged_sam = bam2sam(newbam)
        assert '@RG    ID:group1 LB:group1' in merged_sam



if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
