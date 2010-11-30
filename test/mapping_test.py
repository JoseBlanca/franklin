'''
Created on 2010 aza 30

@author: peio

It test the mapping module of franklin
'''

import unittest, os, StringIO
from os.path import join, exists
from tempfile import NamedTemporaryFile

from franklin.utils.misc_utils import DATA_DIR, NamedTemporaryDir
from franklin.mapping import map_reads_with_gmap, gmap_gff_to_sam
from franklin.sam import bam2sam

class CmapTest(unittest.TestCase):
    'It test the cmap mapper'

    @staticmethod
    def test_cmap_mapper():
        'It test the cmap mapper'
        cmap_dir = join(DATA_DIR, 'mappers', 'gmap')
        work_dir = NamedTemporaryDir()
        temp_genome = join(work_dir.name, 'genome.fa')
        os.symlink(join(cmap_dir, 'genome.fa'), temp_genome)
        os.chdir(work_dir.name)

        reads_fpath = join(cmap_dir, 'lb_lib1.pl_sanger.sm_sam1.fa')
        out_bam_fhand = NamedTemporaryFile(suffix='.bam')
        parameters = {'threads':None}
        map_reads_with_gmap(temp_genome, reads_fpath, out_bam_fhand.name,
                            parameters)

        sam_fhand = NamedTemporaryFile(suffix='.sam')
        bam2sam(out_bam_fhand.name, sam_fhand.name, header=True)
        result = open(sam_fhand.name).read()
        assert exists(out_bam_fhand.name)
        assert '36M2I204M' in result
        assert 'SN:SL2.30ch00' in result
        assert 'seq9_rev_MOD' in result

        os.chdir('/tmp')
        work_dir.close()
        out_bam_fhand.close()
        sam_fhand.close()

    @staticmethod
    def test_cmap_gff3_writer():
        'It tests the cmap gff3 parser'
        cmap_dir = join(DATA_DIR, 'mappers', 'gmap')
        reads_fpath = join(cmap_dir, 'lb_lib1.pl_sanger.sm_sam1.fa')
        ref_fpath   = join(cmap_dir, 'genome.fa')
        out_sam_fhand = StringIO.StringIO()
        cmap_gff3_fpath = join(cmap_dir, 'gmap_output.gff3')
        gmap_gff_to_sam(open(cmap_gff3_fpath), open(ref_fpath),
                        open(reads_fpath), out_sam_fhand)

        print out_sam_fhand.getvalue()


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_gff_writer']
    unittest.main()


