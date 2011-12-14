'''
Created on 2010 aza 30

@author: peio

It test the mapping module of franklin
'''

import unittest, os, StringIO
from os.path import join, exists
from tempfile import NamedTemporaryFile

from franklin.utils.misc_utils import TEST_DATA_DIR, NamedTemporaryDir
from franklin.mapping import map_reads_with_gmap, map_reads_with_bwa
from franklin.sam import bam2sam

SOLEXA = '@seq1\n'
SOLEXA += 'TCATTGAAAGTTGAAACTGATAGTAGCAGAGTTTTTTCCTCTGTTTGG\n'
SOLEXA += '+\n'
SOLEXA += 'IIIIIIHIIIIIIIIIIIIIIIIIIUJUAUGJUUJUDFAOUDJOFSUD\n'
SOLEXA += '@seq2\n'
SOLEXA += 'ATATGATTGAAGATATTTCTGGGCTTTAAGGGTTCTTGAGGATTTATA\n'
SOLEXA += '+\n'
SOLEXA += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
SOLEXA += '@seq14\n'
SOLEXA += 'ATATGATTGAAGATATTTCTGGGCTTTAAGGGTTCTTGAGGATTTATA\n'
SOLEXA += '+\n'
SOLEXA += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
SOLEXA += '@seq15\n'
SOLEXA += 'ATATGATTGAAGATATTTCTGGGCTTTAAGGGTTCTTGAGGATTTATA\n'
SOLEXA += '+\n'
SOLEXA += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
SOLEXA += '@seq12\n'
SOLEXA += 'ATATGATTGAAGATATTTCTGGACTTTAAGGGTTCTTGAGGATTTATA\n'
SOLEXA += '+\n'
SOLEXA += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
SOLEXA += '@seq13\n'
SOLEXA += 'ATATGATTGAAGATATTTCTGGACTTTAAGGGTTCTTGAGGATTTATA\n'
SOLEXA += '+\n'
SOLEXA += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
SOLEXA += '@seq16\n'
SOLEXA += 'ATATGATTGAAGATATTTCTGGACTTTAAGGGTTCTTGAGGATTTATA\n'
SOLEXA += '+\n'
SOLEXA += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
SOLEXA += '@seq17\n'
SOLEXA += 'atgtcgtacatattggcattgcagtcagcggtatctagtgctaggtaa\n'
SOLEXA += '+\n'
SOLEXA += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'


class GmapTest(unittest.TestCase):
    'It test the gmap mapper'

    @staticmethod
    def test_gmap_mapper():
        'It test the cmap mapper'
        mappers_dir = join(TEST_DATA_DIR, 'mappers')
        gmap_dir = join(TEST_DATA_DIR, 'mappers', 'gmap')
        work_dir = NamedTemporaryDir()
        temp_genome = join(work_dir.name, 'genome.fa')
        os.symlink(join(mappers_dir, 'genome.fa'), temp_genome)

        reads_fpath = join(gmap_dir, 'lb_lib1.pl_sanger.sm_sam1.fa')

        out_bam_fhand = NamedTemporaryFile(suffix='.bam')
        parameters = {'threads':None, 'kmer':13}
        map_reads_with_gmap(temp_genome, reads_fpath, out_bam_fhand.name,
                            parameters)

        sam_fhand = NamedTemporaryFile(suffix='.sam')
        bam2sam(out_bam_fhand.name, sam_fhand.name, header=True)
        result = open(sam_fhand.name).read()
        assert exists(out_bam_fhand.name)
        assert '36M2I204M' in result
        assert 'SN:SL2.30ch00' in result
        assert 'seq9_rev_MOD' in result

        work_dir.close()
        out_bam_fhand.close()
        sam_fhand.close()

        work_dir = NamedTemporaryDir()
        temp_genome = join(work_dir.name, 'genome.fa')
        os.symlink(join(mappers_dir, 'genome.fa'), temp_genome)

        reads_fpath = join(gmap_dir, 'lb_lib1.pl_sanger.sm_sam1.sfastq')
        out_bam_fhand = NamedTemporaryFile(suffix='.bam')
        unmapped_fhand = StringIO.StringIO()
        parameters = {'threads':None, 'kmer':13,
                      'unmapped_fhand':unmapped_fhand}
        map_reads_with_gmap(temp_genome, reads_fpath, out_bam_fhand.name,
                            parameters)

        sam_fhand = NamedTemporaryFile(suffix='.sam')
        bam2sam(out_bam_fhand.name, sam_fhand.name, header=True)
        result = open(sam_fhand.name).read()
        assert exists(out_bam_fhand.name)
        assert '36M2I204M' in result
        assert 'SN:SL2.30ch00' in result
        assert 'seq9_rev_MOD' in result
        assert '?????????????????' in result
        work_dir.close()
        out_bam_fhand.close()
        sam_fhand.close()

    @staticmethod
    def test_gmap_without_mapping_output():
        '''It test that the gmap doesn't map anything'''

        mappers_dir = join(TEST_DATA_DIR, 'mappers')
        cmap_dir = join(TEST_DATA_DIR, 'mappers', 'gmap')
        work_dir = NamedTemporaryDir()
        temp_genome = join(work_dir.name, 'genome.fa')
        os.symlink(join(mappers_dir, 'genome.fa'), temp_genome)

        reads_fhand = NamedTemporaryFile()
        reads_fhand.write('>seq\natgtgatagat\n')
        reads_fhand.flush()


        out_bam_fhand = NamedTemporaryFile()
        out_bam_fpath = out_bam_fhand.name
        out_bam_fhand.close()
        parameters = {'threads':None, 'kmer':13}
        map_reads_with_gmap(temp_genome, reads_fhand.name, out_bam_fpath,
                            parameters)
        reads_fhand.close()
        temp_sam_fhand = NamedTemporaryFile(suffix='.sam')
        bam2sam(out_bam_fpath, temp_sam_fhand.name, True)
        result = open(temp_sam_fhand.name).read()
        assert 'seq\t4\t*\t0\t0' in result

class BwaTest(unittest.TestCase):
    'It test the bwa mapper'
    @staticmethod
    def test_bwa_mapping():
        '''It test that the gmap doesn't map anything'''
        reference = join(TEST_DATA_DIR, 'blast/arabidopsis_genes')
        work_dir = NamedTemporaryDir()
        reference_fpath = join(work_dir.name, 'arabidopsis_genes')
        os.symlink(reference, reference_fpath)

        reads_fhand = NamedTemporaryFile(suffix='.sfastq')
        reads_fhand.write(SOLEXA)
        reads_fhand.flush()

        out_bam_fhand = NamedTemporaryFile()
        out_bam_fpath = out_bam_fhand.name
        out_bam_fhand.close()

        parameters = {'colorspace': False, 'reads_length':'short',
                      'threads':None, 'java_conf':None}
        map_reads_with_bwa(reference_fpath, reads_fhand.name, out_bam_fpath,
                           parameters)
        test_sam_fhand = NamedTemporaryFile(suffix='sam')
        bam2sam(out_bam_fpath, test_sam_fhand.name)
        result = open(test_sam_fhand.name).read()
        assert 'seq17' in result

        unmapped_fhand = StringIO.StringIO()
        parameters = {'colorspace': False, 'reads_length':'short',
                      'threads':None, 'java_conf':None,
                      'unmapped_fhand':unmapped_fhand}
        map_reads_with_bwa(reference_fpath, reads_fhand.name, out_bam_fpath,
                           parameters)
        assert 'seq17' in unmapped_fhand.getvalue()
        test_sam_fhand = NamedTemporaryFile(suffix='sam')
        bam2sam(out_bam_fpath, test_sam_fhand.name)
        result = open(test_sam_fhand.name).read()
        assert 'seq17' not in result


if __name__ == "__main__":
    import sys;sys.argv = ['', 'BwaTest.test_bwa_mapping']
    unittest.main()
