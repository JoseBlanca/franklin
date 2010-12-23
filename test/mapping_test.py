'''
Created on 2010 aza 30

@author: peio

It test the mapping module of franklin
'''

import unittest, os, StringIO
from os.path import join, exists
from tempfile import NamedTemporaryFile

from franklin.utils.misc_utils import DATA_DIR, NamedTemporaryDir
from franklin.mapping import (map_reads_with_gmap, gmap_gff_to_sam,
                              map_reads_with_boina, _correct_cigar)
from franklin.sam import bam2sam

class GmapTest(unittest.TestCase):
    'It test the cmap mapper'

    @staticmethod
    def test_gmap_mapper():
        'It test the cmap mapper'
        mappers_dir = join(DATA_DIR, 'mappers')
        cmap_dir = join(DATA_DIR, 'mappers', 'gmap')
        work_dir = NamedTemporaryDir()
        temp_genome = join(work_dir.name, 'genome.fa')
        os.symlink(join(mappers_dir, 'genome.fa'), temp_genome)
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
    def test_gmap_gff3_writer():
        'It test the cmap mapper'
        mappers_dir = join(DATA_DIR, 'mappers')
        cmap_dir = join(DATA_DIR, 'mappers', 'gmap')

        genome_fpath = join(mappers_dir, 'genome.fa')
        in_gmap_gff3 = join(cmap_dir, 'gmap_output.gff3')
        reads_fpath  = join(cmap_dir, 'lb_lib1.pl_sanger.sm_sam2.fa')
        out_sam_fhand = NamedTemporaryFile(suffix='.bam')
        gmap_gff_to_sam(open(in_gmap_gff3), open(genome_fpath),
                        open(reads_fpath), out_sam_fhand)

        result = open(out_sam_fhand.name).read()
        out_sam_fhand.close()
        assert 'seq2\t0\tSL2.30ch00\t241\t255\t36M2I204M' in result
        assert 'seq3\t0\tSL2.30ch00\t241\t255\t48M1D191M' in result
        assert 'seq8_rev\t16\tSL2.30ch00\t3441\t255\t240M' in result
        assert 'seq10_repeated\t4\tSL2.30ch01\t1\t255\t240M' not in result


        out_sam_fhand = NamedTemporaryFile(suffix='.bam')
        gmap_gff_to_sam(open(in_gmap_gff3), open(genome_fpath),
                        open(reads_fpath), out_sam_fhand, keep_unmapped=True)
        result = open(out_sam_fhand.name).read()
        #print result
        out_sam_fhand.close()
        assert 'seq10_repeated\t4\tSL2.30ch01\t1\t255\t240M' in result

    @staticmethod
    def test_correct_cigar():
        'It tests cigar gmap bug corrector'
        cigar = ['2M', '1I', '1D', '3S']
        assert _correct_cigar(cigar) == ['3M', '3S']
        cigar = ['2M', '1I', '1D', '3M']
        assert _correct_cigar(cigar) == ['6M']
        cigar = ['2S', '1I', '1D', '3S']
        assert _correct_cigar(cigar) == ['2S', '1M', '3S']
        cigar = ['2S', '1I', '1M', '3S']
        assert _correct_cigar(cigar) == ['2S', '1I', '1M', '3S']
        cigar = ['90M', '1D', '375M', '1I', '1D', '42M']
        assert _correct_cigar(cigar) == ['90M', '1D', '418M']
        cigar = ['218M', '3D', '102M',  '1I', '4M', '2I', '12M', '6D', '15M',
                 '2I', '16M',  '1I', '27M']
        assert _correct_cigar(cigar) == cigar

        cigar = ['218M', '3D', '96M',  '1D', '24M', '1I', '3D', '15M', '2I',
                 '16M', '1I', '27M']
        assert _correct_cigar(cigar) == cigar

class BoinaTest(unittest.TestCase):
    'It test the gmap mapper'

    @staticmethod
    def test_gmap_mapper():
        'It test the cmap mapper'
        mappers_dir = join(DATA_DIR, 'mappers')
        boina_dir = join(DATA_DIR, 'mappers', 'boina')
        work_dir = NamedTemporaryDir()
        temp_genome = join(work_dir.name, 'genome.fa')
        os.symlink(join(mappers_dir, 'genome.fa'), temp_genome)
        os.chdir(work_dir.name)

        reads_fpath = join(boina_dir, 'lb_lib2.pl_illumina.sm_sam2.sfastq')
        out_bam_fhand = NamedTemporaryFile(suffix='.bam')
        parameters = {'threads':None,
                      'colorspace':False,
                      'reads_length':'short',
                      'java_conf':{}}
        map_reads_with_boina(temp_genome, reads_fpath, out_bam_fhand.name,
                            parameters)
        sam_fhand = NamedTemporaryFile(suffix='.sam')
        bam2sam(out_bam_fhand.name, sam_fhand.name, header=True)
        result = open(sam_fhand.name).read()
        #print result

if __name__ == "__main__":
#    import sys;sys.argv = ['', 'GmapTest.test_correct_cigar']
    unittest.main()


