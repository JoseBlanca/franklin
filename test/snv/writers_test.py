'''
Created on 05/03/2010

@author: peio
'''
import unittest
from franklin.seq.seqs import  SeqFeature, SeqWithQuality, Seq
from Bio.SeqFeature import FeatureLocation
from franklin.snv.snv_annotation import SNP, INVARIANT
from tempfile import NamedTemporaryFile
from franklin.snv.writers import VariantCallFormatWriter

class VariantCallFormatWriterTest(unittest.TestCase):
    'VariantCallFormatWrite tests'

    def test_basic(self):
        seq_str = 'AAA'
        alleles = {('A', SNP):{'read_groups': ['hola_illumina'],
                               'samples': ['individual1', 'individual2'],
                               'read_names': ['seq16', 'seq17'],
                               'orientations': [True, False],
                               'qualities': [57.0, 35.0],
                               'quality': 57.0,
                               'mapping_qualities': [37, 0]},
                   ('T', INVARIANT):{'read_groups': ['hola_illumina'],
                                     'samples': ['individual2'],
                                     'read_names': ['seq16'],
                                     'orientations': [True],
                                     'qualities': [57.0],
                                     'quality': 57.0,
                                     'mapping_qualities': [37]},
                   ('C', SNP):{'read_groups': ['hola_illumina'],
                               'samples': ['individual3'],
                               'read_names': ['seq16', 'seq17'],
                               'orientations': [True, False],
                               'qualities': [57.0, 35.0],
                               'quality': 57.0,
                               'mapping_qualities': [37, 0]},
                   }
        snv1 = SeqFeature(type='snv', location=FeatureLocation(50, 50),
                         qualifiers={'alleles':alleles,
                                     'filters':{'by_kind':{SNP:True}},
                                     'reference_allele':'T'})
        alleles = {('T', INVARIANT):{'read_groups': ['hola_illumina'],
                                     'samples': ['individual4'],
                                     'read_names': ['seq16'],
                                     'orientations': [True],
                                     'qualities': [57.0],
                                     'quality': 57.0,
                                     'mapping_qualities': [37]},
                   }
        snv2 = SeqFeature(type='snv', location=FeatureLocation(50, 50),
                         qualifiers={'alleles':alleles,
                                     'filters':{'by_kind':{SNP:True}},
                                     'reference_allele':'T'})
        seq = SeqWithQuality(seq=Seq(seq_str), qual=[30, 30, 30],
                             name='AT1G55265.1', features=[snv1, snv2])
        fhand = NamedTemporaryFile(mode='a')
        writer = VariantCallFormatWriter(fhand, 'ref1')
        writer.write(seq)
        writer.close()
        vcf = open(fhand.name).read()
        assert 'vks' in vcf
        assert '##FILTER=vks' in vcf
        assert '0|2:1,1' in vcf
        assert '.:.' in vcf

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
