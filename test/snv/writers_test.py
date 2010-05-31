'''
Created on 05/03/2010

@author: peio
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

import unittest
from franklin.seq.seqs import  SeqFeature, SeqWithQuality, Seq
from Bio.SeqFeature import FeatureLocation, ExactPosition
from franklin.snv.snv_annotation import SNP, INVARIANT
from tempfile import NamedTemporaryFile
from franklin.snv.writers import VariantCallFormatWriter, SnvIlluminaWriter

class VariantCallFormatWriterTest(unittest.TestCase):
    'VariantCallFormatWrite tests'

    def test_basic(self):
        seq_str = 'AAA'
        alleles = {('A', SNP):{'read_groups': ['hola_illumina', 'hola2'],
                               'samples': ['individual1', 'individual2'],
                               'read_names': ['seq16', 'seq17'],
                               'orientations': [True, False],
                               'qualities': [57.0, 35.0],
                               'quality': 57.0,
                               'mapping_qualities': [37, 0]},
                   ('T', INVARIANT):{'read_groups': ['sanger'],
                                     'samples': ['individual2'],
                                     'read_names': ['seq16'],
                                     'orientations': [True],
                                     'qualities': [57.0],
                                     'quality': 57.0,
                                     'mapping_qualities': [37]},
                   ('C', SNP):{'read_groups': ['hola_illumina'],
                               'samples': ['individual3'],
                               'read_names': ['seq6'],
                               'orientations': [True, False],
                               'qualities': [57.0, 35.0],
                               'quality': 57.0,
                               'mapping_qualities': [37, 0]},
                   }
        snv1 = SeqFeature(type='snv', location=FeatureLocation(20, 20),
                         qualifiers={'alleles':alleles,
                                     'filters':{'by_kind':{SNP:False}},
                                     'reference_allele':'T'})
        alleles = {('T', INVARIANT):{'read_groups': ['hola_illumina'],
                                     'samples': ['individual4'],
                                     'read_names': ['seq16'],
                                     'orientations': [True],
                                     'qualities': [57.0],
                                     'quality': 57.0,
                                     'mapping_qualities': [37]},
                   }
        filters = {'by_kind':{SNP:True},
                   'high_variable_reg':{(0.8, 0):True},
                   'is_variable':{('libraries', ('lib1',), False):True,
                                  ('libraries', ('lib2',), False):True}
                    }
        snv2 = SeqFeature(type='snv', location=FeatureLocation(50, 50),
                         qualifiers={'alleles':alleles,
                                     'filters':filters,
                                     'reference_allele':'T'})
        seq = SeqWithQuality(seq=Seq(seq_str), qual=[30, 30, 30],
                             name='AT1G55265.1', features=[snv1, snv2])
        fhand = NamedTemporaryFile(mode='a')
        writer = VariantCallFormatWriter(fhand, 'ref1')
        writer.write(seq)
        writer.close()
        vcf = open(fhand.name).read()
        #print vcf

        assert 'GC=2|1|0:2,1,1'
        assert 'VKS' in vcf
        assert '##FILTER=VKS' in vcf
        assert '##FILTER=VLB1' in vcf
        assert '##FILTER=VLB2' in vcf
        assert '1|2:1,1' in vcf
        assert '.:.' in vcf
        assert 'HVR80;VLB1;VLB2;VKS' in vcf
        assert 'AF=0.2,0.5' in vcf
        assert 'GP=0.50' in vcf

        seq_str = 'ATATATATATATATATATATATATAT' * 50
        seq = SeqWithQuality(seq=Seq(seq_str), qual=[30] * len(seq_str),
                             name='AT1G55265.1', features=[snv1, snv2])
        # test snv_illumina writer using the same seqs
        fhand = NamedTemporaryFile(mode='a')
        writer = SnvIlluminaWriter(fhand)
        writer.write(seq)
        illumina_snv = open(fhand.name).read()
        assert 'AT1G55265.1_21,SNP,[A/C/T]TATATATATATA' in illumina_snv

        fhand = NamedTemporaryFile(mode='a')
        writer = SnvIlluminaWriter(fhand, write_short=True)
        writer.write(seq)
        illumina_snv = open(fhand.name).read()
        assert 'ATATATAT[A/C/T]TATATATATA' in illumina_snv

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
