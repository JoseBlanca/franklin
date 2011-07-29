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
from franklin.snv.snv_annotation import SNP, INVARIANT, INSERTION, DELETION
from tempfile import NamedTemporaryFile
from franklin.snv.writers import (VariantCallFormatWriter, SnvIlluminaWriter,
                                  SnvSequenomWriter)

class VariantCallFormatWriterTest(unittest.TestCase):
    'VariantCallFormatWrite tests'

    def test_basic(self):
        seq_str = 'AAA'
        alleles = {('A', SNP):{'read_groups': {'hola_illumina':2, 'hola2':3},
                               'quality': 57.0},
                   ('T', INVARIANT):{'read_groups': {'sanger':5, 'hola2':4},
                                     'quality': 57.0},
                   ('C', SNP):{'read_groups': {'hola_illumina':6},
                               'quality': 57.0},}
        read_groups = {'hola_illumina':{'LB':'lib1'},
                       'hola2':{'LB':'lib2'},
                       'sanger': {'LB': 'lib3'}}
        snv1 = SeqFeature(type='snv', location=FeatureLocation(20, 20),
                          qualifiers={'alleles':alleles,
                                      'filters':{'by_kind':{SNP:False}},
                                      'reference_allele':'T',
                                      'mapping_quality': 47.9,
                                      'quality': 30.1,
                                      'cap_enzymes':['EcoRI'],
                                      'read_groups':read_groups})
        alleles = {('T', INVARIANT):{'read_groups': {'hola_illumina':1},
                                     'quality': 57.0}}
        filters = {'by_kind':{SNP:True},
                   'high_variable_reg':{(0.8, 0):True},
                   'is_variable':{('libraries', ('lib1',), False, False,
                                   None, None):True,
                                  ('libraries', ('lib2',), False, False,
                                   None, None):True},
                   'is_not_variable':{('libraries', ('lib1',), False, False,
                                       None, None):True},
                    }
        snv2 = SeqFeature(type='snv', location=FeatureLocation(50, 50),
                          qualifiers={'alleles':alleles,
                                      'filters':filters,
                                      'reference_allele':'T',
                                      'mapping_quality': 35.6,
                                      'quality': 29.4,
                                      'read_groups':read_groups})
        seq = SeqWithQuality(seq=Seq(seq_str), qual=[30, 30, 30],
                             name='AT1G55265.1', features=[snv1, snv2])
        fhand = NamedTemporaryFile(mode='a')
        writer = VariantCallFormatWriter(fhand, 'ref1')
        writer.write(seq)
        writer.close()
        vcf = open(fhand.name).read()

        assert 'GC=0|1|2:2,2,1' in vcf or 'GC=0|2|1:2,2,1' in vcf
        assert 'VKS' in vcf
        assert '##FILTER=VKS' in vcf
        assert '##FILTER=VLB1' in vcf
        assert '##FILTER=VLB2' in vcf
        assert '1|2:6,2' in vcf
        assert '.:.' in vcf
        assert 'HVR80;VLB1;VLB2;NVLB1;VKS' in vcf
        assert 'AF=0.3,0.2' in vcf
        assert 'GP=0.60' in vcf
        assert 'EZ=EcoRI' in vcf

        #test writer with a non read_groups grouping
        fhand = NamedTemporaryFile(mode='a')
        writer = VariantCallFormatWriter(fhand, 'ref1', grouping='libraries')
        writer.write(seq)
        writer.close()
        vcf = open(fhand.name).read()
        assert 'LIB1' in vcf

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

class SequenomWriterTest(unittest.TestCase):
    'VariantCallFormatWrite tests'

    def test_basic(self):
        seq_str = 'ATGCATGCATGCACTG'

        filters = 'hola'
        alleles1 = {('A', SNP):{}, ('T', INVARIANT):{}, ('C', SNP):{}}
        snv1 = SeqFeature(type='snv', location=FeatureLocation(4, 4),
                        qualifiers={'alleles':alleles1, 'reference_allele':'T',
                        'filters':filters})

        alleles2 = {('T', INVARIANT):{}, ('---', DELETION):{}}
        snv2 = SeqFeature(type='snv', location=FeatureLocation(8, 8),
                        qualifiers={'alleles':alleles2, 'reference_allele':'T'})

        alleles3 = {('AT', INSERTION):{}, ('T', INVARIANT):{}}
        snv3 = SeqFeature(type='snv', location=FeatureLocation(12, 12),
                        qualifiers={'alleles':alleles3, 'reference_allele':'T'})

        seq = SeqWithQuality(seq=Seq(seq_str), qual=[30] * len(seq_str),
                             name='AT1G55265.1', features=[snv1, snv2, snv3])

        # try with SNP
        fhand = NamedTemporaryFile(mode='a')
        writer = SnvSequenomWriter(fhand)
        position = 4
        writer.write(seq, position)
        sequenom_snv = open(fhand.name).read()
        assert 'ATGC[T/A/C]TGCNNNCANNCTG' in sequenom_snv

        # try with Deletion
        fhand = NamedTemporaryFile(mode='a')
        writer = SnvSequenomWriter(fhand)
        position = 8
        writer.write(seq, position)
        sequenom_snv = open(fhand.name).read()
        assert 'ATGCHTGC[ATG/-]CANNCTG' in sequenom_snv

        # try with Insertion
        fhand = NamedTemporaryFile(mode='a')
        writer = SnvSequenomWriter(fhand)
        position = 12
        writer.write(seq, position)
        sequenom_snv = open(fhand.name).read()
        assert 'ATGCHTGCNNNC[-/AT]CTG' in sequenom_snv


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
