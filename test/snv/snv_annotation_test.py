'''
Created on 16/02/2010

@author: jose
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

from __future__ import division
from tempfile import NamedTemporaryFile
import unittest, os
from Bio.Seq import UnknownSeq
from Bio.SeqFeature import FeatureLocation

from franklin.utils.misc_utils import DATA_DIR
from franklin.seq.readers import seqs_in_file
from franklin.seq.seqs import SeqWithQuality, SeqFeature, Seq
from franklin.seq.writers import SequenceWriter
from franklin.snv.snv_annotation import (SNP, INSERTION, DELETION, INVARIANT,
                                         INDEL, COMPLEX, TRANSITION,
                                         TRANSVERSION,
                                         _remove_bad_quality_alleles,
                                         _add_default_sanger_quality,
                                         calculate_snv_kind,
                                         sorted_alleles,
                                         calculate_maf_frequency,
                                         calculate_snv_variability,
                                         calculate_cap_enzymes,
                                         create_snv_annotator,
                                         variable_in_groupping, UNKNOWN)
from franklin.snv.writers import VariantCallFormatWriter

from franklin.pipelines.pipelines import seq_pipeline_runner

class TestSnvAnnotation(unittest.TestCase):
    'It tests the annotation of SeqRecords with snvs'

    @staticmethod
    def test_snv_annotation():
        'It tests the annotation of SeqRecords with snvs'
        bam_fhand = open(os.path.join(DATA_DIR, 'samtools', 'seqs.bam'))
        seq_fhand = open(os.path.join(DATA_DIR, 'samtools', 'reference.fasta'))

        annotator = create_snv_annotator(bam_fhand=bam_fhand, min_quality=30,
                                         min_num_alleles=2)

        expected_snvs = [1, 3]
        for seq, expected in zip(seqs_in_file(seq_fhand), expected_snvs):
            seq = annotator(seq)
            assert expected == len(seq.features)

    def test_snv_remove_edges(self):
        'It test that ww do not annotate snv in the edges'
        edge_remove_settings = {'454':(30, None),
                                'sanger':(None, None)}

        bam_fhand = open(os.path.join(DATA_DIR, 'samtools', 'seqs.bam'))
        seq_fhand = open(os.path.join(DATA_DIR, 'samtools', 'reference.fasta'))

        annotator = create_snv_annotator(bam_fhand=bam_fhand, min_quality=30,
                                         read_edge_conf=edge_remove_settings)
        seqs = list(seqs_in_file(seq_fhand))
        seq = seqs[1]
        seq = annotator(seq)
        assert len(seq.features) == 2

    @staticmethod
    def test_snv_kind():
        'It tests that we can infer the snv kind given a seqrecord with an snv'
        seq = SeqWithQuality(UnknownSeq(100))
        alleles = {('A', INVARIANT): {'read_names':['r1']}}
        feat = SeqFeature(location=FeatureLocation(3, 3), type='snv',
                          qualifiers={'alleles':alleles})
        seq.features.append(feat)

        assert calculate_snv_kind(feat) == INVARIANT

        alleles = {('A', DELETION): {'read_names':['r1']}}
        feat = SeqFeature(location=FeatureLocation(3, 3), type='snv',
                          qualifiers={'alleles':alleles})
        seq.features.append(feat)
        assert calculate_snv_kind(feat) == DELETION

        alleles = {('A', INVARIANT): {'read_names':['r1']},
                   ('A', INSERTION): {'read_names':['r1']},
                   ('A', INSERTION): {'read_names':['r1']},
                   ('A', DELETION): {'read_names':['r1']}}
        feat = SeqFeature(location=FeatureLocation(3, 3), type='snv',
                          qualifiers={'alleles':alleles})
        seq.features.append(feat)
        assert calculate_snv_kind(feat) == INDEL

        alleles = {('A', INVARIANT): {'read_names':['r1']},
                   ('A', SNP): {'read_names':['r1']},
                   ('A', INSERTION): {'read_names':['r1']},
                   ('A', DELETION): {'read_names':['r1']}}
        feat = SeqFeature(location=FeatureLocation(3, 3), type='snv',
                          qualifiers={'alleles':alleles})
        seq.features.append(feat)
        assert calculate_snv_kind(feat) == COMPLEX

        alleles = {('A', INVARIANT): {},
                   ('T', SNP):{}}
        feat = SeqFeature(location=FeatureLocation(3, 3), type='snv',
                          qualifiers={'alleles':alleles})
        assert calculate_snv_kind(feat) == SNP

        #detailed
        assert calculate_snv_kind(feat, detailed=True) == TRANSVERSION

        alleles = {('A', INVARIANT): {},
                   ('C', SNP):{}}
        feat = SeqFeature(location=FeatureLocation(3, 3), type='snv',
                          qualifiers={'alleles':alleles})
        assert calculate_snv_kind(feat, detailed=True) == TRANSVERSION

        alleles = {('A', SNP): {}}
        feat = SeqFeature(location=FeatureLocation(3, 3), type='snv',
                          qualifiers={'alleles':alleles})
        assert calculate_snv_kind(feat, detailed=True) == UNKNOWN
    @staticmethod
    def test_bad_allele_removal():
        'It tests that we can get rid of the alleles with not enough quality'
        min_quality = 45
        default_sanger_quality = 25

        alleles = {('A', SNP): {'qualities':[29, 22],
                                'orientations':[True, True]}}
        _remove_bad_quality_alleles(alleles, min_quality)
        assert len(alleles) == 0

        alleles = {('A', SNP): {'qualities':[29, 22],
                                'orientations':[True, True]}}
        _remove_bad_quality_alleles(alleles, min_quality=32)
        assert len(alleles) == 1

        alleles = {('A', SNP): {'qualities':[23, 22],
                                'orientations':[True, False]}}
        _remove_bad_quality_alleles(alleles, min_quality)
        assert len(alleles) == 0

        alleles = {('A', SNP): {'qualities':[20, 22, None],
                                'orientations':[True, True, False],
                                'read_groups':['rg1', 'rg1', 'rg1']}}
        _add_default_sanger_quality(alleles, default_sanger_quality,
                                    read_groups_info={'rg1':{'PL':'sanger'}})
        _remove_bad_quality_alleles(alleles, min_quality)
        assert len(alleles) == 0

    @staticmethod
    def test_sort_alleles():
        'It tests that we can get the alleles'
        alleles = {('A', INVARIANT): {'read_names':['r1']},
                   ('A', DELETION): {'read_names':['r2', 'r3']}}
        feat = SeqFeature(location=FeatureLocation(3, 3), type='snv',
                          qualifiers={'alleles':alleles})
        alleles = sorted_alleles(feat)
        assert alleles[0]['kind'] == DELETION

    @staticmethod
    def test_maf_alleles():
        'It tests that we can calculate the major allele frequency'
        alleles = {('A', INVARIANT): {'read_groups':{'r1':1}},
                   ('A', DELETION): {'read_groups':{'r2':1, 'r3':1}}}
        feat = SeqFeature(location=FeatureLocation(3, 3), type='snv',
                          qualifiers={'alleles':alleles,
                                      'read_groups':{}})
        maf = calculate_maf_frequency(feat)
        assert maf == 2 / 3

        #with groups
        alleles = {('A', INVARIANT): {'read_groups':{'r1':1},},
                   ('A', DELETION):  {'read_groups':{'r1':1, 'r2':1}}}
        feat = SeqFeature(location=FeatureLocation(3, 3), type='snv',
                          qualifiers={'alleles':alleles,
                                      'read_groups':{'r1':{'LB':'l1'},
                                                     'r2':{'LB':'l2'}}})
        maf = calculate_maf_frequency(feat, groups=['l1'],
                                      group_kind='libraries')
        assert maf == 1 / 2
        maf = calculate_maf_frequency(feat, groups=['l1', 'l2'],
                                      group_kind='libraries')
        assert maf == 2 / 3
        maf = calculate_maf_frequency(feat, groups=['l4'],
                                      group_kind='libraries')
        assert maf is None


    @staticmethod
    def test_snv_variability():
        'It tests that we can calculate the snv variability in a sequence'
        seq = SeqWithQuality(UnknownSeq(100))
        feat1 = SeqFeature(location=FeatureLocation(3, 3), type='snv',
                          qualifiers={'alleles':{}})
        feat2 = SeqFeature(location=FeatureLocation(30, 30), type='snv',
                          qualifiers={'alleles':{}})
        seq.features.extend([feat1, feat2])
        maf = calculate_snv_variability(seq)
        assert maf == 2 / 100


    @staticmethod
    def test_snv_cap():
        'It tests that we can get the enzymes that are a cap'
        reference = SeqWithQuality(seq=Seq('Actgacttactgtca'), name='ref')
        alleles = {('C', SNP): None,
                   ('T', INVARIANT): None}

        feat1 = SeqFeature(location=FeatureLocation(7, 7), type='snv',
                           qualifiers={'alleles':alleles})
        enzymes = calculate_cap_enzymes(feat1, reference, True)
        assert 'HinfI' in enzymes

        #only with the common enzymes
        feat1 = SeqFeature(location=FeatureLocation(7, 7), type='snv',
                           qualifiers={'alleles':alleles})
        enzymes = calculate_cap_enzymes(feat1, reference, False)
        assert 'HinfI' in enzymes

        # With a deletion
        reference = SeqWithQuality(seq=Seq('ATGATGATGgaaattcATGATGATGTGGGAT'),
                                   name='ref')
        alleles = {('A', DELETION): None,
                   ('A', INVARIANT): None}

        feat1 = SeqFeature(location=FeatureLocation(11, 11), type='snv',
                           qualifiers={'alleles':alleles})
        enzymes = calculate_cap_enzymes(feat1, reference, True)
        assert 'EcoRI'  in enzymes

        #with an insertion
        seq = 'ATGATGATG' + 'gaattc' + 'ATGATGATGTGGGAT'
        reference = SeqWithQuality(seq=Seq(seq), name='ref')
        alleles = {('A', INSERTION): None,
                   ('A', INVARIANT): None}

        feat1 = SeqFeature(location=FeatureLocation(11, 11), type='snv',
                           qualifiers={'alleles':alleles})
        enzymes = calculate_cap_enzymes(feat1, reference, True)
        assert 'EcoRI' in enzymes

        #with only one allele
        reference = SeqWithQuality(seq=Seq('Actgacttactgtca'), name='ref')
        alleles = {('C', SNP): None}

        feat1 = SeqFeature(location=FeatureLocation(7, 7), type='snv',
                           qualifiers={'alleles':alleles})
        enzymes = calculate_cap_enzymes(feat1, reference, True)
        assert not enzymes

class TestSnvPipeline(unittest.TestCase):
    'It tests the annotation of SeqRecords with snvs using the pipeline'

    @staticmethod
    def test_snv_annotation_bam():
        'It tests the annotation of SeqRecords with snvs'
        pipeline = 'snv_bam_annotator'

        bam_fhand = open(os.path.join(DATA_DIR, 'samtools', 'seqs.bam'))
        seq_fhand = open(os.path.join(DATA_DIR, 'samtools', 'reference.fasta'))

        configuration = {'snv_bam_annotator': {'bam_fhand':bam_fhand,
                                               'min_quality':30,
                                               'min_num_alleles':2}}

        in_fhands = {}
        in_fhands['in_seq'] = seq_fhand
        seq_fhand = NamedTemporaryFile()
        seq_writer = SequenceWriter(seq_fhand, file_format='repr')
        snv_fhand = NamedTemporaryFile()
        snv_writer = VariantCallFormatWriter(snv_fhand,
                                             reference_name='reference')

        seq_pipeline_runner(pipeline, configuration, in_fhands,
                            writers={'seq':seq_writer, 'snv':snv_writer})

        sequences = list(seqs_in_file(open(seq_fhand.name)))
        num_alleles = 0
        for seq in sequences:
            num_alleles += len(list(seq.get_features('snv')))
        assert num_alleles == 4

        vcf = open(snv_fhand.name).read()
        assert '66' in vcf
        assert '55' in vcf
        assert 'D2' in vcf
        assert 'IAA' in vcf

        #now taking into account the alleles different than the reference
        configuration = {'snv_bam_annotator': {'bam_fhand':bam_fhand,
                                               'min_quality':30,
                                               'min_num_alleles':1}}
        io_fhands = {}
        io_fhands['in_seq'] = seq_fhand
        seq_fhand = NamedTemporaryFile()
        seq_writer = SequenceWriter(seq_fhand, file_format='repr')
        snv_fhand = NamedTemporaryFile()
        snv_writer = VariantCallFormatWriter(snv_fhand,
                                             reference_name='reference')
        seq_pipeline_runner(pipeline, configuration, io_fhands,
                            writers={'seq':seq_writer, 'snv':snv_writer})

        sequences = list(seqs_in_file(open(seq_fhand.name)))
        num_alleles = 0
        for seq in sequences:
            num_alleles += len(list(seq.get_features('snv')))

        assert num_alleles == 8

    @staticmethod
    def test_variable_in_read_group():
        'It test variable_in_reaggroups function'

        alleles = {('A', SNP): {'read_groups':{'rg1':1, 'rg2':1}},
                   ('T', INVARIANT): {'read_groups':{'rg1':1, 'rg3':1}}}

        feature = SeqFeature(location=FeatureLocation(3, 3), type='snv',
                             qualifiers={'alleles':alleles,
                                         'read_groups':{}})
        assert variable_in_groupping('read_groups', feature, ['rg1'])
        assert not variable_in_groupping('read_groups',
                                         feature, ['rg2', 'rg3'])
        assert variable_in_groupping('read_groups', feature, ['rg2', 'rg3'],
                                     in_union=True)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'TestSnvAnnotation.test_snv_remove_edges']
    unittest.main()
