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
from StringIO import StringIO
import pysam

from franklin.utils.misc_utils import TEST_DATA_DIR, DATA_DIR
from franklin.seq.readers import seqs_in_file
from franklin.seq.seqs import SeqWithQuality, SeqFeature, Seq
from franklin.seq.writers import SequenceWriter, write_seqs_in_file
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
                                         UNKNOWN,
                                         variable_in_groupping,
                                         invariant_in_groupping,
                                         _remove_alleles_by_read_number,
                                         _get_segments_from_cigar,
                                         _locate_segment, IN_FIRST_AND_LAST,
                                         IN_FIRST_POS, IN_LAST_POS,
                                         _get_alleles_from_read,
                                         annotate_pic,
                                         annotate_heterozygosity,
                                         snvs_in_window,
                                         _get_alignment_section,
                                         _make_multiple_alignment)

from franklin.sam import create_bam_index, sam2bam, add_header_and_tags_to_sam,\
    bam2sam
from franklin.snv.writers import VariantCallFormatWriter

from franklin.pipelines.pipelines import seq_pipeline_runner
from franklin.mapping import map_reads_with_bwa

class TestSnvAnnotation(unittest.TestCase):
    'It tests the annotation of SeqRecords with snvs'

    @staticmethod
    def test_snv_annotation_basic():
        'It tests the annotation of SeqRecords with snvs'
        bam_fhand = open(os.path.join(TEST_DATA_DIR, 'samtools', 'seqs.bam'))
        seq_fhand = open(os.path.join(TEST_DATA_DIR, 'samtools',
                                      'reference.fasta'))

        annotator = create_snv_annotator(bam_fhand=bam_fhand, min_quality=30,
                                         min_num_alleles=2)

        expected_snvs = [1, 3]
        for seq, expected in zip(seqs_in_file(seq_fhand), expected_snvs):
            seq = annotator(seq)
            assert expected == len(seq.features)
    @staticmethod
    def test_snv_annotation_massive():
        'It tests the annotation of SeqRecords with snvs'
        #another example

        bam_fhand = NamedTemporaryFile(suffix='.bam')
        sam_fhand = NamedTemporaryFile(prefix='pl_illumina.sm_test.lb_test.',
                                       suffix='.sam')
        seq_fhand = open(os.path.join(TEST_DATA_DIR, 'snv_annotator',
                                      'seqs.sfastq'))
        ref_fhand = open(os.path.join(TEST_DATA_DIR, 'snv_annotator',
                                      'reference.fasta'))

        # prepare the sam file
        parameters = {'colorspace':False, 'reads_length':'short',
                      'threads':None, 'java_conf':None }
        map_reads_with_bwa(ref_fhand.name, seq_fhand.name, bam_fhand.name,
                           parameters)
        bam2sam(bam_fhand.name, sam_fhand.name, header=True)
        new_sam_fhand = NamedTemporaryFile(suffix='.sam')
        add_header_and_tags_to_sam(sam_fhand, new_sam_fhand)
        sam2bam(new_sam_fhand.name, bam_fhand.name)
        #################################################################

        annotator = create_snv_annotator(bam_fhand=bam_fhand, min_quality=20,
                                         min_num_alleles=1,
                                         min_num_reads_for_allele=1)

        for seq in seqs_in_file(ref_fhand):
            seq = annotator(seq)


    @staticmethod
    def test_snv_annotation_with_pic_and_heterozygosity():
        'It tests the pic and heterozygosity annotation of SeqRecords with snvs'
        bam_fhand = open(os.path.join(TEST_DATA_DIR, 'samtools', 'seqs.bam'))
        seq_fhand = open(os.path.join(TEST_DATA_DIR, 'samtools',
                                      'reference.fasta'))

        annotator = create_snv_annotator(bam_fhand=bam_fhand, min_quality=30,
                                         min_num_alleles=2)
        seqs = seqs_in_file(seq_fhand)
        seq = seqs.next()
        seq = annotator(seq)
        assert  round(seq.features[0].qualifiers['pic'], 2) == 0.44
        assert  round(seq.features[0].qualifiers['heterozygosity'], 2) == 0.47

    @staticmethod
    def test_snv_annotation_without_rg():
        'It tests the annotation of SeqRecords with snvs'
        bam_fhand = open(os.path.join(TEST_DATA_DIR, 'samtools',
                                      'seqs_without_rg.bam'))
        seq_fhand = open(os.path.join(TEST_DATA_DIR, 'samtools', 'reference.fasta'))

        # this test should fail because we did not provide the read group
        # parameters
        annotator = create_snv_annotator(bam_fhand=bam_fhand,
                                             min_quality=30,
                                             min_num_alleles=2)

        expected_snvs = [1, 3]
        for seq, expected in zip(seqs_in_file(seq_fhand), expected_snvs):
            try:
                seq = annotator(seq)
                raise ValueError('test_snv_annotation_without_rg Test failed')
            except ValueError:
                pass
        # now with the RG option
        default_bam_platform = 'sanger'
        annotator = create_snv_annotator(bam_fhand=bam_fhand,
                                         min_quality=30,
                                         min_num_alleles=2,
                                        default_bam_platform=default_bam_platform)

        expected_snvs = [1, 3]
        for seq, expected in zip(seqs_in_file(seq_fhand), expected_snvs):
            seq = annotator(seq)
            assert expected == len(seq.features)

    def test_snv_remove_edges(self):
        'It test that we do not annotate snv in the edges'

        REF = 'AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT'

        #Coor     01234567890123  4567890123456789012345678901234
        #ref      AGCATGTTAGATAA**GATAGCTGTGCTAGTAGGCAGTCAGCGCCAT
        #                 012345
        #+r003       gcctaAGATAA
        #                          012345              67890
        #+r004                     ATAGCT..............TCAGC
        #                                       43210
        #-r005                            ttagctTAGGC
        #                                               876543210
        #-r001/2                                        CAGCGCCAT

        SAM = '''@HD\tVN:1.3\tSO:coordinate
@SQ\tSN:ref\tLN:45
r003\t0\tref\t9\t30\t5H6M\t*\t0\t0\tgagccc\t*\tNM:i:1
r004\t0\tref\t16\t30\t6M14N5M\t*\t0\t0\tATAGCaaCAGC\t*
r005\t16\tref\t29\t30\t6H5M\t*\t0\t0\tTAGGa\t*\tNM:i:0
r001/2\t83\tref\t37\t30\t9M\t=\t7\t-39\tCAcCGCCAT\t*
'''

        sam = NamedTemporaryFile(suffix='.sam')
        sam.write(SAM)
        sam.flush()
        bam_fhand = NamedTemporaryFile()
        sam2bam(sam.name, bam_fhand.name)
        create_bam_index(bam_fhand.name)

        reference = SeqWithQuality(seq=Seq(REF), name='ref')
        edge_remove_settings = {'sanger':(None, None)}

        annotator = create_snv_annotator(bam_fhand=bam_fhand, min_quality=30,
                                         read_edge_conf=edge_remove_settings,
                                         default_bam_platform='sanger',
                                         default_sanger_quality=60,
                                         min_num_reads_for_allele=1)
        seq = annotator(reference)
        assert len(seq.features) == 10

        reference = SeqWithQuality(seq=Seq(REF), name='ref')
        edge_remove_settings = {'sanger':(1, 1)}

        annotator = create_snv_annotator(bam_fhand=bam_fhand, min_quality=30,
                                         read_edge_conf=edge_remove_settings,
                                         default_bam_platform='sanger',
                                         default_sanger_quality=60,
                                         min_num_reads_for_allele=1)
        seq = annotator(reference)
        assert len(seq.features) == 7


        reference = SeqWithQuality(seq=Seq(REF), name='ref')
        edge_remove_settings = {'sanger':(2, 1)}

        annotator = create_snv_annotator(bam_fhand=bam_fhand, min_quality=30,
                                         read_edge_conf=edge_remove_settings,
                                         default_bam_platform='sanger',
                                         default_sanger_quality=60,
                                         min_num_reads_for_allele=1)
        seq = annotator(reference)
        assert len(seq.features) == 6

        reference = SeqWithQuality(seq=Seq(REF), name='ref')
        edge_remove_settings = {'sanger':(10, 10)}

        annotator = create_snv_annotator(bam_fhand=bam_fhand, min_quality=30,
                                         read_edge_conf=edge_remove_settings,
                                         default_bam_platform='sanger',
                                         default_sanger_quality=60,
                                         min_num_reads_for_allele=1)
        seq = annotator(reference)
        assert len(seq.features) == 0

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
    def test_remove_alleles_by_read_number():
        'it removes the alleles by read number'
        alleles = {('A', SNP): {'qualities':[29, 22],
                                'read_groups':[True, True]}}
        _remove_alleles_by_read_number(alleles, 2)
        assert len(alleles) == 1

        _remove_alleles_by_read_number(alleles, 1)
        assert len(alleles) == 1

        _remove_alleles_by_read_number(alleles, 3)
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

    @staticmethod
    def test_snvs_in_window():
        'It tests the snvs_in_window with maf and type'
        alleles1 = {('A', DELETION): {'read_groups':{'r1':3}},
                    ('T', INVARIANT): {'read_groups':{'r1':8}}}

        snv1 = SeqFeature(location=FeatureLocation(3, 3), type='snv',
                          qualifiers={'alleles':alleles1,
                                      'read_groups':{'r1':{'LB':'l1'},
                                                     'r2':{'LB':'l2'}}})
        alleles2 = {('A', SNP): {'read_groups':{'r1':2}},
                    ('T', INVARIANT): {'read_groups':{'r1':2}}}

        snv2 = SeqFeature(location=FeatureLocation(7, 7), type='snv',
                          qualifiers={'alleles':alleles2,
                                      'read_groups':{'r1':{'LB':'l1'},
                                                     'r2':{'LB':'l2'}}})
        alleles3 = {('A', SNP): {'read_groups':{'r1':8}},
                    ('T', INVARIANT): {'read_groups':{'r1':8}}}
        snv3 = SeqFeature(location=FeatureLocation(20, 20), type='snv',
                          qualifiers={'alleles':alleles3,
                                      'read_groups':{'r1':{'LB':'l1'},
                                                     'r2':{'LB':'l2'}}})
        alleles4 = {('A', SNP): {'read_groups':{'r1':1}},
                    ('T', INVARIANT): {'read_groups':{'r1':10}}}
        snv4 = SeqFeature(location=FeatureLocation(25, 25), type='snv',
                          qualifiers={'alleles':alleles4,
                                      'read_groups':{'r1':{'LB':'l1'},
                                                     'r2':{'LB':'l2'}}})
        alleles5 = {('A', SNP): {'read_groups':{'r1':1}},
                    ('T', INVARIANT): {'read_groups':{'r1':100}}}
        snv5 = SeqFeature(location=FeatureLocation(31, 31), type='snv',
                          qualifiers={'alleles':alleles5,
                                      'read_groups':{'r1':{'LB':'l1'},
                                                     'r2':{'LB':'l2'}}})

        snvs = [snv1, snv2, snv3, snv4, snv5]

        assert snvs_in_window(snv2, snvs, 8) == 1
        assert snvs_in_window(snv2, snvs, 30, snv_type=SNP) == 1
        assert snvs_in_window(snv2, snvs, 30, snv_type=SNP, maf=0.7) == 1
        assert snvs_in_window(snv3, snvs, 30, maf=0.7) == 1

class TestSnvPipeline(unittest.TestCase):
    'It tests the annotation of SeqRecords with snvs using the pipeline'

    @staticmethod
    def test_snv_annotation_bam():
        'It tests the annotation of SeqRecords with snvs'
        pipeline = 'snv_bam_annotator'

        bam_fhand = open(os.path.join(TEST_DATA_DIR, 'samtools', 'seqs.bam'))
        seq_fhand = open(os.path.join(TEST_DATA_DIR, 'samtools', 'reference.fasta'))

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

        assert num_alleles == 9

    @staticmethod
    def test_variable_in_read_group():
        'It test variable_in_reaggroups function'

        alleles = {('A', SNP): {'read_groups':{'rg1':1, 'rg2':2}},
                   ('T', INVARIANT): {'read_groups':{'rg1':2, 'rg3':2}}}

        snv = SeqFeature(location=FeatureLocation(3, 3), type='snv',
                             qualifiers={'alleles':alleles, 'read_groups':{}})


        assert not variable_in_groupping(snv, 'read_groups', ['rg5'],
                                         reference_free=True)
        assert not variable_in_groupping(snv, 'read_groups', ['rg5'],
                                         reference_free=False)

        assert variable_in_groupping(snv, 'read_groups', ['rg1'],
                                     reference_free=False, maf=0.7)

        assert not variable_in_groupping(snv, 'read_groups', ['rg1'],
                                     reference_free=False, maf=0.6)

        assert variable_in_groupping(snv, 'read_groups', ['rg1'],
                                     reference_free=True, maf=0.7)

        assert not variable_in_groupping(snv, 'read_groups', ['rg1'],
                                     reference_free=True, maf=0.6)

        assert not variable_in_groupping(snv, 'read_groups', ['rg1'],
                                     min_reads_per_allele=2)


        assert not variable_in_groupping(snv, 'read_groups', ['rg2'],
                                     reference_free=True)

        assert variable_in_groupping(snv, 'read_groups', ['rg2'],
                                     reference_free=False)

        assert variable_in_groupping(snv, 'read_groups', ['rg2'],
                                     reference_free=False,
                                     min_reads_per_allele=2)

        assert not variable_in_groupping(snv, 'read_groups', ['rg2'],
                                     reference_free=False,
                                     min_reads_per_allele=3)

        assert variable_in_groupping(snv, 'read_groups', ['rg2', 'rg3'])

        assert not variable_in_groupping(snv, 'read_groups', ['rg2', 'rg3'],
                                         in_union=False)

        assert variable_in_groupping(snv, 'read_groups', ['rg2', 'rg3'],
                                     reference_free=False)

        assert not variable_in_groupping(snv, 'read_groups', ['rg2', 'rg3'],
                                         reference_free=False, in_union=False)

    @staticmethod
    def test_not_variable_in_read_group():
        'It test variable_in_reaggroups function'

        alleles = {('A', SNP): {'read_groups':{'rg1':1, 'rg2':1}},
                   ('T', INVARIANT): {'read_groups':{'rg1':2, 'rg3':2}}}

        snv = SeqFeature(location=FeatureLocation(3, 3), type='snv',
                             qualifiers={'alleles':alleles, 'read_groups':{}})


        assert not invariant_in_groupping(snv, 'read_groups', ['rg5'],
                                          reference_free=True)
        assert invariant_in_groupping(snv, 'read_groups', ['rg5'],
                                      reference_free=False)

        assert not invariant_in_groupping(snv, 'read_groups', ['rg1'],
                                     reference_free=False)

        assert not invariant_in_groupping(snv, 'read_groups', ['rg1'],
                                     reference_free=True)

        assert invariant_in_groupping(snv, 'read_groups', ['rg2'],
                                      reference_free=True)

        assert not invariant_in_groupping(snv, 'read_groups', ['rg2'],
                                     reference_free=False)

        assert invariant_in_groupping(snv, 'read_groups', ['rg3'],
                                     reference_free=False)

        assert invariant_in_groupping(snv, 'read_groups', ['rg3'],
                                      reference_free=True)

        assert not invariant_in_groupping(snv, 'read_groups', ['rg2', 'rg3'],
                                          reference_free=False)
        assert not invariant_in_groupping(snv, 'read_groups', ['rg2', 'rg3'],
                                          reference_free=True)

        assert invariant_in_groupping(snv, 'read_groups', ['rg2', 'rg3'],
                                          reference_free=True, in_union=False)

        assert not invariant_in_groupping(snv, 'read_groups', ['rg2', 'rg3'],
                                          reference_free=False, in_union=False)


        assert invariant_in_groupping(snv, 'read_groups', ['rg5', 'rg6'],
                                          reference_free=False)
        assert not invariant_in_groupping(snv, 'read_groups', ['rg5', 'rg6'],
                                          reference_free=True)

REF = 'AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT'

SAM = '''@HD\tVN:1.3\tSO:coordinate
@SQ\tSN:ref\tLN:45
r001/1\t163\tref\t7\t30\t8M2I4M1D3M\t=\t37\t39\tTAAGATAAAGGATACTG\t*
r002\t0\tref\t9\t30\t3S6M1P1I4M\t*\t0\t0\tAAAAGATAAGGATA\t*
r003\t0\tref\t9\t30\t5H6M\t*\t0\t0\tAGCTAA\t*\tNM:i:1
r004\t0\tref\t16\t30\t6M14N5M\t*\t0\t0\tATAGCTTCAGC\t*
r005\t16\tref\t29\t30\t6H5M\t*\t0\t0\tTAGGC\t*\tNM:i:0
r001/2\t83\tref\t37\t30\t9M\t=\t7\t-39\tCAGCGCCAT\t*
'''


class TestReadPos(unittest.TestCase):
    'It tests the reference to read pos coordinate conversion'

    @staticmethod
    def test_get_segments_from_cigar():
        '''It tests the obtention of reference and reads segments.
        It tests the obtention of the reference limits.
        It tests the obtention of segment types and segment lengths'''

        read = 'TTAGATAAAGGATACTG'
        ref_pos = 6
        cigar = [(0, 8), (1, 2), (0, 4), (2, 1), (0, 3)]
        read_len = len(read)

        (ref_segments, read_segments, ref_limits, read_limits, segment_type,
         segment_lens) = _get_segments_from_cigar(ref_pos, cigar, read_len)

        assert ref_segments == [6, None, 14, 18, 19]
        assert read_segments == [0, 8, 10, None, 14]
        assert ref_limits == [6, 21]
        assert read_limits == [0, 16]
        assert segment_type == [0, 1, 0, 2, 0]
        assert segment_lens == [8, 2, 4, 1, 3]


        read = 'aaaAGATAAGGATA'
        ref_pos = 8
        cigar = '3S6M1P1I4M'
        cigar = [(4, 3), (0, 6),(6, 1), (1, 1), (0, 4)]
        read_len = len(read)

        (ref_segments, read_segments, ref_limits, read_limits, segment_type,
         segment_lens) = _get_segments_from_cigar(ref_pos, cigar, read_len)

        assert ref_segments == [8, None, 14]
        assert read_segments == [3, 9, 10]
        assert ref_limits == [8, 17]
        assert read_limits == [3, 13]
        assert segment_type == [0, 1, 0]
        assert segment_lens == [6, 1, 4]


        read = 'TAGGCTTTACgta'
        ref_pos = 28
        cigar = '6H5M3I2M3S'
        cigar = [(0, 5),(1, 3), (0, 2), (4, 3)]
        read_len = len(read)

        (ref_segments, read_segments, ref_limits, read_limits, segment_type,
         segment_lens) = _get_segments_from_cigar(ref_pos, cigar, read_len)

        assert ref_segments == [28, None, 33]
        assert read_segments == [0, 5, 8]
        assert ref_limits == [28, 34]
        assert read_limits == [0, 9]
        assert segment_type == [0, 1, 0]
        assert segment_lens == [5, 3, 2]


        read = 'ATCGTAG'
        ref_pos = 14
        cigar = '3M4N2M3D2M'
        cigar = [(0, 3),(3, 4), (0, 2), (2, 3), (0, 2)]
        read_len = len(read)

        (ref_segments, read_segments, ref_limits, read_limits, segment_type,
         segment_lens) = _get_segments_from_cigar(ref_pos, cigar, read_len)

        assert ref_segments == [14, 17, 21, 23, 26]
        assert read_segments == [0, None, 3, None, 5]
        assert ref_limits == [14, 27]
        assert read_limits ==  [0, 6]
        assert segment_type == [0, 3, 0, 2, 0]
        assert segment_lens == [3, 4, 2, 3, 2]


        read = 'ATCGTAG'
        ref_pos = 8945400
        cigar = '41M4907N1I33M1D93M4861N34M'
        cigar = [(0, 3),(3, 4), (0, 2), (2, 3), (0, 2)]
        read_len = len(read)
        (ref_segments, read_segments, ref_limits, read_limits, segment_type,
         segment_lens) = _get_segments_from_cigar(ref_pos, cigar, read_len)

    @staticmethod
    def test_locate_segment():
        'It tests the correct location of a read position in the segments'

        ref_segments = [6, 8, 10, None, 13, 15, 16]
        ref_limits = [6, 16]
        segment_lens = [2, 2, 3, 2, 2, 1, 1]

        ref_pos = 8
        segment_index, segment_pos = _locate_segment(ref_pos, ref_segments,
                                                    segment_lens, ref_limits)
        assert segment_index == 1
        assert segment_pos == IN_FIRST_POS

        ref_pos = 12
        segment_index, segment_pos = _locate_segment(ref_pos, ref_segments,
                                                    segment_lens, ref_limits)
        assert segment_index == 2
        assert segment_pos == IN_LAST_POS

        ref_pos = 16
        segment_index, segment_pos = _locate_segment(ref_pos, ref_segments,
                                                    segment_lens, ref_limits)
        assert segment_index == 6
        assert segment_pos == IN_FIRST_AND_LAST

    @staticmethod
    def test_get_alleles_from_read():
        'It tests the correct extraction of information in each read position'
        sam = NamedTemporaryFile(suffix='.sam')
        sam.write(SAM)
        sam.flush()
        bam_fhand = NamedTemporaryFile()
        sam2bam(sam.name, bam_fhand.name)
        create_bam_index(bam_fhand.name)

        bam = pysam.Samfile(bam_fhand.name, "rb")

        reference = SeqWithQuality(seq=Seq(REF), name='ref')
        reference_id = reference.name
        reference_seq = reference.seq

        #Coor     01234567890123  4567890123456789012345678901234
        #ref      AGCATGTTAGATAA**GATAGCTGTGCTAGTAGGCAGTCAGCGCCAT
        #               01234567890123 456
        #+r001/1        TAAGATAAAGGATA*CTG
        #              012345678 90123
        #+r002         aaaAGATAA*GGATA
        #                 012345
        #+r003       gcctaAGCTAA
        #                          012345              67890
        #+r004                     ATAGCT..............TCAGC
        #                                       43210
        #-r005                            ttagctTAGGC
        #                                               876543210
        #-r001/2                                        CAGCGCCAT

        expected = {('r001/1', 5): [],
                    ('r001/1', 6): [('T', INVARIANT, None, False)],
                    ('r001/1', 7): [('A', SNP, None, False)],
                    ('r001/1', 13): [('A', INVARIANT, None, False),
                                     ('AG', INSERTION, None, False)],
                    ('r001/1', 14): [('G', INVARIANT, None, False)],
                    ('r001/1', 17): [('A', INVARIANT, None, False)],
                    ('r001/1', 18): [('-', DELETION, None, False)],
                    ('r001/1', 22): [],
                    ('r002', 5): [],
                    ('r002', 7): [],
                    ('r002', 8): [('A', INVARIANT, None, False)],
                    ('r002', 14): [('G', INVARIANT, None, False)],
                    ('r002', 18): [],
                    ('r003', 7): [],
                    ('r003', 8): [('A', INVARIANT, None, False)],
                    ('r003', 14): [],
                    ('r004', 20): [('T', INVARIANT, None, False)],
                    ('r004', 21): [],
                    ('r004', 34): [],
                    ('r004', 35): [('T', INVARIANT, None, False)],
                    ('r005', 27): [],
                    ('r005', 28): [('T', INVARIANT, None, True)],
                    ('r005', 29): [('A', INVARIANT, None, True)],
                    ('r005', 32): [('C', INVARIANT, None, True)],
                    }

        for column in bam.pileup(reference=reference_id):
            ref_pos = column.pos
            ref_allele = reference_seq[ref_pos].upper()
            for pileup_read in column.pileups:
                read_name = pileup_read.alignment.qname
                if (read_name, ref_pos) == ('r005', 28):
                    pass
                alleles, read_limits = _get_alleles_from_read(ref_allele,
                                                              ref_pos,
                                                              pileup_read)
                if (read_name, ref_pos) in expected:
                    if alleles != expected[(read_name, ref_pos)]:
                        print repr(alleles)
                        print expected[(read_name, ref_pos)]
                    assert alleles == expected[(read_name, ref_pos)]


        # It should work without errors
        sam_fpath = os.path.join(TEST_DATA_DIR, 'samtools', 'fail.sam')
        bam_fhand = NamedTemporaryFile()
        sam2bam(sam_fpath, bam_fhand.name)

        create_bam_index(bam_fhand.name)
        bam = pysam.Samfile(bam_fhand.name, "rb")

        seq_fpath     = os.path.join(TEST_DATA_DIR, 'samtools', 'ref2.fasta')
        seqs          = seqs_in_file(open(seq_fpath))
        reference     = seqs.next()
        reference_id  = reference.name
        reference_seq = reference.seq
        for column in bam.pileup(reference=reference_id):
            ref_pos = column.pos
            ref_allele = reference_seq[ref_pos].upper()
            for pileup_read in column.pileups:
                read_name = pileup_read.alignment.qname
                alleles, read_limits = _get_alleles_from_read(ref_allele,
                                                              ref_pos,
                                                              pileup_read)
    @staticmethod
    def test_get_aligned_read_section():
        'It checks that we can get aligned sections from a read'
        #Coor               1111  1111112222222222333333333344444
        #Coor     01234567890123  4567890123456789012345678901234
        #ref      AGCATGTTAGATAA**GATAGCTGTGCTAGTAGGCAGTCAGCGCCAT
        #               01234567890123 456
        #+r001/1        TAAGATAAAGGATA*CTG
        #              012345678 90123
        #+r002         aaaAGATAA*GGATA
        #                 012345
        #+r003       gcctaAGCTAA
        #                          012345              67890
        #+r004                     ATAGCT..............TCAGC
        #                                       43210
        #-r005                            ttagctTAGGC
        #                                               876543210
        #-r001/2                                        CAGCGCCAT
        sam = NamedTemporaryFile(suffix='.sam')
        sam.write(SAM)
        sam.flush()
        bam_fhand = NamedTemporaryFile()
        sam2bam(sam.name, bam_fhand.name)
        create_bam_index(bam_fhand.name)

        bam = pysam.Samfile(bam_fhand.name, "rb")

        reference = SeqWithQuality(seq=Seq(REF), name='ref')
        reference_id = reference.name

        for column in bam.pileup(reference=reference_id):
            for pileup_read in column.pileups:
                if pileup_read.alignment.qname == 'r001/1':

                    alignment = _get_alignment_section(pileup_read, 13, 21)
                    assert alignment[0] == 'N--NNNNNNNN'
                    assert alignment[1] == 'AAGGATA-CTG'

                    alignment = _get_alignment_section(pileup_read, 13, 21, reference)
                    assert alignment[0] == 'A--GATAGCTG'
                    assert alignment[1] == 'AAGGATA-CTG'

                    alignment = _get_alignment_section(pileup_read, 14, 21)
                    assert alignment[0] == 'NNNNNNNN'
                    assert alignment[1] == 'GATA-CTG'

                    alignment = _get_alignment_section(pileup_read, 14, 21, reference)
                    assert alignment[0] == 'GATAGCTG'
                    assert alignment[1] == 'GATA-CTG'

                    alignment = _get_alignment_section(pileup_read, 15, 21)
                    assert alignment[0] == 'NNNNNNN'
                    assert alignment[1] == 'ATA-CTG'

                    alignment = _get_alignment_section(pileup_read, 15, 21, reference)
                    assert alignment[0] == 'ATAGCTG'
                    assert alignment[1] == 'ATA-CTG'

                    alignment = _get_alignment_section(pileup_read, 5, 22)
                    assert alignment[0] == 'NNNNNNNNN--NNNNNNNNN'
                    assert alignment[1] == ' TAAGATAAAGGATA-CTG '

                    alignment = _get_alignment_section(pileup_read, 5, 22, reference)
                    assert alignment[0] == 'GTTAGATAA--GATAGCTGT'
                    assert alignment[1] == ' TAAGATAAAGGATA-CTG '

                    alignment = _get_alignment_section(pileup_read, 10, 16)
                    assert alignment[0] == 'NNNN--NNN'
                    assert alignment[1] == 'ATAAAGGAT'

                    alignment = _get_alignment_section(pileup_read, 10, 16, reference)
                    assert alignment[0] == 'ATAA--GAT'
                    assert alignment[1] == 'ATAAAGGAT'

                elif pileup_read.alignment.qname == 'r003':
                    alignment = _get_alignment_section(pileup_read, 5, 13)
                    assert alignment[0] == 'NNNNNNNNN'
                    assert alignment[1] == '   AGCTAA'

                    alignment = _get_alignment_section(pileup_read, 5, 13, reference)
                    assert alignment[0] == 'GTTAGATAA'
                    assert alignment[1] == '   AGCTAA'
                elif pileup_read.alignment.qname == 'r004':
                    alignment = _get_alignment_section(pileup_read, 18, 38)
                    assert alignment[0] == 'NNNNNNNNNNNNNNNNNNNNN'
                    assert alignment[1] == 'GCT--------------TCAG'

                    alignment = _get_alignment_section(pileup_read, 18, 38, reference)
                    assert alignment[0] == 'GCTGTGCTAGTAGGCAGTCAG'
                    assert alignment[1] == 'GCT--------------TCAG'
    @staticmethod
    def test_join_alignments():
        'It test that we can join single alignments into a multiple alignment'

        # check alignment without inserts
        # ref     TGTCGTATGTAGTGGTAGTCTAGTAGTA
        # READ1           GT--TGGT
        # READ2           GTAGTGGT
        alignments = {'read1':('GTAGTGGT', 'GT--TGGT'),
                      'read2':('GTAGTGGT', 'GTAGTGGT')}
        alig = _make_multiple_alignment(alignments)
        assert alig['reference'] == ['G', 'T', 'A', 'G', 'T', 'G', 'G', 'T']
        assert alig['reads']['read1'] == ['G', 'T', '-', '-', 'T', 'G', 'G',
                                          'T']
        # ref     TGTCGTATGTAGTG---GTAGTCTAGTAGTA
        # READ1           GT--TG---GT
        # READ2           GTAGTG---GT
        # read3           GTAGTGAA-GT
        # read4           GTAGTGAACGT
        alignments = {'read1':('GTAGTGGT', 'GT--TGGT'),
                      'read2':('GTAGTGGT', 'GTAGTGGT'),
                      'read3':('GTAGTG--GT', 'GTAGTGAAGT'),
                      'read4':('GTAGTG---GT', 'GTAGTGAACGT')}
        alignment = _make_multiple_alignment(alignments)
        assert alignment == {'reference': ['G','T','A','G','T','G','---','G','T'],
            'reads':{'read1': ['G', 'T', '-', '-', 'T', 'G', '---', 'G', 'T'],
                     'read2': ['G', 'T', 'A', 'G', 'T', 'G', '---', 'G', 'T'],
                     'read3': ['G', 'T', 'A', 'G', 'T', 'G', 'AA-', 'G', 'T'],
                     'read4': ['G', 'T', 'A', 'G', 'T', 'G', 'AAC', 'G', 'T']}}

class PoblationCalculationsTest(unittest.TestCase):
    'It checks the calculations of the poblations'

    @staticmethod
    def test_annotate_pic():
        'It tests the calculation of PIC(UMVU)'

        alleles = {('A', SNP): {'read_groups':{'rg1':1, 'rg2':1, 'rg4':2}},
                   ('T', INVARIANT): {'read_groups':{'rg1':1, 'rg3':2}}}
        snv = SeqFeature(type='snv', location=FeatureLocation(11, 11),
                         qualifiers={'alleles':alleles,
                                     'read_groups':{}})

        annotate_pic(snv)
        assert round(snv.qualifiers['pic'], 2) == 0.49

        alleles = {('T', INVARIANT): {'read_groups':{'rg1':1, 'rg3':2}}}
        snv = SeqFeature(type='snv', location=FeatureLocation(11, 11),
                         qualifiers={'alleles':alleles,
                                     'read_groups':{}})

        annotate_pic(snv)
        assert round(snv.qualifiers['pic'], 2) == 0

        alleles = {('A', SNP): {'read_groups':{'rg1':1}},
                   ('T', INVARIANT): {'read_groups':{'rg1':1, 'rg3':1}}}
        snv = SeqFeature(type='snv', location=FeatureLocation(11, 11),
                         qualifiers={'alleles':alleles,
                                     'read_groups':{}})

        annotate_pic(snv)
        assert snv.qualifiers['pic'] == None

    @staticmethod
    def test_annotate_heterozygosity():
        'It tests the calculation of heterozygosity'

        alleles = {('A', SNP): {'read_groups':{'rg1':1, 'rg2':1, 'rg4':2}},
                   ('T', INVARIANT): {'read_groups':{'rg1':1, 'rg3':2}}}
        snv = SeqFeature(type='snv', location=FeatureLocation(11, 11),
                         qualifiers={'alleles':alleles,
                                     'read_groups':{}})

        annotate_heterozygosity(snv, ploidy=2)
        assert round(snv.qualifiers['heterozygosity'], 2) == 0.53

        alleles = {('T', INVARIANT): {'read_groups':{'rg1':1, 'rg3':2}}}
        snv = SeqFeature(type='snv', location=FeatureLocation(11, 11),
                         qualifiers={'alleles':alleles,
                                     'read_groups':{}})

        annotate_heterozygosity(snv, ploidy=2)
        assert round(snv.qualifiers['heterozygosity'], 2) == 0

        alleles = {('A', SNP): {'read_groups':{'rg1':100, 'rg2':150}},
                   ('T', INVARIANT): {'read_groups':{'rg1':50, 'rg3':100}}}
        snv = SeqFeature(type='snv', location=FeatureLocation(11, 11),
                         qualifiers={'alleles':alleles,
                                     'read_groups':{}})
        annotate_heterozygosity(snv, ploidy=2)
        assert round(snv.qualifiers['heterozygosity'], 2) == 0.47

if __name__ == "__main__":
    import sys;sys.argv = ['', 'TestSnvAnnotation.test_snv_annotation_massive']
#    import sys;sys.argv = ['', 'TestReadPos.test_join_alignments']
    unittest.main()
