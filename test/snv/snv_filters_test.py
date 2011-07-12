'''
Created on 19/02/2010

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

import unittest, os

from Bio.SeqFeature import FeatureLocation

from franklin.utils.misc_utils import TEST_DATA_DIR
from franklin.seq.seqs import SeqWithQuality, Seq, SeqFeature
from franklin.snv.snv_annotation import (INVARIANT, SNP, INDEL, DELETION,
                                         INSERTION)
from franklin.snv.snv_filters import (create_unique_contiguous_region_filter,
                                      create_close_to_intron_filter,
                                      create_high_variable_region_filter,
                                      create_close_to_snv_filter,
                                      create_snv_close_to_limit_filter,
                                      create_major_allele_freq_filter,
                                      create_kind_filter,
                                      create_cap_enzyme_filter,
                                      create_is_variable_filter,
                                      create_reference_in_list_filter,
                                      create_not_variable_in_group_filter,
                                      create_min_groups_filter,
                                      SnvNamer,
                                      create_in_segment_filter)

class SeqVariationFilteringTest(unittest.TestCase):
    'It checks the filtering methods.'

    @staticmethod
    def test_ref_in_list_filter():
        'We filter out the snv close to an intron'
        snv = SeqFeature(type='snv', location=FeatureLocation(100, 100),
                          qualifiers={})
        seq1 = SeqWithQuality(name='seq1', seq=Seq('A'), features=[snv])
        snv1 = SeqFeature(type='snv', location=FeatureLocation(100, 100),
                          qualifiers={})
        seq2 = SeqWithQuality(name='seq2', seq=Seq('A'), features=[snv1])
        seq_list = ['seq1']
        filter_ = create_reference_in_list_filter(seq_list)

        filter_(seq1)
        for snv, expected in zip(seq1.get_features(kind='snv'), [True]):
            result = snv.qualifiers['filters']['ref_not_in_list'][None]
            assert result == expected

        filter_(seq2)
        for snv, expected in zip(seq2.get_features(kind='snv'), [False]):
            result = snv.qualifiers['filters']['ref_not_in_list'][None]
            assert result == expected

    @staticmethod
    def test_unique_contiguous_region():
        '''It tests that we remove snv located in regions that are not unique in
        the genome'''

        filter_id = 'uniq_contiguous'
        arabidopsis_genes = 'arabidopsis_genes+'

        genomic_db = os.path.join(TEST_DATA_DIR, 'blast', arabidopsis_genes)
        distance = 60
        filter_ = create_unique_contiguous_region_filter(distance=distance,
                                                         genomic_db=genomic_db,
                                            genomic_seqs_fpath=genomic_db)
        #an snv in an unique region
        seq = 'CTGGAATCTCTGAGTTTCTGGGTTCAAGTTGCACTGACCATTGTTGGATTTGTAGATTGTTTC'
        seq += 'TTCATTTCATTAGGCATTGATTATGGGTAAATGCGTGGGTACATATAATATATATCTGTTGAA'
        seq += 'TGCAATTTACACATTGACTGAGGAACAACATGAACATGGCAGCTTTCTCAAAATTGAACCACA'
        seq += 'GAAGGCTTAAAAGCAAAGTCTTTGGAGAATCAGACTAAGCTTGAGAAAGGGTTGAAATTTATC'
        seq += 'TATCTAAGACTTGTGAAAAGAAACTGTATTTAATGCTTTCACATTGCTGGTTTTTGAGGGACA'
        seq += 'TACATCATTGAAAGTTGAAACTGATAGTAGCAGAGTTTTTTCCTCTGTTTGGTTCTCTTTTGT'
        seq += 'TCTTTTAGTTCGTAGTGATAGCTTAATGCGTAGCTGAGAAAAGTGGAATGCTGTTTGTTTATA'
        seq += 'ATACTAAAGTAGGTCATTCTGCATTTCTTGCAATGTCAGTTGTA'
        alleles = {('A', SNP): None, ('T', INVARIANT):None}
        snv1 = SeqFeature(type='snv', location=FeatureLocation(50, 50),
                          qualifiers={'alleles':alleles})
        seq = SeqWithQuality(seq=Seq(seq), features=[snv1])
        filter_(seq)
        assert not seq.features[0].qualifiers['filters'][filter_id][distance]

        #an snv in a region with two matches
        seq = 'CCACTACAAGAGGTGGAAGAGCGAAAACTCTGTTTATTACTAGCTAGGGTTTCTATTAATGAA'
        seq += 'AGGTTCATGTAAATATATGAAGATGGGAAGCAAGAGGTGTTCAAGGAGAAGAGGGAGTTAGAC'
        seq += 'GACCAGAAGATGGCGGTGACGTTTAGAGGGCTGGATGGTCATGTGATGGAGCAGCTTAAGGTG'
        seq += 'TATGACGTCATCTTCCAATTCGTCCCTAAGTCTCAGGAAGGTTGCGTCTGCAAAGTCACTATG'
        seq += 'TTTTGGGAGAAGCGCTA'
        alleles = {('A', SNP): None, ('T', INVARIANT):None}
        snv1 = SeqFeature(type='snv', location=FeatureLocation(100, 100),
                          qualifiers={'alleles':alleles})
        seq = SeqWithQuality(seq=Seq(seq), features=[snv1])
        seq = filter_(seq)
        assert seq.features[0].qualifiers['filters'][filter_id][distance]

        #a sequence with one hit but two hsps, but a contiguous region according
        #to est2genome
        seq = 'TAACATATGTATCGTTTGCTAACTGTATATCAGGAGAAGGAGTAATGTAAAAATTCGAGA'
        seq += 'TTATTGTAATTTAAGAAACTCTTTAGGTAGAATGTAATTAATACCAAAATCATTAAAATT'
        seq += 'TCTGTACTTTATACTTGTCATATTCTCAAGATACGAGTCAAAAATGTCTGATTAAATATG'
        seq += 'AACATATTTTTATATATTCATTCCATCATAATCTTCGAGATTTAAAAAGCTCTGATTATC'
        alleles = {('A', SNP): None, ('T', INVARIANT):None}
        snv1 = SeqFeature(type='snv', location=FeatureLocation(50, 50),
                          qualifiers={'alleles':alleles})
        seq = SeqWithQuality(seq=Seq(seq), features=[snv1])
        seq = filter_(seq)
        seq = filter_(seq)
        assert not seq.features[0].qualifiers['filters'][filter_id][distance]

    @staticmethod
    def test_close_to_intron_filter():
        'We filter out the snv close to an intron'
        intron = SeqFeature(location=FeatureLocation(478, 478), type='intron')

        snv1 = SeqFeature(type='snv', location=FeatureLocation(100, 100),
                          qualifiers={})
        snv2 = SeqFeature(type='snv', location=FeatureLocation(450, 450),
                          qualifiers={})
        snv3 = SeqFeature(type='snv', location=FeatureLocation(640, 640),
                          qualifiers={})
        snv4 = SeqFeature(type='snv', location=FeatureLocation(700, 700),
                          qualifiers={})

        seq = SeqWithQuality(seq=Seq('A' * 1000), features=[intron, snv1, snv2,
                                                            snv3, snv4])

        filter_ = create_close_to_intron_filter(distance=60)
        filter_(seq)
        filter_(seq)

        for snv, expected in zip(seq.get_features(kind='snv'),
                                 [False, True, False, False]):
            result = snv.qualifiers['filters']['close_to_intron'][60]
            assert result == expected

    @staticmethod
    def test_high_variable_region_filter():
        'It test  high_variable_region_filter'
        snv1 = SeqFeature(type='snv', location=FeatureLocation(1, 1),
                          qualifiers={})
        snv2 = SeqFeature(type='snv', location=FeatureLocation(4, 4),
                          qualifiers={})
        snv3 = SeqFeature(type='snv', location=FeatureLocation(6, 6),
                          qualifiers={})

        seq_str = 'AATATA'
        seq = SeqWithQuality(seq=Seq(seq_str), qual=[30] * len(seq_str),
                             features=[snv1, snv2, snv3])
        max_variability = 0.4
        filter_ = create_high_variable_region_filter(max_variability)
        filter_(seq)
        threshold = (max_variability, None)
        for snv, expected in zip(seq.get_features(kind='snv'),
                                 [True, True, True]):
            result = snv.qualifiers['filters']['high_variable_reg'][threshold]
            assert result == expected

        max_variability = 0.6
        filter_ = create_high_variable_region_filter(max_variability)
        filter_(seq)
        filter_(seq)
        threshold = (max_variability, None)
        for snv, expected in zip(seq.get_features(kind='snv'),
                                 [False, False, False]):
            result = snv.qualifiers['filters']['high_variable_reg'][threshold]
            assert result == expected
        max_variability = 0.25
        window = 6
        threshold = (max_variability, window)
        filter_ = create_high_variable_region_filter(max_variability,
                                                       window=window)
        filter_(seq)
        for snv, expected in zip(seq.get_features(kind='snv'),
                                 [False, True, False]):
            result = snv.qualifiers['filters']['high_variable_reg'][threshold]
            assert result == expected

    @staticmethod
    def test_close_to_seqvar_filter():
        'It tests that we can detect snvs by its proximity to another snv'
        alleles = {('A', SNP): None, ('T', INVARIANT):None}
        snv1 = SeqFeature(type='snv', location=FeatureLocation(1, 1),
                          qualifiers={'alleles':alleles})
        snv2 = SeqFeature(type='snv', location=FeatureLocation(4, 4),
                          qualifiers={'alleles':alleles})
        snv3 = SeqFeature(type='snv', location=FeatureLocation(6, 6),
                          qualifiers={'alleles':alleles})
        seq_str = 'AATATA'
        seq = SeqWithQuality(seq=Seq(seq_str), qual=[30] * len(seq_str),
                             features=[snv1, snv2, snv3])
        proximity = 3

        filter_ = create_close_to_snv_filter(proximity)
        filter_(seq)
        for snv, expected in zip(seq.get_features(kind='snv'),
                                 [False, True, False]):
            result = snv.qualifiers['filters']['close_to_snv'][(proximity, None, None)]
            assert result == expected

        snv.qualifiers['filters']['close_to_snv']
        alleles2 = {('A', DELETION): None, ('T', INVARIANT):None}
        snv1 = SeqFeature(type='snv', location=FeatureLocation(1, 1),
                          qualifiers={'alleles':alleles2})
        snv2 = SeqFeature(type='snv', location=FeatureLocation(4, 4),
                          qualifiers={'alleles':alleles2})
        snv3 = SeqFeature(type='snv', location=FeatureLocation(6, 6),
                          qualifiers={'alleles':alleles2})
        alleles3 = {('A', INSERTION): None, ('T', INVARIANT):None}
        snv4 = SeqFeature(type='snv', location=FeatureLocation(9, 9),
                          qualifiers={'alleles':alleles3})
        seq = SeqWithQuality(seq=Seq(seq_str), qual=[30] * len(seq_str),
                             features=[snv1, snv2, snv3, snv4])
        filter_ = create_close_to_snv_filter(proximity, INDEL)
        filter_(seq)
        for snv, expected in zip(seq.get_features(kind='snv'),
                                 [False, True, True, False]):
            result = snv.qualifiers['filters']['close_to_snv'][(proximity, INDEL, None)]
            assert result == expected


    @staticmethod
    def test_close_to_limit_filter():
        'It tests that we can detect snvs close to the limit'
        snv1 = SeqFeature(type='snv', location=FeatureLocation(1, 1),
                          qualifiers={})
        snv2 = SeqFeature(type='snv', location=FeatureLocation(4, 4),
                          qualifiers={})
        snv3 = SeqFeature(type='snv', location=FeatureLocation(6, 6),
                          qualifiers={})
        seq_str = 'AATATA'
        seq = SeqWithQuality(seq=Seq(seq_str), qual=[30] * len(seq_str),
                             features=[snv1, snv2, snv3])
        distance = 2

        filter_ = create_snv_close_to_limit_filter(distance)
        filter_(seq)
        for snv, expected in zip(seq.get_features(kind='snv'),
                                 [True, False, True]):
            result = snv.qualifiers['filters']['close_to_limit'][distance]
            assert result == expected

    @staticmethod
    def test_major_allele_freq_filter_snv():
        'It test the first allele percent filter'
        read_groups = {'g1':{}, 'g2':{}}
        alleles = {('A', INVARIANT):{'read_groups':{'g1':4}},
                   ('T', SNP)      :{'read_groups':{'g2':2}}}
        snv1 = SeqFeature(type='snv', location=FeatureLocation(1, 1),
                          qualifiers={'alleles':alleles,
                                      'read_groups':read_groups})
        alleles = {('A', INVARIANT):{'read_groups':{'g1':3}},
                   ('T', SNP)      :{'read_groups':{'g2':2}}}

        snv2 = SeqFeature(type='snv', location=FeatureLocation(3, 3),
                          qualifiers={'alleles':alleles,
                                      'read_groups':read_groups})

        seq_str = 'AATATA'
        seq = SeqWithQuality(seq=Seq(seq_str), qual=[30] * len(seq_str),
                             features=[snv1, snv2])
        frecuency = 0.59999999999999998
        filter_ = create_major_allele_freq_filter(frecuency)
        filter_(seq)
        for snv, expected in zip(seq.get_features(kind='snv'), [True, False]):
            result = snv.qualifiers['filters']['maf'][(frecuency, )]
            assert result == expected

        #now we do it only for one read group
        snv1 = SeqFeature(type='snv', location=FeatureLocation(1, 1),
                          qualifiers={'alleles':alleles,
                                      'read_groups':read_groups})
        snv2 = SeqFeature(type='snv', location=FeatureLocation(3, 3),
                          qualifiers={'alleles':alleles,
                                      'read_groups':read_groups})
        seq = SeqWithQuality(seq=Seq(seq_str), qual=[30] * len(seq_str),
                             features=[snv1, snv2])
        frecuency = 0.59999999999999998
        filter_ = create_major_allele_freq_filter(frecuency,
                                                   groups=['g1'],
                                                   group_kind='read_groups')
        filter_(seq)
        for snv, expected in zip(seq.get_features(kind='snv'), [True, True]):
            parameters = (0.59999999999999998, ('g1',), 'read_groups')
            result = snv.qualifiers['filters']['maf'][parameters]
            assert result == expected

    @staticmethod
    def test_kind_filter():
        'It test the kind filter'
        alleles = {('A', INVARIANT):{},
                   ('T', SNP)      :{}}

        snv1 = SeqFeature(type='snv', location=FeatureLocation(1, 1),
                          qualifiers={'alleles':alleles})
        alleles = {('A', INVARIANT):{},
                   ('T', INDEL)    :{}}

        snv2 = SeqFeature(type='snv', location=FeatureLocation(3, 3),
                          qualifiers={'alleles':alleles})

        seq_str = 'AATATA'
        seq = SeqWithQuality(seq=Seq(seq_str), qual=[30] * len(seq_str),
                             features=[snv1, snv2])
        kind = SNP
        filter_ = create_kind_filter(kind)
        filter_(seq)

        for snv, expected in zip(seq.get_features(kind='snv'), [False, True]):
            result = snv.qualifiers['filters']['by_kind'][kind]
            assert result == expected

    @staticmethod
    def test_cap_enzyme_filter():
        'It test the cap enzyme filter'
        seq = 'ATGATGATG' + 'gaaattc' + 'ATGATGATGTGGGAT'

        alleles = {('A', INVARIANT):{},
                   ('A', DELETION) :{}}
        snv = SeqFeature(type='snv', location=FeatureLocation(11, 11),
                         qualifiers={'alleles':alleles})
        seq = SeqWithQuality(seq=Seq(seq), name='ref', features=[snv])

        all_enzymes = True
        filter_ = create_cap_enzyme_filter(all_enzymes)
        filter_(seq)
        for snv, expected in zip(seq.get_features(kind='snv'), [False]):
            result = snv.qualifiers['filters']['cap_enzymes'][all_enzymes]
            assert result == expected

        #No cap
        seq = 'ATGATGATG' + 'ATGATGATGTGGGAT'

        alleles = {('A', INVARIANT):{},
                   ('A', DELETION) :{}}
        snv = SeqFeature(type='snv', location=FeatureLocation(11, 11),
                         qualifiers={'alleles':alleles})
        seq = SeqWithQuality(seq=Seq(seq), name='ref', features=[snv])

        all_enzymes = True
        filter_ = create_cap_enzyme_filter(all_enzymes)
        filter_(seq)
        for snv, expected in zip(seq.get_features(kind='snv'), [True]):
            result = snv.qualifiers['filters']['cap_enzymes'][all_enzymes]
            assert result == expected

    @staticmethod
    def test_is_not_variable_filter():
        'It tests variable filter function'
        alleles = {('A', SNP): {'read_groups':{'rg1':1, 'rg2':1, 'rg4':2}},
                   ('T', INVARIANT): {'read_groups':{'rg1':1, 'rg3':2}}}
        snv = SeqFeature(type='snv', location=FeatureLocation(11, 11),
                         qualifiers={'alleles':alleles,
                                     'read_groups':{}})
        seq = 'ATGATGATGgaaattcATGATGATGTGGGAT'

        alleles2 = {('A', SNP): {'read_groups':{'rg1':2}}}
        snv2 = SeqFeature(type='snv', location=FeatureLocation(21, 21),
                         qualifiers={'alleles':alleles2,
                                     'read_groups':{}})
        seq = SeqWithQuality(seq=Seq(seq), name='ref', features=[snv, snv2])

        filter_ = create_not_variable_in_group_filter(group_kind= 'read_groups',
                                                      groups=['rg1'],
                                                      in_union=True)
        filter_(seq)

        filter_ = create_not_variable_in_group_filter(group_kind= 'read_groups',
                                                      groups=['rg1'],
                                                      in_union=True,
                                                      reference_free=False)
        filter_(seq)

        filter_ = create_not_variable_in_group_filter(group_kind= 'read_groups',
                                                      groups=['rg1'],
                                                      in_union=False,
                                                      reference_free=False)
        filter_(seq)

        filter_ = create_not_variable_in_group_filter(group_kind= 'read_groups',
                                                      groups=['rg2', 'rg4'],
                                                      in_union=True)
        filter_(seq)

        filter_ = create_not_variable_in_group_filter(group_kind= 'read_groups',
                                                      groups=['rg3', 'rg4'],
                                                      in_union=True)
        filter_(seq)

        filter_ = create_not_variable_in_group_filter(group_kind= 'read_groups',
                                                      groups=['rg3', 'rg4'],
                                                      in_union=False)
        filter_(seq)

        filter_ = create_not_variable_in_group_filter(group_kind= 'read_groups',
                                                      groups=['rg5'],
                                                      in_union=True)
        filter_(seq)

        snv_1 = seq.features[0].qualifiers['filters']['is_not_variable']
        snv_2 = seq.features[1].qualifiers['filters']['is_not_variable']

        assert snv_1['read_groups', ('rg1',), True, True, None, None]
        assert snv_1['read_groups', ('rg1',), True, False, None, None]
        assert snv_1['read_groups', ('rg1',), False, False, None, None]
        assert not snv_2['read_groups', ('rg1',), True, True, None, None]


        assert not snv_1['read_groups', ('rg2', 'rg4'), True, True, None, None]
        assert snv_2['read_groups', ('rg2', 'rg4'), True, True, None, None]


        assert  snv_1['read_groups', ('rg3', 'rg4'), True, True, None, None]
        assert  not snv_1['read_groups', ('rg3', 'rg4'), False, True, None, None]

        assert snv_1['read_groups', ('rg5',), True, True, None, None]

        alleles = {('A', SNP): {'read_groups':{'rg1':10}},
                   ('T', INVARIANT): {'read_groups':{'rg1':90}}}
        snv = SeqFeature(type='snv', location=FeatureLocation(11, 11),
                         qualifiers={'alleles':alleles,
                                     'read_groups':{}})

        seq = 'ATGATGATGgaaattcATGATGATGTGGGAT'
        seq = SeqWithQuality(seq=Seq(seq), name='ref', features=[snv])
        maf = 0.95
        min_num_reads = 50

        filter_ = create_not_variable_in_group_filter(group_kind= 'read_groups',
                                                      groups=['rg1'],
                                                      in_union=True,
                                                      maf=maf,
                                                    min_num_reads=min_num_reads)
        filter_(seq)

        maf = 0.6
        min_num_reads = 50

        filter_ = create_not_variable_in_group_filter(group_kind= 'read_groups',
                                                      groups=['rg1'],
                                                      in_union=True,
                                                      maf=maf,
                                                    min_num_reads=min_num_reads)
        filter_(seq)

        maf = 0.95
        min_num_reads = 200

        filter_ = create_not_variable_in_group_filter(group_kind= 'read_groups',
                                                      groups=['rg1'],
                                                      in_union=True,
                                                      maf=maf,
                                                    min_num_reads=min_num_reads)
        filter_(seq)

        maf = 0.6
        min_num_reads = 200

        filter_ = create_not_variable_in_group_filter(group_kind= 'read_groups',
                                                      groups=['rg1'],
                                                      in_union=True,
                                                      maf=maf,
                                                    min_num_reads=min_num_reads)
        filter_(seq)

        snv_1 = seq.features[0].qualifiers['filters']['is_not_variable']

        assert snv_1['read_groups', ('rg1',), True, True, 0.95, 50]
        assert not snv_1['read_groups', ('rg1',), True, True, 0.6, 50]
        assert not snv_1['read_groups', ('rg1',), True, True, 0.95, 200]
        assert not snv_1['read_groups', ('rg1',), True, True, 0.6, 200]


        alleles = {('A', SNP): {'read_groups':{'rg1':1, 'rg2':1, 'rg4':2}},
                   ('T', INVARIANT): {'read_groups':{'rg1':1, 'rg3':2}}}
        snv = SeqFeature(type='snv', location=FeatureLocation(11, 11),
                         qualifiers={'alleles':alleles,
                                     'read_groups':{'rg1':{'SM':'sm1'},
                                                    'rg2':{'SM':'sm2'},
                                                    'rg3':{'SM':'sm3'},
                                                    'rg4':{'SM':'sm4'}}})

        seq = 'ATGATGATGgaaattcATGATGATGTGGGAT'
        seq = SeqWithQuality(seq=Seq(seq), name='ref', features=[snv])

        filter_ = create_not_variable_in_group_filter(group_kind= 'samples',
                                                      groups=['sm1'],
                                                      in_union=True)
        filter_(seq)
        snv_1 = seq.features[0].qualifiers['filters']['is_not_variable']
        assert snv_1['samples', ('sm1',), True, True, None, None]

    @staticmethod
    def test_is_variable_filter():
        'It tests variable filter function'
        alleles = {('A', SNP): {'read_groups':{'rg1':1, 'rg2':2, 'rg4':2}},
                   ('T', INVARIANT): {'read_groups':{'rg1':1, 'rg3':2}}}
        snv = SeqFeature(type='snv', location=FeatureLocation(11, 11),
                         qualifiers={'alleles':alleles,
                                     'read_groups':{}})
        seq = 'ATGATGATGgaaattcATGATGATGTGGGAT'
        seq = SeqWithQuality(seq=Seq(seq), name='ref', features=[snv])

        alleles2 = {('A', SNP): {'read_groups':{'rg1':2}}}
        snv2 = SeqFeature(type='snv', location=FeatureLocation(11, 11),
                         qualifiers={'alleles':alleles2,
                                     'read_groups':{}})
        seq2 = 'ATGATGATGgaaattcATGATGATGTGGGAT'
        seq2 = SeqWithQuality(seq=Seq(seq2), name='ref2', features=[snv2])

        filters = []
        parameters = []
        results = []
        reference_free = True
        maf= None
        in_all_groups = True
        min_num_reads = None
        min_reads_per_allele = None


        kind = 'read_groups'
        groups = ('rg1',)
        in_union = False
        params = (kind, groups, in_union, in_all_groups, reference_free, maf,
                  min_num_reads, min_reads_per_allele)
        parameters.append(params)
        filter_ = create_is_variable_filter(*params)
        filters.append(filter_)
        results.append(False)

        kind = 'read_groups'
        groups = ('rg1',)
        in_union = True
        params = (kind, groups, in_union, in_all_groups, reference_free, maf,
                  min_num_reads, min_reads_per_allele)
        parameters.append(params)
        filter_ = create_is_variable_filter(*params)
        filters.append(filter_)
        results.append(False)

        kind = 'read_groups'
        groups = ('rg1',)
        in_union = True
        min_reads_per_allele = 2
        params = (kind, groups, in_union, in_all_groups, reference_free, maf,
                  min_num_reads, min_reads_per_allele)
        parameters.append(params)
        filter_ = create_is_variable_filter(*params)
        filters.append(filter_)
        results.append(True)

        min_reads_per_allele = None


        kind = 'read_groups'
        groups = ('rg2', 'rg4')
        in_union = True
        params = (kind, groups, in_union, in_all_groups, reference_free, maf,
                  min_num_reads, min_reads_per_allele)
        parameters.append(params)
        filters.append(create_is_variable_filter(*params))
        results.append(True)

        kind = 'read_groups'
        groups = 'fake'
        in_union = True
        params = (kind, (groups,), in_union, in_all_groups, reference_free, maf,
                  min_num_reads, min_reads_per_allele)
        parameters.append(params)
        filters.append(create_is_variable_filter(*params))
        results.append(True)

        kind = 'read_groups'
        groups = ('rg2', 'rg3')
        in_union = False
        params = (kind, groups, in_union, in_all_groups, reference_free, maf,
                  min_num_reads, min_reads_per_allele)
        parameters.append(params)
        filters.append(create_is_variable_filter(*params))
        results.append(True)

        kind = 'read_groups'
        groups = ('rg2', 'rg3')
        in_union = True
        params = (kind, groups, in_union, in_all_groups, reference_free, maf,
                  min_num_reads, min_reads_per_allele)
        parameters.append(params)
        filters.append(create_is_variable_filter(*params))
        results.append(False)

        kind = 'read_groups'
        groups = ('rg5',)
        in_union = True
        params = (kind, groups, in_union, in_all_groups, reference_free, maf,
                  min_num_reads, min_reads_per_allele)
        parameters.append(params)
        filters.append(create_is_variable_filter(*params))
        results.append(True)

        kind = 'read_groups'
        groups = ('rg2',)
        in_union = True
        reference_free = False
        params = (kind, groups, in_union, in_all_groups, reference_free, maf,
                  min_num_reads, min_reads_per_allele)
        parameters.append(params)
        filters.append(create_is_variable_filter(*params))
        results.append(False)

        kind = 'read_groups'
        groups = ('rg2',)
        in_union = True
        reference_free = True
        params = (kind, groups, in_union, in_all_groups, reference_free, maf,
                  min_num_reads, min_reads_per_allele)
        parameters.append(params)
        filters.append(create_is_variable_filter(*params))
        results.append(True)


        for filter_ in filters:
            filter_(seq)
            filter_(seq2)

        for params, expected in zip(parameters, results):
            for snv, expected in zip(seq.get_features(kind='snv'), [expected]):
                result = snv.qualifiers['filters']['is_variable'][params]
                assert result == expected

        for params in parameters:
            for snv in seq2.get_features(kind='snv'):
                assert snv.qualifiers['filters']['is_variable'][params]


    @staticmethod
    def test_min_num_groups_filter():
        'It tests the min num groups svn filter'
        alleles = {('A', SNP): {'read_groups':{'rg1':1, 'rg2':1, 'rg4':2}},
                   ('T', INVARIANT): {'read_groups':{'rg1':1, 'rg3':2}}}
        snv = SeqFeature(type='snv', location=FeatureLocation(11, 11),
                         qualifiers={'alleles':alleles,
                                     'read_groups':{'rg1':{'LB':1},
                                                    'rg2':{'LB':2},
                                                    'rg3':{'LB':3},
                                                    'rg4':{'LB':4}}})
        seq = 'ATGATGATGgaaattcATGATGATGTGGGAT'
        seq = SeqWithQuality(seq=Seq(seq), name='ref', features=[snv])

        filter_ = create_min_groups_filter(min_groups=4)
        filter_(seq)
        filter_ = create_min_groups_filter(min_groups=5)
        filter_(seq)
        filter_ = create_min_groups_filter(min_groups=4, group_kind='libraries')
        filter_(seq)
        results = seq.features[0].qualifiers['filters']['min_groups']
        assert not results[(4, 'read_groups')]
        assert results[(5, 'read_groups')]
        assert not results[(4, 'libraries')]

    @staticmethod
    def test_get_filter_description():
        'It tests get_filter_description function'

        filter_name = 'close_to_intron'
        parameters = 30
        filter_descriptions = {}
        namer = SnvNamer()
        name, desc = namer.get_filter_description(filter_name, parameters,
                                            filter_descriptions)
        assert name == 'I30'
        assert desc == 'An intron is located closer than 30 base pairs'

        filter_name = 'maf'
        parameters = 0.6
        filter_descriptions = {}
        name, desc = namer.get_filter_description(filter_name, parameters,
                                            filter_descriptions)
        assert name == 'maf1'
        assert desc == 'The most frequent allele in All: All. frequency greater than 0.60'

        filter_name = 'maf'
        parameters = (0.6, ('1', '2'), 'read_group')
        filter_descriptions = {}
        name, desc = namer.get_filter_description(filter_name, parameters,
                                            filter_descriptions)
        assert name == 'maf2'
        assert desc == 'The most frequent allele in read_group: 1,2. frequency greater than 0.60'

        filter_name = 'by_kind'
        parameters = SNP
        filter_descriptions = {}
        name, desc = namer.get_filter_description(filter_name, parameters,
                                            filter_descriptions)

        filter_name = 'is_variable'
        kind = 'read_groups'
        groups = ['rg1', 'rg2']
        maf = 0.6
        min_num_reads = 3
        in_union = True
        in_all_groups = True
        reference_free = True
        parameters = (kind, tuple(groups), in_union, in_all_groups,
                      reference_free, maf, min_num_reads)
        filter_descriptions = {}
        name, desc = namer.get_filter_description(filter_name, parameters,
                                            filter_descriptions)
        assert name[:3] == 'vrg'
        descrip = "It is not variable, or no data, in the read_groups : rg1,rg2."
        descrip += ' All together: True. maf:0.600000. min_num_reads:3'
        assert desc == descrip

        name, desc = namer.get_filter_description(filter_name, parameters,
                                            filter_descriptions)
        assert name[:3] == 'vrg'

        filter_name = 'is_not_variable'
        kind = 'read_groups'
        groups = ['rg1', 'rg2']
        maf = 0.7
        min_num_reads = 2
        in_union = True
        in_all_groups = True
        reference_free = True
        parameters = (kind, tuple(groups), in_union, in_all_groups,
                      reference_free, maf, min_num_reads)
        filter_descriptions = {}
        name, desc = namer.get_filter_description(filter_name, parameters,
                                            filter_descriptions)
        assert name[:4] == 'nvrg'
        descrip = "It is variable, or no data, in the read_groups : rg1,rg2."
        descrip += ' All together: True. maf:0.700000. min_num_reads:2'
        assert desc == descrip

        filter_name = 'is_not_variable'
        kind = 'read_groups'
        groups = ['rg1', 'rg2']
        maf = None
        min_num_reads = None
        in_union = True
        in_all_groups = True
        reference_free = True
        parameters = (kind, tuple(groups), in_union, in_all_groups,
                      reference_free, maf, min_num_reads)
        filter_descriptions = {}
        name, desc = namer.get_filter_description(filter_name, parameters,
                                            filter_descriptions)
        assert name[:4] == 'nvrg'
        descrip = "It is variable, or no data, in the read_groups : rg1,rg2."
        descrip += ' All together: True'
        assert desc == descrip

        parameters = (4, 'read_groups')
        name, desc = namer.get_filter_description('min_groups',
                                            parameters,
                                            filter_descriptions)
        assert name == 'mr4'
        assert desc == 'SNV read in less than 4 read_groups'

        parameters = (True)
        name, desc = namer.get_filter_description('cap_enzymes',
                                            parameters,
                                            filter_descriptions)
        assert name == 'cet'
        assert desc == 'SNV is not a CAP detectable by the enzymes: all'

        parameters = (False)
        name, desc = namer.get_filter_description('cap_enzymes',
                                            parameters,
                                            filter_descriptions)
        assert name == 'cef'
        assert desc == 'SNV is not a CAP detectable by the enzymes: cheap ones'

    @staticmethod
    def test_create_in_segment_filter():
        'It tests create_in_segment_filter function'

        snv = SeqFeature(type='snv', location=FeatureLocation(0, 0),
                         qualifiers={})
        snv2 = SeqFeature(type='snv', location=FeatureLocation(11, 11),
                          qualifiers={})
        snv3 = SeqFeature(type='snv', location=FeatureLocation(25, 25),
                          qualifiers={})
        snv4 = SeqFeature(type='snv', location=FeatureLocation(35, 35),
                          qualifiers={})
        snv5 = SeqFeature(type='snv', location=FeatureLocation(50, 75),
                          qualifiers={})
        snv6 = SeqFeature(type='snv', location=FeatureLocation(65, 65),
                          qualifiers={})
        snv7 = SeqFeature(type='snv', location=FeatureLocation(75, 75),
                          qualifiers={})
        snv8 = SeqFeature(type='snv', location=FeatureLocation(80, 80),
                          qualifiers={})
        snv9 = SeqFeature(type='snv', location=FeatureLocation(90, 110),
                          qualifiers={})
        snv10 = SeqFeature(type='snv', location=FeatureLocation(110, 110),
                           qualifiers={})
        snv11 = SeqFeature(type='snv', location=FeatureLocation(125, 125),
                           qualifiers={})
        snv12 = SeqFeature(type='snv', location=FeatureLocation(125, 126),
                           qualifiers={})

        seq1 = 'ATGATGATGgaaattcATGATGATGTGGGAT'
        seq = SeqWithQuality(seq=Seq(seq1), name='seq1',
                             features=[snv, snv2, snv3, snv4, snv5, snv6, snv7,
                                       snv8, snv9, snv10, snv11, snv12])

        segments = {'seq1':[(5,25), (50, 75), (100, 125)]}

        filter_ = create_in_segment_filter(segments, edge_avoidance=None)
        filter_(seq)

        assert seq.features[0].qualifiers['filters']['in_segment'][0]
        assert not seq.features[1].qualifiers['filters']['in_segment'][0]
        assert not seq.features[2].qualifiers['filters']['in_segment'][0]
        assert seq.features[3].qualifiers['filters']['in_segment'][0]
        assert not seq.features[4].qualifiers['filters']['in_segment'][0]
        assert not seq.features[5].qualifiers['filters']['in_segment'][0]
        assert not seq.features[6].qualifiers['filters']['in_segment'][0]
        assert seq.features[7].qualifiers['filters']['in_segment'][0]
        assert seq.features[8].qualifiers['filters']['in_segment'][0]
        assert not seq.features[9].qualifiers['filters']['in_segment'][0]
        assert not seq.features[10].qualifiers['filters']['in_segment'][0]
        assert seq.features[11].qualifiers['filters']['in_segment'][0]

        for i in range(12):
            del seq.features[i].qualifiers['filters']

        filter_ = create_in_segment_filter(segments, edge_avoidance=3)
        filter_(seq)

        assert seq.features[0].qualifiers['filters']['in_segment'][3]
        assert not seq.features[1].qualifiers['filters']['in_segment'][3]
        assert seq.features[2].qualifiers['filters']['in_segment'][3]
        assert seq.features[3].qualifiers['filters']['in_segment'][3]
        assert seq.features[4].qualifiers['filters']['in_segment'][3]
        assert seq.features[5].qualifiers['filters']['in_segment'][3]
        assert seq.features[6].qualifiers['filters']['in_segment'][3]
        assert seq.features[7].qualifiers['filters']['in_segment'][3]
        assert seq.features[8].qualifiers['filters']['in_segment'][3]
        assert not seq.features[9].qualifiers['filters']['in_segment'][3]
        assert seq.features[10].qualifiers['filters']['in_segment'][3]
        assert seq.features[11].qualifiers['filters']['in_segment'][3]


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'SeqVariationFilteringTest.test_is_not_variable_filter']
    unittest.main()
