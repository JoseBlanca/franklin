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

from franklin.utils.misc_utils import DATA_DIR
from franklin.seq.seqs import SeqWithQuality, Seq, SeqFeature
from franklin.snv.snv_annotation import INVARIANT, SNP
from franklin.snv.snv_filters import (create_unique_contiguous_region_filter,
                                      create_close_to_intron_filter,
                                      create_high_variable_region_filter)

class SeqVariationFilteringTest(unittest.TestCase):
    'It checks the filtering methods.'

    @staticmethod
    def test_unique_contiguous_region():
        '''It tests that we remove snv located in regions that are not unique in
        the genome'''

        filter_id = 'uniq_contiguous'
        genomic_db = os.path.join(DATA_DIR, 'blast', 'arabidopsis_genes')
        distance = 60
        filter_ = create_unique_contiguous_region_filter(distance=distance,
                                                         genomic_db=genomic_db,
                                            genomic_seqs_fhand=open(genomic_db))
        #an snv in an unique region
        seq  = 'CTGGAATCTCTGAGTTTCTGGGTTCAAGTTGCACTGACCATTGTTGGATTTGTAGATTGTTTC'
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
        assert seq.features[0].qualifiers['filters'][filter_id][distance]

        #an snv in a region with two matches
        seq  = 'CCACTACAAGAGGTGGAAGAGCGAAAACTCTGTTTATTACTAGCTAGGGTTTCTATTAATGAA'
        seq += 'AGGTTCATGTAAATATATGAAGATGGGAAGCAAGAGGTGTTCAAGGAGAAGAGGGAGTTAGAC'
        seq += 'GACCAGAAGATGGCGGTGACGTTTAGAGGGCTGGATGGTCATGTGATGGAGCAGCTTAAGGTG'
        seq += 'TATGACGTCATCTTCCAATTCGTCCCTAAGTCTCAGGAAGGTTGCGTCTGCAAAGTCACTATG'
        seq += 'TTTTGGGAGAAGCGCTA'
        alleles = {('A', SNP): None, ('T', INVARIANT):None}
        snv1 = SeqFeature(type='snv', location=FeatureLocation(100, 100),
                          qualifiers={'alleles':alleles})
        seq = SeqWithQuality(seq=Seq(seq), features=[snv1])
        filter_(seq)
        assert not seq.features[0].qualifiers['filters'][filter_id][distance]

        #a sequence with one hit but two hsps, but a contiguous region according
        #to est2genome
        seq  = 'TAACATATGTATCGTTTGCTAACTGTATATCAGGAGAAGGAGTAATGTAAAAATTCGAGA'
        seq += 'TTATTGTAATTTAAGAAACTCTTTAGGTAGAATGTAATTAATACCAAAATCATTAAAATT'
        seq += 'TCTGTACTTTATACTTGTCATATTCTCAAGATACGAGTCAAAAATGTCTGATTAAATATG'
        seq += 'AACATATTTTTATATATTCATTCCATCATAATCTTCGAGATTTAAAAAGCTCTGATTATC'
        alleles = {('A', SNP): None, ('T', INVARIANT):None}
        snv1 = SeqFeature(type='snv', location=FeatureLocation(50, 50),
                          qualifiers={'alleles':alleles})
        seq = SeqWithQuality(seq=Seq(seq), features=[snv1])
        filter_(seq)
        filter_(seq)
        assert seq.features[0].qualifiers['filters'][filter_id][distance]
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
                                [True, False, True, True]):
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
        seq = SeqWithQuality(seq=seq_str, qual=[30] * len(seq_str),
                             features = [snv1, snv2, snv3])
        max_variability = 40
        filter_ = create_high_variable_region_filter(max_variability)
        filter_(seq)
        threshold = (max_variability, None)
        for snv, expected in zip(seq.get_features(kind='snv'),
                                 [True, True, True]):
            result = snv.qualifiers['filters']['high_variable_reg'][threshold]
            assert result == expected

        max_variability = 60
        filter_ = create_high_variable_region_filter(max_variability)
        filter_(seq)
        filter_(seq)
        threshold = (max_variability, None)
        for snv, expected in zip(seq.get_features(kind='snv'),
                                 [False, False, False]):
            result = snv.qualifiers['filters']['high_variable_reg'][threshold]
            assert result == expected
        max_variability = 25
        window          = 6
        threshold = (max_variability, window)
        filter_   = create_high_variable_region_filter(max_variability,
                                                       window=window)
        filter_(seq)
        for snv, expected in zip(seq.get_features(kind='snv'),
                                 [False, True, True]):
            result = snv.qualifiers['filters']['high_variable_reg'][threshold]
            assert result == expected




if __name__ == "__main__":
    #import sys;sys.argv = ['', 'SeqVariationFilteringTest.test_svn_pipeline']
    unittest.main()
