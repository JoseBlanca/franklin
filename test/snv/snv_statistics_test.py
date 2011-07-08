'''
Created on 29/06/2011

@author: dani
'''

import unittest

from Bio.SeqFeature import FeatureLocation

from franklin.utils.misc_utils import NamedTemporaryDir
from franklin.snv.snv_statistics import (create_pic_distribution,
                                         create_het_distribution,
                                         create_maf_distribution)
from franklin.seq.seqs import SeqWithQuality, Seq, SeqFeature
from franklin.snv.snv_annotation import (INVARIANT, SNP)
from os.path import join


class SnvAnalysisTest(unittest.TestCase):
    'It checks the analysis of snvs'

    @staticmethod
    def _create_test_seqs():
        'create test seqs with snvs'

        alleles = {('A', SNP): {'read_groups':{'rg1':1, 'rg2':1, 'rg4':2}},
               ('T', INVARIANT): {'read_groups':{'rg1':1, 'rg3':2}}}
        read_groups = {'rg1':{'LB':'l1'}, 'rg2':{'LB':'l2'}, 'rg3':{'LB':'l3'},
                       'rg4':{'LB':'l4'}}
        snv = SeqFeature(type='snv', location=FeatureLocation(11, 11),
                         qualifiers={'alleles':alleles,
                                     'read_groups':read_groups,
                                     'heterozygosity': 0.57,
                                     'pic': 0.48})
        seq = 'ATGATGATGgaaattcATGATGATGTGGGAT'

        alleles2 = {('A', SNP): {'read_groups':{'rg1':2}},
                    ('T', INVARIANT): {'read_groups':{'rg1':1, 'rg3':2}}}
        read_groups2 = {'rg1':{'LB':'l1'}, 'rg3':{'LB':'l3'}}
        snv2 = SeqFeature(type='snv', location=FeatureLocation(21, 21),
                         qualifiers={'alleles':alleles2,
                                     'read_groups':read_groups2,
                                     'heterozygosity': 0,
                                     'pic': 0})
        seq1 = SeqWithQuality(seq=Seq(seq), name='seq1', features=[snv, snv2])

        alleles = {('A', SNP): {'read_groups':{'rg1':1, 'rg2':1}},
               ('T', INVARIANT): {'read_groups':{'rg1':1, 'rg2':2}}}
        read_groups3 = {'rg1':{'LB':'l1'}, 'rg2':{'LB':'l2'}}
        snv3 = SeqFeature(type='snv', location=FeatureLocation(11, 11),
                         qualifiers={'alleles':alleles,
                                     'read_groups':read_groups3,
                                     'heterozygosity': 0.80,
                                     'pic': 0.75})
        seq = 'ATGATGATGgaaattcATGATGATGTGGGAT'

        seq2 = SeqWithQuality(seq=Seq(seq), name='seq2', features=[snv3])
        seqs = [seq1, seq2]
        return seqs

    @staticmethod
    def _create_test_seqs_with_more_reads():
        'create test seqs with snvs'

        alleles = {('A', SNP): {'read_groups':{'rg1':5, 'rg2':1, 'rg4':2}},
               ('T', INVARIANT): {'read_groups':{'rg1':1, 'rg3':2}}}
        read_groups = {'rg1':{'LB':'l1'}, 'rg2':{'LB':'l2'}, 'rg3':{'LB':'l3'},
                       'rg4':{'LB':'l4'}}
        snv = SeqFeature(type='snv', location=FeatureLocation(11, 11),
                         qualifiers={'alleles':alleles,
                                     'read_groups':read_groups,
                                     'heterozygosity': 0.57,
                                     'pic': 0.48})
        seq = 'ATGATGATGgaaattcATGATGATGTGGGAT'

        alleles2 = {('A', SNP): {'read_groups':{'rg1':2}},
                    ('T', INVARIANT): {'read_groups':{'rg1':6, 'rg3':2}}}
        read_groups2 = {'rg1':{'LB':'l1'}, 'rg3':{'LB':'l3'}}
        snv2 = SeqFeature(type='snv', location=FeatureLocation(21, 21),
                         qualifiers={'alleles':alleles2,
                                     'read_groups':read_groups2,
                                     'heterozygosity': 0,
                                     'pic': 0})
        seq1 = SeqWithQuality(seq=Seq(seq), name='seq1', features=[snv, snv2])

        alleles = {('A', SNP): {'read_groups':{'rg1':5, 'rg2':7}},
               ('T', INVARIANT): {'read_groups':{'rg1':5, 'rg2':7}}}
        read_groups3 = {'rg1':{'LB':'l1'}, 'rg2':{'LB':'l2'}}
        snv3 = SeqFeature(type='snv', location=FeatureLocation(11, 11),
                         qualifiers={'alleles':alleles,
                                     'read_groups':read_groups3,
                                     'heterozygosity': 0.80,
                                     'pic': 0.75})
        seq = 'ATGATGATGgaaattcATGATGATGTGGGAT'

        seq2 = SeqWithQuality(seq=Seq(seq), name='seq2', features=[snv3])
        seqs = [seq1, seq2]
        return seqs

    def test_create_pic_distribution(self):
        'It tests the calculation of the pic stats'
        test_dir = NamedTemporaryDir()
        seqs = self._create_test_seqs_with_more_reads()

        distrib_fhand = open(join(test_dir.name, 'distrib_fhand'), 'w')
        plot_fhand = open(join(test_dir.name, 'plot_fhand'), 'w')
        summary_fhand = open(join(test_dir.name, 'summary_fhand'), 'w')

        create_pic_distribution(seqs, distrib_fhand=distrib_fhand,
                                plot_fhand=plot_fhand,
                                summary_fhand=summary_fhand)
        distrib_fhand.close()
        plot_fhand.close()
        summary_fhand.close()

        result_fpath = join(test_dir.name, 'distrib_fhand')
        results = open(result_fpath).read()
        expected = '''values: 3
bin_edges: 0,1'''
        assert expected in results

        result_fpath = join(test_dir.name, 'summary_fhand')
        results = open(result_fpath).read()
        expected = '''Statistics for pic
-------------------
minimum: 0
maximum: 0.7500
average: 0.4100
variance: 0.0962
sum: 1.2300
items: 3'''
        assert expected in results

        distrib_fhand = open(join(test_dir.name, 'distrib_fhand'), 'w')
        plot_fhand = open(join(test_dir.name, 'plot_fhand'), 'w')
        summary_fhand = open(join(test_dir.name, 'summary_fhand'), 'w')

        create_pic_distribution(seqs, distrib_fhand=distrib_fhand,
                                plot_fhand=plot_fhand,
                                summary_fhand=summary_fhand,
                                groups=['l1'],
                                group_kind='libraries')
        distrib_fhand.close()
        plot_fhand.close()
        summary_fhand.close()

        result_fpath = join(test_dir.name, 'summary_fhand')
        results = open(result_fpath).read()
        expected = '''Statistics for pic (libraries: l1)
-----------------------------------
minimum: 0.3333
maximum: 0.4762
average: 0.4008
variance: 0.0034
sum: 1.2024
items: 3'''
        assert expected in results

        distrib_fhand = open(join(test_dir.name, 'distrib_fhand'), 'w')
        plot_fhand = open(join(test_dir.name, 'plot_fhand'), 'w')
        summary_fhand = open(join(test_dir.name, 'summary_fhand'), 'w')

        create_pic_distribution(seqs, distrib_fhand=distrib_fhand,
                                plot_fhand=plot_fhand,
                                summary_fhand=summary_fhand,
                                groups=['l2'],
                                group_kind='libraries')
        distrib_fhand.close()
        plot_fhand.close()
        summary_fhand.close()

        result_fpath = join(test_dir.name, 'summary_fhand')
        results = open(result_fpath).read()
        expected = '''Statistics for pic (libraries: l2)
-----------------------------------
minimum: 0.4650
maximum: 0.4650
average: 0.4650
variance: 0.0000
sum: 0.4650
items: 1'''
        assert expected in results

        distrib_fhand = open(join(test_dir.name, 'distrib_fhand'), 'w')
        plot_fhand = open(join(test_dir.name, 'plot_fhand'), 'w')
        summary_fhand = open(join(test_dir.name, 'summary_fhand'), 'w')

        create_pic_distribution(seqs, distrib_fhand=distrib_fhand,
                                plot_fhand=plot_fhand,
                                summary_fhand=summary_fhand,
                                groups=['l4'],
                                group_kind='libraries')
        distrib_fhand.close()
        plot_fhand.close()
        summary_fhand.close()

        result_fpath = join(test_dir.name, 'summary_fhand')
        results = open(result_fpath).read()
        expected = ''''''
        assert expected == results


    def test_create_het_distribution(self):
        'It tests the calculation of the heterozygosity stats'
        test_dir = NamedTemporaryDir()
        seqs = self._create_test_seqs()

        distrib_fhand = open(join(test_dir.name, 'distrib_fhand'), 'w')
        plot_fhand = open(join(test_dir.name, 'plot_fhand'), 'w')
        summary_fhand = open(join(test_dir.name, 'summary_fhand'), 'w')

        create_het_distribution(seqs, distrib_fhand=distrib_fhand,
                                plot_fhand=plot_fhand,
                                summary_fhand=summary_fhand)
        distrib_fhand.close()
        plot_fhand.close()
        summary_fhand.close()

        result_fpath = join(test_dir.name, 'distrib_fhand')
        results = open(result_fpath).read()
        expected = '''values: 3
bin_edges: 0,1'''
        assert expected in results

        result_fpath = join(test_dir.name, 'summary_fhand')
        results = open(result_fpath).read()
        expected = '''Statistics for heterozygosity
------------------------------
minimum: 0
maximum: 0.8000
average: 0.4567
variance: 0.1131
sum: 1.3700
items: 3'''
        assert expected in results

        seqs = self._create_test_seqs()
        distrib_fhand = open(join(test_dir.name, 'distrib_fhand'), 'w')
        plot_fhand = open(join(test_dir.name, 'plot_fhand'), 'w')
        summary_fhand = open(join(test_dir.name, 'summary_fhand'), 'w')

        create_het_distribution(seqs, distrib_fhand=distrib_fhand,
                                plot_fhand=plot_fhand,
                                summary_fhand=summary_fhand,
                                groups=['l1'],
                                group_kind='libraries')
        distrib_fhand.close()
        plot_fhand.close()
        summary_fhand.close()

        result_fpath = join(test_dir.name, 'summary_fhand')
        results = open(result_fpath).read()
        expected = '''Statistics for heterozygosity (libraries: l1)
----------------------------------------------
minimum: 0.5333
maximum: 0.6667
average: 0.6222
variance: 0.0040
sum: 1.8667
items: 3'''
        assert expected in results

        seqs = self._create_test_seqs()
        distrib_fhand = open(join(test_dir.name, 'distrib_fhand'), 'w')
        plot_fhand = open(join(test_dir.name, 'plot_fhand'), 'w')
        summary_fhand = open(join(test_dir.name, 'summary_fhand'), 'w')

        create_het_distribution(seqs, distrib_fhand=distrib_fhand,
                                plot_fhand=plot_fhand,
                                summary_fhand=summary_fhand,
                                groups=['l2'],
                                group_kind='libraries')
        distrib_fhand.close()
        plot_fhand.close()
        summary_fhand.close()

        result_fpath = join(test_dir.name, 'summary_fhand')
        results = open(result_fpath).read()
        expected = '''Statistics for heterozygosity (libraries: l2)
----------------------------------------------
minimum: 0.0000
maximum: 0.5333
average: 0.2667
variance: 0.0711
sum: 0.5333
items: 2'''
        assert expected in results

        seqs = self._create_test_seqs()
        distrib_fhand = open(join(test_dir.name, 'distrib_fhand'), 'w')
        plot_fhand = open(join(test_dir.name, 'plot_fhand'), 'w')
        summary_fhand = open(join(test_dir.name, 'summary_fhand'), 'w')

        create_het_distribution(seqs, distrib_fhand=distrib_fhand,
                                plot_fhand=plot_fhand,
                                summary_fhand=summary_fhand,
                                groups=['l4'],
                                group_kind='libraries')
        distrib_fhand.close()
        plot_fhand.close()
        summary_fhand.close()

        result_fpath = join(test_dir.name, 'summary_fhand')
        results = open(result_fpath).read()
        expected = '''Statistics for heterozygosity (libraries: l4)
----------------------------------------------
minimum: 0.0000
maximum: 0.0000
average: 0.0000
variance: 0.0000
sum: 0.0000
items: 1'''
        assert expected in results

    def test_create_maf_distribution(self):
        'It tests the calculation of the maf stats'
        test_dir = NamedTemporaryDir()
        seqs = self._create_test_seqs()

        distrib_fhand = open(join(test_dir.name, 'distrib_fhand'), 'w')
        plot_fhand = open(join(test_dir.name, 'plot_fhand'), 'w')
        summary_fhand = open(join(test_dir.name, 'summary_fhand'), 'w')

        create_maf_distribution(seqs, distrib_fhand=distrib_fhand,
                                plot_fhand=plot_fhand,
                                summary_fhand=summary_fhand)
        distrib_fhand.close()
        plot_fhand.close()
        summary_fhand.close()

        result_fpath = join(test_dir.name, 'summary_fhand')
        results = open(result_fpath).read()
        expected = '''Statistics for maf
-------------------
minimum: 0.5714
maximum: 0.6000
average: 0.5905
variance: 0.0002
sum: 1.7714
items: 3'''
        assert expected in results

        seqs = self._create_test_seqs()
        distrib_fhand = open(join(test_dir.name, 'distrib_fhand'), 'w')
        plot_fhand = open(join(test_dir.name, 'plot_fhand'), 'w')
        summary_fhand = open(join(test_dir.name, 'summary_fhand'), 'w')

        create_maf_distribution(seqs, distrib_fhand=distrib_fhand,
                                plot_fhand=plot_fhand,
                                summary_fhand=summary_fhand,
                                groups=['l1'],
                                group_kind='libraries')
        distrib_fhand.close()
        plot_fhand.close()
        summary_fhand.close()

        result_fpath = join(test_dir.name, 'summary_fhand')
        results = open(result_fpath).read()
        expected = '''Statistics for maf (libraries: l1)
-----------------------------------
minimum: 0.5000
maximum: 0.6667
average: 0.5556
variance: 0.0062
sum: 1.6667
items: 3'''
        assert expected in results

        seqs = self._create_test_seqs()
        distrib_fhand = open(join(test_dir.name, 'distrib_fhand'), 'w')
        plot_fhand = open(join(test_dir.name, 'plot_fhand'), 'w')
        summary_fhand = open(join(test_dir.name, 'summary_fhand'), 'w')

        create_maf_distribution(seqs, distrib_fhand=distrib_fhand,
                                plot_fhand=plot_fhand,
                                summary_fhand=summary_fhand,
                                groups=['l2'],
                                group_kind='libraries')
        distrib_fhand.close()
        plot_fhand.close()
        summary_fhand.close()

        result_fpath = join(test_dir.name, 'summary_fhand')
        results = open(result_fpath).read()
        expected ='''Statistics for maf (libraries: l2)
-----------------------------------
minimum: 0.6667
maximum: 1.0000
average: 0.8333
variance: 0.0278
sum: 1.6667
items: 2'''
        assert expected in results

        seqs = self._create_test_seqs()
        distrib_fhand = open(join(test_dir.name, 'distrib_fhand'), 'w')
        plot_fhand = open(join(test_dir.name, 'plot_fhand'), 'w')
        summary_fhand = open(join(test_dir.name, 'summary_fhand'), 'w')

        create_maf_distribution(seqs, distrib_fhand=distrib_fhand,
                                plot_fhand=plot_fhand,
                                summary_fhand=summary_fhand,
                                groups=['l4'],
                                group_kind='libraries')
        distrib_fhand.close()
        plot_fhand.close()
        summary_fhand.close()

        result_fpath = join(test_dir.name, 'summary_fhand')
        results = open(result_fpath).read()
        expected ='''Statistics for maf (libraries: l4)
-----------------------------------
minimum: 1.0000
maximum: 1.0000
average: 1.0000
variance: 0.0000
sum: 1.0000
items: 1'''
        assert expected in results

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'SeqVariationFilteringTest.test_create_pic_distribution']
    unittest.main()
