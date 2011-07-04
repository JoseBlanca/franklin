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
maximum: 0.75
average: 0.41
variance: 0.10
sum: 1.23
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
minimum: 0.33
maximum: 0.48
average: 0.40
variance: 0.00
sum: 1.20
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
minimum: 0.47
maximum: 0.47
average: 0.47
variance: 0.00
sum: 0.47
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
maximum: 0.80
average: 0.46
variance: 0.11
sum: 1.37
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
minimum: 0.53
maximum: 0.67
average: 0.62
variance: 0.00
sum: 1.87
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
minimum: 0.00
maximum: 0.53
average: 0.27
variance: 0.07
sum: 0.53
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
minimum: 0.00
maximum: 0.00
average: 0.00
variance: 0.00
sum: 0.00
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
minimum: 0.57
maximum: 0.60
average: 0.59
variance: 0.00
sum: 1.77
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
minimum: 0.50
maximum: 0.67
average: 0.56
variance: 0.01
sum: 1.67
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
minimum: 0.67
maximum: 1.00
average: 0.83
variance: 0.03
sum: 1.67
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
minimum: 1.00
maximum: 1.00
average: 1.00
variance: 0.00
sum: 1.00
items: 1'''
        assert expected in results

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'SeqVariationFilteringTest.test_create_pic_distribution']
    unittest.main()
