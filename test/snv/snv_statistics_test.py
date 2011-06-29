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
from franklin.seq.readers import seqs_in_file
from franklin.seq.writers import write_seqs_in_file
from franklin.seq.seqs import SeqWithQuality, Seq, SeqFeature
from franklin.snv.snv_annotation import (INVARIANT, SNP)
from os.path import join


class SnvAnalysisTest(unittest.TestCase):
    'It checks the analysis of snvs'

    @staticmethod
    def test_create_pic_distribution():

        test_dir = NamedTemporaryDir()

        alleles = {('A', SNP): {'read_groups':{'rg1':1, 'rg2':1, 'rg4':2}},
               ('T', INVARIANT): {'read_groups':{'rg1':1, 'rg3':2}}}
        snv = SeqFeature(type='snv', location=FeatureLocation(11, 11),
                         qualifiers={'alleles':alleles,
                                     'read_groups':{},
                                     'heterozygosity': 0.57,
                                     'pic': 0.48})
        seq = 'ATGATGATGgaaattcATGATGATGTGGGAT'

        alleles2 = {('A', SNP): {'read_groups':{'rg1':2}},
                    ('T', INVARIANT): {'read_groups':{'rg1':1, 'rg3':2}}}
        snv2 = SeqFeature(type='snv', location=FeatureLocation(21, 21),
                         qualifiers={'alleles':alleles2,
                                     'read_groups':{},
                                     'heterozygosity': 0,
                                     'pic': 0})
        seq1 = SeqWithQuality(seq=Seq(seq), name='seq1', features=[snv, snv2])

        alleles = {('A', SNP): {'read_groups':{'rg1':1, 'rg2':1}},
               ('T', INVARIANT): {'read_groups':{'rg1':1, 'rg2':2}}}
        snv = SeqFeature(type='snv', location=FeatureLocation(11, 11),
                         qualifiers={'alleles':alleles,
                                     'read_groups':{},
                                     'heterozygosity': 0.80,
                                     'pic': 0.75})
        seq = 'ATGATGATGgaaattcATGATGATGTGGGAT'

        seq2 = SeqWithQuality(seq=Seq(seq), name='seq2', features=[snv])

        fhand = open(join(test_dir.name, 'seqs.pickle'), 'a')
        write_seqs_in_file([seq1, seq2], fhand, format='pickle')
        seqs = seqs_in_file(open(join(test_dir.name, 'seqs.pickle'), 'r'))

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

    @staticmethod
    def test_create_het_distribution():

        test_dir = NamedTemporaryDir()

        alleles = {('A', SNP): {'read_groups':{'rg1':1, 'rg2':1, 'rg4':2}},
               ('T', INVARIANT): {'read_groups':{'rg1':1, 'rg3':2}}}
        snv = SeqFeature(type='snv', location=FeatureLocation(11, 11),
                         qualifiers={'alleles':alleles,
                                     'read_groups':{},
                                     'heterozygosity': 0.57,
                                     'pic': 0.48})
        seq = 'ATGATGATGgaaattcATGATGATGTGGGAT'

        alleles2 = {('A', SNP): {'read_groups':{'rg1':2}},
                    ('T', INVARIANT): {'read_groups':{'rg1':1, 'rg3':2}}}
        snv2 = SeqFeature(type='snv', location=FeatureLocation(21, 21),
                         qualifiers={'alleles':alleles2,
                                     'read_groups':{},
                                     'heterozygosity': 0,
                                     'pic': 0})
        seq1 = SeqWithQuality(seq=Seq(seq), name='seq1', features=[snv, snv2])

        alleles = {('A', SNP): {'read_groups':{'rg1':1, 'rg2':1}},
               ('T', INVARIANT): {'read_groups':{'rg1':1, 'rg2':2}}}
        snv = SeqFeature(type='snv', location=FeatureLocation(11, 11),
                         qualifiers={'alleles':alleles,
                                     'read_groups':{},
                                     'heterozygosity': 0.80,
                                     'pic': 0.75})
        seq = 'ATGATGATGgaaattcATGATGATGTGGGAT'

        seq2 = SeqWithQuality(seq=Seq(seq), name='seq2', features=[snv])

        fhand = open(join(test_dir.name, 'seqs.pickle'), 'a')
        write_seqs_in_file([seq1, seq2], fhand, format='pickle')
        seqs = seqs_in_file(open(join(test_dir.name, 'seqs.pickle'), 'r'))

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

    @staticmethod
    def test_create_maf_distribution():

        test_dir = NamedTemporaryDir()

        alleles = {('A', SNP): {'read_groups':{'rg1':1, 'rg2':1, 'rg4':2}},
               ('T', INVARIANT): {'read_groups':{'rg1':1, 'rg3':2}}}
        snv = SeqFeature(type='snv', location=FeatureLocation(11, 11),
                         qualifiers={'alleles':alleles,
                                     'read_groups':{},
                                     'heterozygosity': 0.57,
                                     'pic': 0.48})
        seq = 'ATGATGATGgaaattcATGATGATGTGGGAT'

        alleles2 = {('A', SNP): {'read_groups':{'rg1':2}},
                    ('T', INVARIANT): {'read_groups':{'rg1':1, 'rg3':2}}}
        snv2 = SeqFeature(type='snv', location=FeatureLocation(21, 21),
                         qualifiers={'alleles':alleles2,
                                     'read_groups':{},
                                     'heterozygosity': 0,
                                     'pic': 0})
        seq1 = SeqWithQuality(seq=Seq(seq), name='seq1', features=[snv, snv2])

        alleles = {('A', SNP): {'read_groups':{'rg1':1, 'rg2':1}},
               ('T', INVARIANT): {'read_groups':{'rg1':1, 'rg2':2}}}
        snv = SeqFeature(type='snv', location=FeatureLocation(11, 11),
                         qualifiers={'alleles':alleles,
                                     'read_groups':{},
                                     'heterozygosity': 0.80,
                                     'pic': 0.75})
        seq = 'ATGATGATGgaaattcATGATGATGTGGGAT'

        seq2 = SeqWithQuality(seq=Seq(seq), name='seq2', features=[snv])

        fhand = open(join(test_dir.name, 'seqs.pickle'), 'a')
        write_seqs_in_file([seq1, seq2], fhand, format='pickle')
        seqs = seqs_in_file(open(join(test_dir.name, 'seqs.pickle'), 'r'))

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

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'SeqVariationFilteringTest.test_create_pic_distribution']
    unittest.main()
