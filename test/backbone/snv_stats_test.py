'''
Created on 27/06/2011

@author: dani
'''

from Bio.SeqFeature import FeatureLocation

from franklin.utils.misc_utils import NamedTemporaryDir
from franklin.backbone.create_project import create_project
from franklin.backbone.backbone_runner import do_analysis
from franklin.seq.writers import write_seqs_in_file
from franklin.seq.seqs import SeqWithQuality, Seq, SeqFeature
from franklin.snv.snv_annotation import (INVARIANT, SNP)

from os.path import join
import unittest, os

THREADS = False

class SnvStatsTest(unittest.TestCase):
    'it tests the snv stats analysis'
    test_dir = NamedTemporaryDir()
    project_name = 'backbone'

    config = {'General_settings':{'threads':THREADS}}

    settings_path = create_project(directory=test_dir.name,
                                       name=project_name,
                                       configuration=config)
    project_dir = join(test_dir.name, project_name)

    db_dir = join(project_dir, 'annotations', 'db')
    os.makedirs(db_dir)

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
    seq = SeqWithQuality(seq=Seq(seq), name='ref', features=[snv, snv2])

    fhand = open(join(db_dir, 'snv_stats.pickle'), 'a')
    write_seqs_in_file([seq], fhand, format='pickle')

    do_analysis(project_settings=settings_path, kind='snv_stats',
                    silent=True)

    stats_fpath = join(project_dir, 'annotations', 'features', 'stats', 'snv')

    result = open(join(stats_fpath, 'distrib.sum')).read()
    expected = '''Statistics for heterozygosity
------------------------------
minimum: 0
maximum: 0.57
average: 0.28
variance: 0.08
sum: 0.57
items: 2

Statistics for pic
-------------------
minimum: 0
maximum: 0.48
average: 0.24
variance: 0.06
sum: 0.48
items: 2

Statistics for maf
-------------------
minimum: 0.57
maximum: 0.60
average: 0.59
variance: 0.00
sum: 1.17
items: 2'''
    assert expected in result
    test_dir.close()

if    __name__ == "__main__":
    #import sys;sys.argv = ['', '']
    unittest.main()
