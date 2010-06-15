'''
Created on 14/06/2010

@author: jose
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
from StringIO import StringIO

from franklin.seq.seqs import SeqWithQuality, Seq, SeqFeature, FeatureLocation
from franklin.seq.seq_stats import (do_annotation_statistics,
                                    _location_to_orf)


class SeqStatsTest(unittest.TestCase):
    'It tests the statistics for the sequences'

    @staticmethod
    def test_annotation_statistics():
        'It tests the seq annotation statistics'

        orf = SeqFeature(FeatureLocation(1, 3), type='orf',
                         qualifiers={'strand':'reverse'})
        intron = SeqFeature(FeatureLocation(2, 2), type='intron')
        ssr1 = SeqFeature(FeatureLocation(2, 2), type='microsatellite',
                          qualifiers={'unit':'GAA'})
        ssr2 = SeqFeature(FeatureLocation(7, 10), type='microsatellite',
                          qualifiers={'unit':'GAA'})
        snv = SeqFeature(FeatureLocation(7, 7), type='snv',
                          qualifiers={'alleles':{('T',0):'', ('G',0):''}})

        feats1 = [intron, ssr1]
        feats2 = [orf, ssr2, snv]

        annots1 = {'arabidopsis-orthologs':['arab_1']}

        annots2 = {'melon-orthologs':['melon_1'],
                   'GOs':['Go1', 'Go2']}

        seq1 = SeqWithQuality(Seq('ACTG'), description='hola', features=feats1,
                              annotations=annots1)
        seq2 = SeqWithQuality(Seq('ACTG'), features=feats2, annotations=annots2)
        seqs = [seq1, seq2]

        stats_fhand = StringIO()

        do_annotation_statistics(seqs, stats_fhand)
        result = stats_fhand.getvalue()

        expected = '''Annotation statistics
---------------------
Number of sequences: 2
Sequences with description: 1
Sequences with ORF: 1
Number of ORFs: 1
Sequences with intron: 1
Number of introns: 1

Orthologs
_________
Sequences with melon orthologs: 1
Number of melon orthologs: 1
Sequences with arabidopsis orthologs: 1
Number of arabidopsis orthologs: 1

GO terms
________
Sequences with GOs: 1
Number of GOs: 2

SNVs
____
Sequences with SNVs: 2
SNV types:
\ttransversion: 1
SNV locations:
\tunknown: 1

Microsatellites
_______________
Sequences with microsatellites: 2
Microsatellite types:
\ttrinucleotide: 2
Microsatellite locations:
\tin 5 prime: 1
\tunknown: 1

'''
        result = result.splitlines()
        expected = expected.splitlines()
        for index, line in enumerate(result):
            #print line, expected[index]
            assert line == expected[index]

    @staticmethod
    def test_location_to_orf():
        'It test the location respect the orf'
        orf = SeqFeature(FeatureLocation(10, 20), type='orf',
                         qualifiers={'strand':'reverse'})

        feat = SeqFeature(FeatureLocation(5, 5), type='snv')

        assert _location_to_orf([orf], feat) == 'in 3 prime'

        orf = SeqFeature(FeatureLocation(10, 20), type='orf',
                         qualifiers={'strand':'forward'})
        assert _location_to_orf([orf], feat) == 'in 5 prime'

        feat = SeqFeature(FeatureLocation(15, 15), type='snv')
        assert _location_to_orf([orf], feat) == 'in orf'

        feat = SeqFeature(FeatureLocation(25, 25), type='snv')
        assert _location_to_orf([orf], feat) == 'in 3 prime'

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()