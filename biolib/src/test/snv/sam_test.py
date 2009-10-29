'''A samtools pileup parser.'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of biolib.
# biolib is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# biolib is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with biolib. If not, see <http://www.gnu.org/licenses/>.

import unittest, os, StringIO

import biolib
from biolib.seqvar.sam_pileup import (_is_seq_var, _locations_in_pileups,
                                      snvs_in_sam_pileups, _check_pileup)
from biolib.seqvar.seqvariation import INVARIANT, SNP
from biolib.statistics import calculate_read_coverage

DATA_DIR = os.path.join(os.path.split(biolib.__path__[0])[0], 'data')

REFERENCES='''>SGN-U562678
ATATATATATATATATATAT
>SGN-U562679
GCGCGCGCGCGCGG'''

class Test(unittest.TestCase):
    'It tests the samtools pileup parser'

    @staticmethod
    def test_coverage():
        'It test that we can get the coverage from the samtools pileup file'
        sam_fname = os.path.join(DATA_DIR, 'sam.pileup')
        sam_fhand = open(sam_fname)
        coverage = calculate_read_coverage(sam_fhand)
        assert coverage['distrib'] == [16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 20, 0, 0,
                                       0, 0, 0, 0, 0, 0, 30]

    @staticmethod
    def test_is_seq_var():
        'It test is_seq_var function'
        ref_base = 'A'
        quality  = ['~', '~', '~', '~', '~']
        alleles = [{'allele':'A', 'reads':3, 'kind':SNP, 'qualities': quality},
                   {'allele':'T', 'reads':2, 'kind':INVARIANT,
                                                         'qualities': quality}]
        min_num_bases  = 2
        assert _is_seq_var(ref_base, alleles, min_num_bases)

        alleles = [{'allele':'T', 'reads':3, 'kind':SNP, 'qualities': quality},
                   {'allele':'A', 'reads':2, 'kind':INVARIANT,
                                                       'qualities': quality}]
        min_num_bases  = 3
        assert _is_seq_var(ref_base, alleles, min_num_bases)

    @staticmethod
    def test_positions_in_pileups():
        'We can get the equivalent positions from different pileup files'

        pileup1 = '''ref1     1      A      1       ,       ~
ref1     2      A      1       ,       ~
ref1     4      A      1       ,       ~
ref2     2      A      1       ,       ~
ref2     3      A      1       ,       ~'''
        pileup2 = '''ref1     2      A      1       ,       ~
ref1     3      A      1       ,       ~
ref2     1      A      1       ,       ~
ref2     3      A      1       ,       ~
ref2     4      A      1       ,       ~'''
        pileup1 = StringIO.StringIO(pileup1)
        pileup2 = StringIO.StringIO(pileup2)
        locations = list(_locations_in_pileups([pileup1, pileup2]))
        location = locations[0]
        assert location[0][:2] == ['ref1', '1']
        assert not location[1]
        last_loc = locations[7]

        assert not last_loc[0]
        assert last_loc[1][:2] == ['ref2', '4']

    @staticmethod
    def test_seqvars_in_sam_pileups():
        'We can get the Snv from a list of pile ups'
        pileup1 = '''ref1     1      A      2       .t       aa
ref1     2      A      1       ,       ~
ref1     4      A      1       ,       ~
ref2     2      A      1       ,       ~
ref2     3      A      1       ,       a'''
        pileup2 = '''ref1     2      A      1       ,       ~
ref1     3      A      1       ,       ~
ref2     1      A      1       ,       ~
ref2     3      A      1       Tt      AA
ref2     4      A      1       ,       ~'''
        pileup1 = StringIO.StringIO(pileup1)
        pileup2 = StringIO.StringIO(pileup2)

        snvs = list(snvs_in_sam_pileups([pileup1, pileup2],
                                        libraries=['lib1', 'lib2']))
        assert snvs[0].reference == 'ref1'
        assert snvs[0].location == 1
        assert len(snvs[0].per_lib_info[0]['alleles']) == 2 
        assert snvs[0].per_lib_info[0]['alleles'][0]['orientations'] == [True]
        assert snvs[0].per_lib_info[0]['alleles'][1]['orientations'] == [False]

        assert snvs[1].reference == 'ref2'
        assert snvs[1].location == 3
        assert len(snvs[1].per_lib_info[0]['alleles']) == 1
        assert snvs[1].per_lib_info[0]['alleles'][0]['qualities'] == [64]
        assert snvs[1].per_lib_info[1]['alleles'][0]['qualities'] == [32, 32]

    @staticmethod
    def test_pileup_checker():
        'It tests if we can check the pileups'
        pileup1 = '''ref1     1      A      2       ,T       aa
ref1     2      A      1       ,.       ~
ref1     4      A      1       ,.       ~
ref2     2      A      1       ,.       ~
ref2     3      A      1       ,..       ~'''
        pileup2 = '''ref1     2      A      1       ,       ~
ref1     3      A      1       T       ~
ref2     1      A      1       T       ~
ref2     3      A      1       T       ~
ref2     4      A      1       T       ~'''
        pileup3 = '''ref1     2      A      1       ,       ~
ref1     3      N      1       .       ~
ref2     1      N      1       .       ~
ref2     3      N      1       .       ~
ref2     4      N      1       .       ~'''

        pileup1 = StringIO.StringIO(pileup1)
        pileup2 = StringIO.StringIO(pileup2)
        pileup3 = StringIO.StringIO(pileup3)
        assert _check_pileup(pileup1)
        assert not _check_pileup(pileup2)
        assert not _check_pileup(pileup3)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_coverage']
    unittest.main()