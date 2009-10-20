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
from biolib.seqvar.sam_pileup import (seqvars_in_sam_pileup, _is_seq_var,
                                      locations_in_pileups)
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
    def test_seqvars_in_sam_pileup():
        'It tests the parser of the sam pileup file'
        sam_fname = os.path.join(DATA_DIR, 'sam.pileup')
        fhand = open(sam_fname)
        references = StringIO.StringIO(REFERENCES)
        snvs = list(seqvars_in_sam_pileup(fhand, references=references))
        assert len(snvs) == 9
        assert snvs[0][0].reference.name == 'SGN-U562678'
        assert snvs[7][0].reference.name == 'SGN-U562679'
        assert snvs[0][0].reference.seq == 'ATATATATATATATATATAT'
        assert snvs[7][0].reference.seq == 'GCGCGCGCGCGCGG'
        assert snvs[0][0].per_lib_info[0]['alleles'][0]['allele'] == 'G'
        assert snvs[0][0].per_lib_info[0]['alleles'][0]['qualities'] == [None,
                                                                     None, 20]
        assert snvs[2][0].per_lib_info[0]['alleles'][0]['reads'] == 3
        assert snvs[1][0].per_lib_info[0]['alleles'][0]['kind'] == INVARIANT

    @staticmethod
    def test_is_seq_bar():
        'It test is_seq_var function'
        coverage = 5
        ref_base = 'A'
        quality  = ['~', '~', '~', '~', '~']
        alleles = [{'allele':'A', 'reads':3, 'kind':SNP, 'qualities': quality},
                   {'allele':'T', 'reads':2, 'kind':INVARIANT,
                                                         'qualities': quality}]
        min_num_bases  = 2
        assert _is_seq_var(coverage, ref_base, alleles, min_num_bases)

        alleles = [{'allele':'T', 'reads':3, 'kind':SNP, 'qualities': quality},
                   {'allele':'A', 'reads':2, 'kind':INVARIANT,
                                                       'qualities': quality}]
        min_num_bases  = 3
        assert _is_seq_var(coverage, ref_base, alleles, min_num_bases)

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
        locations = list(locations_in_pileups([pileup1, pileup2]))
        location = locations[0]
            loc1, loc2 = location[0:2]
            loc1 = loc1[:2]
            loc2 = loc2[:2]
            print loc1, loc2


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_coverage']
    unittest.main()