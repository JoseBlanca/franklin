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

import unittest, os

import biolib
from biolib.sam import (calculate_read_coverage, seqvars_in_sam_pileup,
                        is_seq_var)
from biolib.seqvariation import INVARIANT, SNP

DATA_DIR = os.path.join(os.path.split(biolib.__path__[0])[0], 'data')

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
        seq_vars = list(seqvars_in_sam_pileup(fhand, 2))
        assert len(seq_vars) == 2
        assert seq_vars[0][0].reference == 'SGN-U562678'
        assert seq_vars[1][0].reference == 'SGN-U562679'

    @staticmethod
    def test_is_seq_bar():
        'It test is_seq_var function'
        coverage = 5
        ref_base = 'A'
        quality  = ['~', '~', '~', '~', '~']
        alleles = [{'allele':'A', 'reads':3, 'kind':SNP, 'quality': quality},
                   {'allele':'T', 'reads':2, 'kind':INVARIANT,
                                                         'quality': quality}]
        min_num_bases  = 2
        assert is_seq_var(coverage, ref_base, alleles, min_num_bases)

        read_bases = 'TTT..'
        alleles = [{'allele':'T', 'reads':3, 'kind':SNP, 'quality': quality},
                   {'allele':'A', 'reads':2, 'kind':INVARIANT,
                                                       'quality': quality}]
        min_num_bases  = 3
        assert is_seq_var(coverage, ref_base, alleles, min_num_bases)





if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_coverage']
    unittest.main()