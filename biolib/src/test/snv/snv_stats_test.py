'''
Created on 08/10/2009

@author: jose
'''

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

import unittest
from StringIO import StringIO

from biolib.snv.snv import Snv, svn_contexts_in_file, snvs_in_file, SNP
from biolib.snv.snv_stats import (calculate_ref_variability_ditrib,
                                  snv_distrib)

class TestSnvStats(unittest.TestCase):
    'It tests the snv statistics'
    @staticmethod
    def test_ref_variability_distrib():
        'It tests the reference variability distribution'
        seq1 = 'A' * 100
        seq2 = 'T' * 100
        reference_fhand = StringIO('>hola\n' + seq1 + '\n>caracola\n' + seq2 +
                                   '\n')
        seq_var = Snv(reference='hola', location=3, per_lib_info=[])
        seq_var2 = Snv(reference='hola', location=4, per_lib_info=[])
        seq_var3 = Snv(reference='caracola', location=2, per_lib_info=[])
        snv_fhand = StringIO(repr(seq_var) + repr(seq_var2) + repr(seq_var3))
        snv_contexts = svn_contexts_in_file(snv_fhand, reference_fhand)
        distrib = calculate_ref_variability_ditrib(snv_contexts)

        assert distrib['distrib'] == [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 0, 0, 2]
    @staticmethod
    def test_maf_distrib():
        'It tests the reference variability distribution'
        libs1 = [{'alleles':[{'allele':'A', 'reads':3, 'kind':SNP},
                             {'allele':'T', 'reads':3, 'kind':SNP}]}]
        libs2 = [{'alleles':[{'allele':'A', 'reads':30, 'kind':SNP},
                             {'allele':'T', 'reads':1, 'kind':SNP}]}]
        libs3 = [{'alleles':[{'allele':'A', 'reads':9, 'kind':SNP},
                             {'allele':'T', 'reads':1, 'kind':SNP}]}]
        seq_var = Snv(reference='hola', location=3, per_lib_info=libs1)
        seq_var2 = Snv(reference='hola', location=4, per_lib_info=libs2)
        seq_var3 = Snv(reference='caracola', location=2, per_lib_info=libs3)
        snv_fhand = StringIO(repr(seq_var) + repr(seq_var2) + repr(seq_var3))
        snvs = snvs_in_file(snv_fhand)

        distrib = snv_distrib(snvs, kind='maf_distrib')
        assert distrib['distrib'] == [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                      0, 0, 0, 1, 0, 1]

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_ref_variability_distrib']
    unittest.main()