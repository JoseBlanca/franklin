'''
Created on 2009 mar 25

@author: peio
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
from biolib.seqvar.seqvariation import (pic, cap_enzymes, SNP,
                                        INSERTION, DELETION, INVARIANT, INDEL,
                                        COMPLEX, Snv, snvs_in_file,
                                        major_allele_frequency,
                                        reference_variability,
                                        calculate_kind,
                                        mayor_frec_allele_per_library)
from biolib.seqs import SeqWithQuality

class SnvTest(unittest.TestCase):
    '''Here we will check if the Snv module works as it should.'''
    @staticmethod
    def test_init():
        '''It tests the init'''
        per_lib_info = [{'alleles':[{'allele':'A', 'reads':3, 'kind':SNP},
                                  {'allele':'T', 'reads':4, 'kind':INVARIANT}]}]

        snv = Snv(per_lib_info=per_lib_info, reference='ref', location=2)
        assert len(snv.per_lib_info) == 1

    @staticmethod
    def test_kind():
        'It test that we can get the kind for a seqvariation'

        seq_var = Snv(reference='hola', location=3,
                      per_lib_info=[{'alleles':[{'allele':'A', 'reads':3,
                                               'kind':INVARIANT}]}])

        assert seq_var.kind == INVARIANT

        seq_var = Snv(reference='hola', location=3,
                      per_lib_info=[{'alleles':[{'allele':'A', 'reads':3,
                                               'kind':DELETION}]}])
        assert seq_var.kind == DELETION
        seq_var = Snv(reference='hola', location=3,
                      per_lib_info=[{'alleles':[{'allele':'A', 'reads':3,
                                               'kind':DELETION},
                                               {'allele':'A', 'reads':3,
                                               'kind':INSERTION}]},
                                   {'alleles':[{'allele':'A', 'reads':3,
                                               'kind':INSERTION}]},
                                   {'alleles':[{'allele':'A', 'reads':3,
                                               'kind':INVARIANT}]}])
        assert seq_var.kind == INDEL

        seq_var = Snv(reference='hola', location=3,
                      per_lib_info=[{'alleles':[{'allele':'A', 'reads':3,
                                               'kind':SNP},
                                               {'allele':'A', 'reads':3,
                                               'kind':INSERTION}]},
                                   {'alleles':[{'allele':'A', 'reads':3,
                                               'kind':DELETION}]},
                                   {'alleles':[{'allele':'A', 'reads':3,
                                               'kind':INVARIANT}]}])
        assert seq_var.kind == COMPLEX

    @staticmethod
    def test_snv_repr():
        'It'
        seq_var = Snv(reference='hola', location=3,
                      per_lib_info=[{'alleles':[{'allele':'A', 'reads':3,
                                               'kind':SNP},
                                               {'allele':'A', 'reads':3,
                                               'kind':INSERTION}]},
                                   {'library':'library',
                                    'alleles':[{'allele':'A', 'reads':3,
                                               'kind':DELETION}]},
                                   {'alleles':[{'allele':'A', 'reads':3,
                                               'kind':INVARIANT}]}])
        snv = eval(repr(seq_var))
        assert seq_var.reference == snv.reference
        assert seq_var.per_lib_info == snv.per_lib_info

    @staticmethod
    def test_sorted_alleles():
        'It checks that we can get the alleles sorted by the number of reads.'
        per_lib_info = [{'alleles':[{'allele':'A', 'reads':2},
                                   {'allele':'T', 'reads':3}]}]
        snp = Snv(per_lib_info=per_lib_info, reference='hola', location=3)
        alleles = snp.per_lib_info[0]['alleles']
        assert alleles[0]['allele'] == 'T'
        assert alleles[1]['allele'] == 'A'

    @staticmethod
    def test_libraries_per_allele():
        'It checks that we can use libraries_per_allele method'
        per_lib_info = [{'library':'lib1',
                         'alleles':[{'allele':'A', 'reads':2},
                                   {'allele':'T', 'reads':3}]},
                        {'library':'lib2',
                         'alleles':[{'allele':'A', 'reads':2}]},
                        {'library':'lib3',
                         'alleles':[{'allele':'A', 'reads':2},
                                   {'allele':'T', 'reads':3}]}]
        snp = Snv(per_lib_info=per_lib_info, reference='hola', location=3)
        lib_per_allele = snp.libraries_per_allele()
        assert lib_per_allele[0]['allele'] == 'A'
        assert lib_per_allele[1]['allele'] == 'T'
        assert len(lib_per_allele[0]['libraries']) == 3
        assert len(lib_per_allele[1]['libraries']) == 2

    @staticmethod
    def test_frec_allele_per_library():
        'It tests frec_allele_per_library'
        per_lib_info = [{'library':'lib1',
                         'alleles':[{'allele':'A', 'reads':2},
                                   {'allele':'T', 'reads':3}]},
                        {'library':'lib2',
                         'alleles':[{'allele':'A', 'reads':2}]},
                        {'library':'lib3',
                         'alleles':[{'allele':'A', 'reads':2},
                                   {'allele':'T', 'reads':3}]}]
        snp = Snv(per_lib_info=per_lib_info, reference='hola', location=3)
        assert mayor_frec_allele_per_library(snp) == 0.5







class SnvCaracterizationTest(unittest.TestCase):
    'It checks that the svns are properly analyzed'
    @staticmethod
    def test_maf():
        'It checks that we can calculate the maf frequency'
        per_lib_info = [{'alleles':[{'allele':'A', 'reads':2},
                                   {'allele':'T', 'reads':2}]}]
        snp = Snv(per_lib_info=per_lib_info, reference='hola', location=3)
        mafs = major_allele_frequency(snp)
        assert mafs[0] == 0.5

    @staticmethod
    def test_variability():
        'It checks that we can calculate the reference variability'
        reference = SeqWithQuality(seq='Actgacttac', name='ref')
        snp = Snv(per_lib_info=[], reference=reference, location=3)
        ref_var = reference_variability(snp, ['snp1'])
        assert ref_var == 10.0

    @staticmethod
    def test_calculate_pic():
        'It checks that we are able to calculate the PIC values'
        per_lib_info = [{'alleles':[{'allele':'A', 'reads':200},
                                   {'allele':'T', 'reads':200}]}]
        snp1 = Snv(per_lib_info=per_lib_info, reference='hola', location=3)
        per_lib_info = [{'alleles':[{'allele':'A', 'reads':300},
                                   {'allele':'T', 'reads':100}]}]
        snp2 = Snv(per_lib_info=per_lib_info, reference='hola', location=3)

        pic_1 = pic(snp1)[0]
        pic_2 = pic(snp2)[0]
        assert pic_1 > pic_2
        pic_1_1 = pic(snp1)[0]
        #we check the buffer
        assert pic_1 == pic_1_1

    @staticmethod
    def test_remap():
        '''It test if the remap external program works '''
        reference = SeqWithQuality(seq='Actgacttactgtca', name='ref')
        per_lib_info = [{'alleles':[{'allele':'C', 'reads':1, 'kind':SNP},
                                  {'allele':'T', 'reads':1, 'kind':INVARIANT}]}]
        snp = Snv(per_lib_info=per_lib_info, reference=reference, location=7)
        enzymes = cap_enzymes(snp, True)
        assert set(['HinfI', 'TscAI']) == set(enzymes)

        # With a deletion
        seq  = 'ATGATGATG' + 'gaaattc' + 'ATGATGATGTGGGAT'
        reference = SeqWithQuality(seq=seq, name='ref')
        per_lib_info = [{'alleles':[{'allele':'A', 'reads':1, 'kind':DELETION},
                                  {'allele':'A', 'reads':1, 'kind':INVARIANT}]}]
        snp = Snv(per_lib_info=per_lib_info, reference=reference, location=11)
        enzymes = cap_enzymes(snp, True)
        assert 'EcoRI' in enzymes

        #with an insertion
        seq  = 'ATGATGATG' + 'gaattc' + 'ATGATGATGTGGGAT'
        reference = SeqWithQuality(seq=seq, name='ref')
        per_lib_info = [{'alleles':[{'allele':'A', 'reads':1, 'kind':INSERTION},
                                  {'allele':'A', 'reads':1, 'kind':INVARIANT}]}]
        snp = Snv(per_lib_info=per_lib_info, reference=reference, location=11)
        enzymes = cap_enzymes(snp, True)
        assert 'EcoRI' in enzymes

        #with only one allele
        reference = SeqWithQuality(seq='Actgacttactgtca', name='ref')
        per_lib_info = [{'alleles':[{'allele':'C', 'reads':1, 'kind':SNP}]}]
        snp = Snv(per_lib_info=per_lib_info, reference=reference, location=7)
        enzymes = cap_enzymes(snp, True)
        assert [] == enzymes
    @staticmethod
    def test_calculate_kind():
        'It calculates the rerulting kind'
        assert calculate_kind(SNP, INVARIANT) == SNP
        assert calculate_kind(INSERTION, DELETION) == INDEL
        assert calculate_kind(COMPLEX, SNP) == COMPLEX
        assert calculate_kind(SNP, INDEL) == COMPLEX


        kind1 = SNP
        kind2 = SNP

class SnvIOTest(unittest.TestCase):
    ''' It checks if we have problems with remaps and it functions'''

    @staticmethod
    def test_svns_in_file():
        'It we can read the svn file'
        seq_var = Snv(reference='hola', location=3,
                      per_lib_info=[{'alleles':[{'allele':'A', 'reads':3,
                                               'kind':SNP},
                                               {'allele':'A', 'reads':3,
                                               'kind':INSERTION}]},
                                   {'library':'library',
                                    'alleles':[{'allele':'A', 'reads':3,
                                               'kind':DELETION}]},
                                   {'alleles':[{'allele':'A', 'reads':3,
                                               'kind':INVARIANT}]}])
        fhand = StringIO(repr(seq_var))
        snv   = snvs_in_file(fhand).next()
        assert seq_var.reference == snv.reference
        assert seq_var.per_lib_info == snv.per_lib_info

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_SeqVariation_init']
    unittest.main()
