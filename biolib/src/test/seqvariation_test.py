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
from biolib.seqvar.seqvariation import (calculate_pic, cap_enzime,
                                 SeqVariation, SNP, INSERTION, DELETION,
                                 INVARIANT, INDEL, COMPLEX, Snv, snvs_in_file )
from biolib.seqs import SeqWithQuality

class SeqVariationTest(unittest.TestCase):
    '''Here we will check if the SeqVariation module works as it should.'''
    @staticmethod
    def test_init():
        '''It tests the init'''
        #required fields
        alleles = [{'allele':'A', 'reads':3, 'kind':SNP},
                   {'allele':'T', 'reads':4, 'kind':INVARIANT }]

        svar = SeqVariation(alleles=alleles, reference='ref')
        assert len(svar.alleles) == 2

    @staticmethod
    def test_kind():
        'It test that we can get the kind for a seqvariation'

        seq_var = SeqVariation(alleles= [{'allele':'A', 'reads':3,
                                          'kind':INVARIANT}],
                                reference='hola')
        assert seq_var.kind == INVARIANT

        seq_var = SeqVariation(alleles= [{'allele':'A', 'reads':3,
                                           'kind':DELETION}],
                                reference='hola')
        assert seq_var.kind == DELETION

        seq_var = SeqVariation(alleles= [{'allele':'A', 'reads':3,
                                           'kind':INSERTION}],
                                reference='hola')
        assert seq_var.kind == INSERTION

        seq_var = SeqVariation(alleles= [{'allele':'A', 'reads':3,
                                           'kind':INSERTION},
                                          {'allele':'T',  'reads':2,
                                          'kind':DELETION}],
                                reference='hola')
        assert seq_var.kind == INDEL

        seq_var = SeqVariation(alleles= [{'allele':'A', 'reads':3,
                                           'kind':INSERTION},
                                          {'allele':'T', 'reads':2,
                                           'kind':SNP}],
                                reference='hola')
        assert seq_var.kind == COMPLEX

        seq_var = SeqVariation(alleles= [{'allele':'A', 'reads':3,
                                           'kind':INVARIANT},
                                          {'allele':'T', 'reads':2,
                                           'kind':SNP}],
                                reference='hola')
        assert seq_var.kind == SNP

    @staticmethod
    def test_sorted_alleles():
        'It checks that we can get the alleles sorted by the number of reads.'
        snp = SeqVariation(alleles=[{'allele':'A', 'reads':2},
                                     {'allele':'T', 'reads':3}],
                            reference='hola')
        alleles = snp.alleles
        assert alleles[0]['allele'] == 'T'
        assert alleles[1]['allele'] == 'A'

class SnvTest(unittest.TestCase):
    '''Here we will check if the Snv module works as it should.'''
    @staticmethod
    def test_init():
        '''It tests the init'''
        lib_alleles = [{'alleles':[{'allele':'A', 'reads':3, 'kind':SNP},
                                  {'allele':'T', 'reads':4, 'kind':INVARIANT}]}]

        snv = Snv(lib_alleles=lib_alleles, reference='ref', location=2)
        assert len(snv.lib_alleles) == 1

    @staticmethod
    def test_kind():
        'It test that we can get the kind for a seqvariation'

        seq_var = Snv(reference='hola', location=3,
                      lib_alleles=[{'alleles':[{'allele':'A', 'reads':3,
                                               'kind':INVARIANT}]}])

        assert seq_var.kind == INVARIANT

        seq_var = Snv(reference='hola', location=3,
                      lib_alleles=[{'alleles':[{'allele':'A', 'reads':3,
                                               'kind':DELETION}]}])
        assert seq_var.kind == DELETION
        seq_var = Snv(reference='hola', location=3,
                      lib_alleles=[{'alleles':[{'allele':'A', 'reads':3,
                                               'kind':DELETION},
                                               {'allele':'A', 'reads':3,
                                               'kind':INSERTION}]},
                                   {'alleles':[{'allele':'A', 'reads':3,
                                               'kind':INSERTION}]},
                                   {'alleles':[{'allele':'A', 'reads':3,
                                               'kind':INVARIANT}]}])
        assert seq_var.kind == INDEL

        seq_var = Snv(reference='hola', location=3,
                      lib_alleles=[{'alleles':[{'allele':'A', 'reads':3,
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
                      lib_alleles=[{'alleles':[{'allele':'A', 'reads':3,
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
        assert seq_var.lib_alleles == snv.lib_alleles







class SeqVariationCaracterization(unittest.TestCase):
    '''It tests seqvar caracterization functions  '''
    @staticmethod
    def test_calculate_pic():
        'It checks that we are able to calculate the PIC values'
        snp1 = SeqVariation(alleles=[{'allele':'A', 'reads':200, 'kind':SNP},
                                      {'allele':'T', 'reads':200, 'kind':SNP}],
                             reference='hola')
        snp2 = SeqVariation(alleles=[{'allele':'A', 'reads':300, 'kind':SNP},
                                      {'allele':'T', 'reads':100, 'kind':SNP}],
                             reference='hola')
        pic_1 = calculate_pic(snp1)
        pic_2 = calculate_pic(snp2)
        assert pic_1 > pic_2

class SeqVariationrEnzime(unittest.TestCase):
    ''' It checks if we have problems with remaps and it functions'''

    @staticmethod
    def test_remap():
        '''It test if the remap external program works '''
        reference = SeqWithQuality(seq='Actgacttactgtca', name='ref')
        snp = SeqVariation(alleles=[{'allele':'C', 'reads':2, 'kind':SNP},
                                     {'allele':'T', 'reads':3,
                                      'kind':INVARIANT}],
                            location=7, reference=reference)
        enzymes = cap_enzime(snp, True)
        assert ['HinfI', 'TscAI'] == enzymes

        # With a deletion
        seq  = 'ATGATGATG' + 'gaaattc' + 'ATGATGATGTGGGAT'
        reference = SeqWithQuality(seq=seq, name='ref')
        snp = SeqVariation(alleles=[{'allele':'A', 'reads':2, 'kind':DELETION},
                                     {'allele':'A', 'reads':3,
                                      'kind':INVARIANT}],
                            location=11, reference=reference)
        enzymes = cap_enzime(snp, True)
        assert 'EcoRI' in enzymes

        #with an insertion
        seq  = 'ATGATGATG' + 'gaattc' + 'ATGATGATGTGGGAT'
        reference = SeqWithQuality(seq=seq, name='ref')
        snp = SeqVariation(alleles=[{'allele':'A', 'reads':2, 'kind':INSERTION},
                                     {'allele':'A', 'reads':3,
                                      'kind':INVARIANT}],
                            location=11, reference=reference)
        enzymes = cap_enzime(snp, True)
        assert 'EcoRI' in enzymes

    @staticmethod
    def test_svns_in_file():
        'It we can read the svn file'
        seq_var = Snv(reference='hola', location=3,
                      lib_alleles=[{'alleles':[{'allele':'A', 'reads':3,
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
        assert seq_var.lib_alleles == snv.lib_alleles

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_SeqVariation_init']
    unittest.main()
