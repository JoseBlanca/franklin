'''
Created on 2009 uzt 30

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

from biolib.seqvar.seqvariation import (SeqVariation, SNP, DELETION, INVARIANT,
                                        Snv)
from biolib.seqvar.snp_cleaner import (#create_major_allele_freq_filter,
                                       create_pic_filter,
                                       create_close_to_seqvar_filter,
                                       create_cap_enzyme_filter,
                                       create_high_variable_region_filter,
                                       create_allele_number_cleaner,
                                       create_bad_quality_reads_cleaner,
                                       create_snv_close_to_limit_filter,
                                       create_major_allele_freq_cleaner,
                                       create_read_number_cleaner)
from biolib.seqs import SeqWithQuality

class SeqVariationFilteringTest(unittest.TestCase):
    'It checks the filtering methods.'

    @staticmethod
    def test_pic_filter():
        '''It test if we can filter by pic '''
        alleles = [{'allele':'A', 'reads':200, 'kind':SNP},
                   {'allele':'T', 'reads':200, 'kind':INVARIANT }]
        snp1 = SeqVariation(alleles=alleles, reference='ref')

        alleles = [{'allele':'A', 'reads':300, 'kind':SNP},
                   {'allele':'T', 'reads':100, 'kind':INVARIANT }]
        snp2 = SeqVariation(alleles=alleles, reference='ref')
        pic_filter = create_pic_filter(0.35)
        snp1 = (snp1, 'context')
        assert pic_filter(snp1) == True
        snp2 = (snp2, 'context')
        assert pic_filter(snp2) == False

    @staticmethod
    def test_enzyme_filter1():
        'It test the enzyme filter'
        seq  = 'ATGATGATG' + 'gaaattc' + 'ATGATGATGTGGGAT'
        reference = SeqWithQuality(seq=seq, name='ref')
        snp = SeqVariation(alleles=[{'allele':'A', 'reads':2, 'kind':DELETION},
                                     {'allele':'A', 'reads':3,
                                      'kind':INVARIANT}],
                            location=11, reference=reference)
        cap_enzime = create_cap_enzyme_filter(True)
        snp = (snp, 'context')
        assert cap_enzime(snp) == True

##############################################################################
    @staticmethod
    def test_major_allele_freq_filter_snv():
        'It test the first allele percent filter'
        snp = Snv(lib_alleles=[
                        {'alleles':[{'allele':'A', 'reads':4, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]},
                        {'alleles':[{'allele':'A', 'reads':2, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]}],
                  location=11, reference='reference')
        snp = (snp, 'context')
        frec_cleaner = create_major_allele_freq_cleaner(0.51)
        snp = frec_cleaner(snp)[0]

        assert len(snp.lib_alleles) == 1

        snp = Snv(location=11, reference='reference', lib_alleles=[
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]},
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]}])
        snp = (snp, 'context')
        frec_cleaner = create_major_allele_freq_cleaner(0.51)
        snp = frec_cleaner(snp)
        assert  snp is None

    @staticmethod
    def test_high_variable_region_filter():
        'It test percent_variations_in_seq_ref filter'
        reference = 'atatat'
        snp = Snv(location=1, reference=reference, lib_alleles=[
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]},
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]}])

        snp1 = Snv(location=4, reference=reference, lib_alleles=[
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]},
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]}])

        snp2 = Snv(location=7, reference=reference, lib_alleles=[
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]},
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]}])

        context = [snp, snp1, snp2]
        snp     = (snp1, context)
        # This reference variability is 0.5
        filter40 = create_high_variable_region_filter(0.4)

        filter60 = create_high_variable_region_filter(0.8)
        assert filter60(snp)

    @staticmethod
    def test_create_filter_close_seqvar():
        'It test the first allele percent filter'
        reference = 'reference'
        snp = Snv(location=10, reference=reference, lib_alleles=[
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]},
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]}])

        snp1 = Snv(location=15, reference=reference, lib_alleles=[
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]},
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]}])

        snp2 = Snv(location=20, reference=reference, lib_alleles=[
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]},
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]}])
        context = [snp, snp1, snp2]
        snp = (snp, context)

        filter4 = create_close_to_seqvar_filter(5)
        filter7 = create_close_to_seqvar_filter(7)
        assert filter4(snp)
        assert not filter7(snp)

    @staticmethod
    def test_allele_number_cleaner():
        'It test percent_variations_in_seq_ref filter'
        reference = 'atatat'
        snp = Snv(location=10, reference=reference, lib_alleles=[
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]},
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]}])

        snp1 = Snv(location=15, reference=reference, lib_alleles=[
                    {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT}]},
                    {'alleles':[{'allele':'T', 'reads':2,'kind':SNP}]}])

        snp2 = Snv(location=20, reference=reference, lib_alleles=[
                      {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT}]},
                      {'alleles':[{'allele':'T', 'reads':2,'kind':INVARIANT}]}])

        snp  = (snp, [snp])
        snp1 = (snp1, [snp1])
        snp2 = (snp2, [snp2])
        allele_number_cleaner = create_allele_number_cleaner(2)
        snp = allele_number_cleaner(snp)
        snp1 = allele_number_cleaner(snp1)
        snp2 = allele_number_cleaner(snp2)

        assert len(snp[0].lib_alleles) == 2
        assert len(snp[0].lib_alleles[0]['alleles']) == 2
        assert len(snp[0].lib_alleles[1]['alleles']) == 2

        assert len(snp1[0].lib_alleles) == 1
        assert len(snp1[0].lib_alleles[0]['alleles']) == 1

        assert snp2 is None

    @staticmethod
    def test_create_seqvar_close_to_limt():
        'It tests the close to limit filter function '
        snp = Snv(location=4, reference='atatatatat', lib_alleles=[
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]},
                        {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}]}])
        snp = (snp, 'context')
        filter_ = create_snv_close_to_limit_filter(3)
        assert filter_(snp)
        filter_ = create_snv_close_to_limit_filter(5)
        assert not filter_(snp)

    @staticmethod
    def test_remove_bad_quality_reads():
        'It test it we can remove the alleles with bad qualities'
        #the good reads in 454 are around a quality of 25 to 40
        #a bad letter is less than 20
        #for sanger the phred quality values are related to the probability
        #of having an error in the sequence
        #http://www.phrap.com/phred/
        #Phred quality score  Probability base is wrong   Accuracy of the base
        #       10                   1 in 10                  90%
        #       20                   1 in 100                 99%
        #       30                   1 in 1,000               99.9%
        #       40                   1 in 10,000              99.99%
        #       50                   1 in 100,000             99.999%
        bad_quality_cleaner = create_bad_quality_reads_cleaner(25)
        snp = Snv(location=4, reference='atatatatat', lib_alleles=[
                        {'alleles':[{'allele':'A', 'reads':4, 'kind':INVARIANT,
                                     'quality':[30, 30, 30, 20]},
                                    {'allele':'T', 'reads':2,'kind':SNP,
                                     'quality':[30, 20]}]}])
        snp = (snp, 'context')
        snp = bad_quality_cleaner(snp)[0]
        assert len(snp.lib_alleles[0]['alleles']) == 2
        assert snp.lib_alleles[0]['alleles'][0]['reads'] == 3

        snp = Snv(location=4, reference='atatatatat', lib_alleles=[
                {'alleles':[{'allele':'A', 'reads':3, 'kind':INVARIANT,
                             'quality':[30, 30, 30, 20]},
                            {'allele':'T', 'reads':2,'kind':SNP,
                             'quality':[20, 20]}]}])
        snp = (snp, 'context')
        snp = bad_quality_cleaner(snp)[0]

        assert len(snp.lib_alleles[0]['alleles']) == 1
        assert snp.lib_alleles[0]['alleles'][0]['reads'] == 3

        snp = Snv(location=4, reference='atatatatat', lib_alleles=[
                {'alleles':[{'allele':'AT', 'reads':4, 'kind':DELETION,
                             'quality':[None, None, None, None]},
                            {'allele':'T', 'reads':2,'kind':SNP,
                             'quality':[30, 30]}]}])
        snp = (snp, 'context')
        snp = bad_quality_cleaner(snp)[0]

        assert len(snp.lib_alleles[0]['alleles']) == 2
        assert snp.lib_alleles[0]['alleles'][0]['reads'] == 4
        assert len(snp.lib_alleles[0]['alleles'][0]['quality']) == 4
    @staticmethod
    def test_read_number_cleaner():
        'Tests read number_cleaner'
        snp = Snv(location=4, reference='atatatatat', lib_alleles=[
                        {'alleles':[{'allele':'A', 'reads':4, 'kind':INVARIANT,
                                     'quality':[30, 30, 30, 20]},
                                    {'allele':'T', 'reads':1,'kind':SNP,
                                     'quality':[30, 20]}]}])
        snp = (snp, 'context')
        read_number_cleaner = create_read_number_cleaner(2)
        snp = read_number_cleaner(snp)[0]
        assert len(snp.lib_alleles[0]['alleles']) == 1

        snp = Snv(location=4, reference='atatatatat', lib_alleles=[
                        {'alleles':[{'allele':'A', 'reads':1, 'kind':INVARIANT,
                                     'quality':[30, 30, 30, 20]},
                                    {'allele':'T', 'reads':1,'kind':SNP,
                                     'quality':[30, 20]}]}])
        snp = (snp, 'context')
        read_number_cleaner = create_read_number_cleaner(2)
        snp = read_number_cleaner(snp)
        assert snp is None

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_SeqVariation_init']
    unittest.main()
