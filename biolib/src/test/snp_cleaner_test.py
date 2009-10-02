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

from biolib.seqvar.seqvariation import (SeqVariation, SNP, DELETION, INVARIANT)
from biolib.seqvar.snp_cleaner import (create_major_allele_freq_filter,
                                       create_pic_filter,
                                       create_close_to_seqvar_filter,
                                       create_cap_enzyme_filter,
                                       create_high_variable_region_filter,
                                       create_allele_number_filter,
                                       create_bad_quality_reads_cleaner,
                                       create_seqvar_close_to_limit_filter)
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

    @staticmethod
    def test_major_allele_freq_filter():
        'It test the first allele percent filter'
        snp = SeqVariation(alleles=[{'allele':'A', 'reads':4, 'kind':INVARIANT},
                                     {'allele':'T', 'reads':2,'kind':SNP}],
                           location=11, reference='reference')
        snp = (snp, 'context')
        first_percent = create_major_allele_freq_filter(0.6)
        assert first_percent(snp) == False
        first_percent = create_major_allele_freq_filter(0.8)
        assert first_percent(snp) == True

    @staticmethod
    def test_create_filter_close_seqvar():
        'It test the first allele percent filter'
        snp = SeqVariation(alleles=[{'allele':'A', 'reads':4, 'kind':INVARIANT},
                                     {'allele':'T', 'reads':2,'kind':SNP}],
                           location=10, reference='reference')
        snp1 = SeqVariation(alleles=[{'allele':'A', 'reads':4,
                                      'kind':INVARIANT},
                                     {'allele':'T', 'reads':2,'kind':SNP}],
                           location=20, reference='reference')

        snp2 = SeqVariation(alleles=[{'allele':'A', 'reads':4,
                                      'kind':INVARIANT},
                                     {'allele':'T', 'reads':2,'kind':SNP}],
                           location=15, reference='reference')
        context = [snp, snp1, snp2]
        snp = (snp, context)

        filter4 = create_close_to_seqvar_filter(5)
        filter7 = create_close_to_seqvar_filter(7)
        assert filter4(snp)
        assert not filter7(snp)

    @staticmethod
    def test_high_variable_region_filter():
        'It test percent_variations_in_seq_ref filter'
        reference = 'atatat'
        snp = SeqVariation(alleles=[{'allele':'A', 'reads':4, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}],
                           location=1, reference=reference)
        snp1 = SeqVariation(alleles=[{'allele':'A', 'reads':4,
                                      'kind':INVARIANT},
                                     {'allele':'T', 'reads':2,'kind':SNP}],
                           location=2, reference=reference)

        snp2 = SeqVariation(alleles=[{'allele':'A', 'reads':4,
                                      'kind':INVARIANT},
                                     {'allele':'T', 'reads':2,'kind':SNP}],
                           location=4, reference=reference)
        context = [snp, snp1, snp2]
        snp     = (snp1, context)
        # This snp variability is 0.75 for window_witdt = 2
        filter40 = create_high_variable_region_filter(0.4, 2)
        assert filter40(snp)
        filter60 = create_high_variable_region_filter(0.8, 2)
        assert not filter60(snp)

    @staticmethod
    def test_allele_number_filter():
        'It test percent_variations_in_seq_ref filter'
        reference = 'atatat'
        snp = SeqVariation(alleles=[{'allele':'A', 'reads':4, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}],
                           location=1, reference=reference)
        snp1 = SeqVariation(alleles=[{'allele':'A', 'reads':4,
                                      'kind':INVARIANT}],
                           location=4, reference=reference)
        filter_ = create_allele_number_filter(2)
        snp = (snp, 'context')
        snp1 = (snp1, 'context')
        assert filter_(snp)
        assert not filter_(snp1)

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

        snp = SeqVariation(alleles=[{'allele':'A', 'reads':4, 'kind':INVARIANT,
                                     'quality':[30, 30, 30, 20]},
                                    {'allele':'T', 'reads':2,'kind':SNP,
                                      'quality':[30, 20]}],
                           location=1, reference='reference')
        snp = (snp, 'context')
        snp = bad_quality_cleaner(snp)[0]
        assert len(snp.alleles)    == 2
        assert snp.alleles[0]['reads'] == 3

        snp = SeqVariation(alleles=[{'allele':'A', 'reads':4, 'kind':INVARIANT,
                                     'quality':[30, 30, 30, 20]},
                                    {'allele':'T', 'reads':2,'kind':SNP,
                                      'quality':[20, 20]}],
                           location=1, reference='reference')
        snp = (snp, 'context')
        snp = bad_quality_cleaner(snp)[0]
        assert len(snp.alleles)    == 1
        assert snp.alleles[0]['reads'] == 3

    @staticmethod
    def test_create_seqvar_close_to_limt():
        'It tests the close to limit filter function '
        snp = SeqVariation(alleles=[{'allele':'A', 'reads':4, 'kind':INVARIANT},
                                    {'allele':'T', 'reads':2,'kind':SNP}],
                           location=5, reference='atatatata')
        snp = (snp, 'context')
        filter_ = create_seqvar_close_to_limit_filter(3)
        assert filter_(snp)
        filter_ = create_seqvar_close_to_limit_filter(5)
        assert not filter_(snp)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_SeqVariation_init']
    unittest.main()
