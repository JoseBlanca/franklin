'''
Created on 2009 api 28

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
from biolib.contig import Contig
from biolib.locatable_sequence import locate_sequence
from biolib.contig_cleaner import (create_contig_read_stripper,
                                   create_non_matched_region_stripper,
                                   create_read_number_contig_filter)
from biolib.seqs import SeqWithQuality
from Bio.Seq import Seq

class StripTest(unittest.TestCase):
    "It test the contig strip functions, in all methods"

    @staticmethod
    def test_contig_strip():
        '''It tests the contig strip function'''
        contig = Contig()
        seq = 'tTCTTAAGGTTATGCGTACGTGCAGTag'
        mask = (1, 26)
        contig.append_to_location(seq, mask=mask)
        seq = 'tTCTTAAGGTTATGCGTACGTGCAGTag'
        mask = (1, 26)
        contig.append_to_location(seq, mask=mask)
        contig_strip = create_contig_read_stripper(3)
        contig1 = contig_strip(contig)
        assert str(contig1[0]).strip() == 'TAAGGTTATGCGTACGTGCA'
        contig_strip = create_contig_read_stripper(5)
        contig2 = contig_strip(contig)
        assert str(contig2[0]).strip() == 'AGGTTATGCGTACGTG'

        contig3 = Contig()
        seq = Seq('tTCTTAAGGTTAag')
        mask = (2, 12)
        contig3.append_to_location(seq, mask=mask, forward=False, strand=-1)
        contig_strip = create_contig_read_stripper(2)
        contig = contig_strip(contig3)
        assert str(contig[0]).strip() == 'AACCTTA'

    @staticmethod
    def test_contig_water_strip():
        '''It test the contig water strip function '''
        seq1 = SeqWithQuality(name='consensus', seq='AATTCCGG')
        contig = Contig(consensus=locate_sequence(sequence=seq1, location=1))
        seq = 'tAATTCCGGt'
        mask = (1, 8)
        contig.append_to_location(seq, mask=mask)

        water_alignment_strip = create_non_matched_region_stripper()

        contig = water_alignment_strip(contig)
        assert str(contig[0]).strip() == 'AATTCCGG'


        contig = Contig(consensus=locate_sequence(sequence=seq1, location=1))
        seq = 'taaTTCCggt'
        mask = (3, 6)
        contig.append_to_location(seq, mask=mask)

        contig = water_alignment_strip(contig)
        assert str(contig[0]).strip() == 'TTCC'

    @staticmethod
    def test_read_number_contig_filter():
        'It test the read number contig filter'
        seq1   = SeqWithQuality(name='consensus', seq='AATTCCGG')
        contig = Contig(consensus=locate_sequence(sequence=seq1, location=1))
        seq = 'tAATTCCGGt'
        contig.append_to_location(seq)
        contig.append_to_location(seq)
        contig.append_to_location(seq)

        filter_function = create_read_number_contig_filter(2)
        assert filter_function(contig) == True

        filter_function = create_read_number_contig_filter(5)
        assert filter_function(contig) == False

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
