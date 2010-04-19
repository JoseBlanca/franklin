'It tests the blast utilities.'

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of franklin.
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
from franklin.pairwise_alignment import bl2seq, water
from franklin.seq.seqs import SeqWithQuality, Seq

class Bl2seqTest(unittest.TestCase):
    'It tests the bl2seq function'
    @staticmethod
    def test_bl2seq():
        'It tests the bl2seq function'
        seq1 = 'tatcgtactgatgctgatgctatgtgcatgtgcgtatgcgtatgctgagtcgtat' + \
               'tcatgcgtagtctatcgtagtcgtagtgctgtatgcgagctgatgctgatgcgtg'
        seq2 = 'aaacgttggtatgctgatgctatgtgcatgtgcgtatgcgtatgctgagtcgtat' + \
               'tcatgcgtagtctatcgtagtcgtagtgctgtatgcgagctgatgatcgtatgtt'
        result = bl2seq(SeqWithQuality(Seq(seq1), name='seq1'),
                       SeqWithQuality(Seq(seq2), name='seq2'))
        #we check only one hsp
        result = result[0]
        assert isinstance(result['evalue'], float)
        name1, name2 = sorted(result['alignment'].keys())
        assert result['alignment'][name1]['start'] < 20
        assert result['alignment'][name2]['start'] < 20
        assert result['alignment'][name1]['end'] > 80
        assert result['alignment'][name2]['end'] > 20

    @staticmethod
    def test_water():
        'It tests the water function'
        seq1 = 'atatcgtactgatgctgatgctatgtgcatgtgcgtatgcgtatgctgagtcgtat' + \
               'tcatgcgtagtctatcgtagtcgtagtgctgtatgcgagctgatgctgatgcgtg'
        seq2 = 'aacgttggtatgctgatgctatgtgcatgtgcgtatgcgtatgctgagtcgtat' + \
               'tcatgcgtagtctatcgtagtcgtagtgctgtatgcgagctgatgatcgtatgtt'
        result = water(SeqWithQuality(Seq(seq1), name='seq1'),
                       SeqWithQuality(Seq(seq2), name='seq2'))
        assert isinstance(result['score'], float)
        name1, name2 = sorted(result['alignment'].keys())
        assert result['alignment'][name1]['start'] > \
                                          result['alignment'][name2]['start']
        assert result['alignment'][name1]['end'] > \
                                            result['alignment'][name2]['end']
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
