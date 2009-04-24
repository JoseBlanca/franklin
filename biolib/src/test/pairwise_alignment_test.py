'It tests the blast utilities.'

import unittest
from biolib.pairwise_alignment import bl2seq, water

class Bl2seqTest(unittest.TestCase):
    'It tests the bl2seq function'
    @staticmethod
    def test_bl2seq():
        'It tests the bl2seq function'
        seq1 = 'tatcgtactgatgctgatgctatgtgcatgtgcgtatgcgtatgctgagtcgtat' + \
               'tcatgcgtagtctatcgtagtcgtagtgctgtatgcgagctgatgctgatgcgtg'
        seq2 = 'aaacgttggtatgctgatgctatgtgcatgtgcgtatgcgtatgctgagtcgtat' + \
               'tcatgcgtagtctatcgtagtcgtagtgctgtatgcgagctgatgatcgtatgtt'
        result = bl2seq(seq1, seq2)
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
        result = water(seq1, seq2)
        assert isinstance(result['score'], float)
        name1, name2 = sorted(result['alignment'].keys())
        assert result['alignment'][name1]['start'] > \
                                          result['alignment'][name2]['start']
        assert result['alignment'][name1]['end'] > \
                                            result['alignment'][name2]['end']
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
