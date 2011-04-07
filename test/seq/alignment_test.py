'''
Created on 06/04/2011

@author: jose
'''
import unittest

from franklin.seq.alignment import BlastAligner

class PairwiseAlignmentTest(unittest.TestCase):
    'It tests the different classes of Pairwise alignments'

    def test_blast_alignment(self):
        'We can align using blast'
        seq1 = 'ACTACGGTTACACACGTGTATCAGTTACACAGTGTTGTCATCACTATCTAGTCAGTAGTCTAG'
        seq2_1 = 'CACGCTAGTCGTAGTCGCTAGT'
        seq2_2 = 'CACGCTAGTCGTAGTCGCTCGT'

        blast_params = {'task': 'blastn-short', 'expect': 0.001}
        aligner = BlastAligner(subject=seq2_1, parameters=blast_params)
        alignments = list(aligner.do_alignment(seq1 + seq2_1))
        match = alignments[0]['matches'][0]
        assert match['start'] == 63
        assert match['end'] == 84

        filters = [{'kind':'score_threshold', 'score_key': 'expect',
                    'max_score': 1e-20}]
        aligner = BlastAligner(subject=seq2_1, parameters=blast_params,
                                       filters=filters)
        alignments = list(aligner.do_alignment(seq1 + seq2_1))
        assert not alignments


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_blast_alignment']
    unittest.main()