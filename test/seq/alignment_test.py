'''
Created on 06/04/2011

@author: jose
'''
import unittest

from franklin.seq.alignment import BlastAligner, sw_align, match_words
from franklin.seq.seqs import SeqWithQuality, Seq

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

    def test_biopython_alignment(self):
        'We can align with biopython local and global'
        query   = SeqWithQuality(Seq('TACTGGCTTT'))
        subject = SeqWithQuality(Seq('TACTGcCTTTCC'))
        sw_align(query, subject)

        query   = SeqWithQuality(Seq('TACTGGCTTT'))
        subject = SeqWithQuality(Seq('CCTACTGGCTTT'))
        sw_align(query, subject)

        query   = SeqWithQuality(Seq('CCTACTGGCTTT'))
        subject = SeqWithQuality(Seq('TACTGGCTTT'))
        sw_align(query, subject)

        query   = SeqWithQuality(Seq('TACTGGCTTTCC'))
        subject = SeqWithQuality(Seq('TACTGGCTTT'))
        sw_align(query, subject)

        query   = SeqWithQuality(Seq('ATGTGGCTTT'))
        subject = SeqWithQuality(Seq('TACTGGCTTT'))
        sw_align(query, subject)

        query   = SeqWithQuality(Seq('ATGTGGCTTT'))
        subject = SeqWithQuality(Seq('TCCTGAGT'))
        sw_align(query, subject)

class WordMatchTest(unittest.TestCase):
    'It test that we can match words against sequences'

    @staticmethod
    def test_forward_words():
        'It test that we can match words against in the same orientation'

        seq = 'gCACAggTGTGggTATAgg'
        seq = SeqWithQuality(seq=Seq(seq))

        result = match_words(seq, ['CACA', 'TATA', 'KK'])[0]
        assert result['query'] == seq

        #The match por CACA
        match = result['matches'][0]
        assert match['subject'] == 'CACA'
        assert match['start'] == 1
        assert match['end'] == 10
        assert len(match['match_parts']) == 2
        #the reverse match part
        assert match['match_parts'][1] == {'query_start':7,
                                           'query_end':10,
                                           'query_strand':1,
                                           'subject_start':0,
                                           'subject_end':3,
                                           'subject_strand':-1}

        #The match por TATA
        match = result['matches'][1]
        assert match['subject'] == 'TATA'
        assert match['start'] == 13
        assert match['end'] == 16
        assert len(match['match_parts']) == 2

        #No matches for KK
        assert len(result['matches']) == 2


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_blast_alignment']
    unittest.main()