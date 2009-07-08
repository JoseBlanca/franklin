'''
Created on 2009 uzt 6

@author: peio
'''
import unittest
from biolib.seq_cleaner import (mask_low_complexity,
                                mask_polya,
                                strip_seq_by_quality_trimpoly,
                                strip_seq_by_quality,
                                strip_seq_by_quality_lucy,
                                strip_vector_by_alignment,
                                mask_repeats_by_repeatmasker)
from biolib.seqs import SeqWithQuality
from biolib.biolib_utils import temp_multi_fasta_file

class SeqCleanerTest(unittest.TestCase):
    'It tests cleaner function from seq_cleaner'

    @staticmethod
    def test_strip_seq_by_quality():
        'test trim_seq_by_quality '
        qual = [20, 20, 20, 60, 60, 60, 60, 60, 20, 20, 20, 20]
        seq  = 'ataataataata'
        new_seq = strip_seq_by_quality(SeqWithQuality(qual=qual, seq=seq),
                                       quality_treshold=40, min_seq_length=2,
                                       min_quality_bases=3)
        assert new_seq.seq == 'ataat'
        qual = [60, 60, 60, 60, 60, 60, 60]
        seq  = 'ataataa'
        new_seq = strip_seq_by_quality(SeqWithQuality(qual=qual, seq=seq),
                                       quality_treshold=40, min_seq_length=2,
                                       min_quality_bases=3)
        assert  new_seq.seq == 'ataataa'
        qual = [60, 60, 60, 60, 60, 60, 0]
        seq  = 'ataataa'
        new_seq = strip_seq_by_quality(SeqWithQuality(qual=qual, seq=seq),
                                       quality_treshold=40, min_seq_length=2,
                                       min_quality_bases=3)
        assert new_seq.seq == 'ataata'
        assert new_seq.qual == [60, 60, 60, 60, 60, 60]

    @staticmethod
    def test_mask_low_complexity():
        'It test mask_low_complexity function'
        seq = 'TCGCATCGATCATCGCAGATCGACTGATCGATCGATCGGGGGGGGGGGGGGGGGGGGGGGG'
        qual = [30] * 61
        seq1 = SeqWithQuality(seq=seq, qual=qual)
        masked_seq = mask_low_complexity(seq1)
        assert masked_seq.seq[-10:] == 'gggggggggg'
        assert len(masked_seq.qual) == 61

    @staticmethod
    def test_mask_polya():
        'It test mask_polyA function'
        seq = 'TCGCATCGATCATCGCAGATCGACTGATCGATCGATCAAAAAAAAAAAAAAAAAAAAAAA'
        seq1 = SeqWithQuality(seq=seq)
        masked_seq = mask_polya(seq1)
        exp_seq = 'TCGCATCGATCATCGCAGATCGACTGATCGATCGATCaaaaaaaaaaaaaaaaaaaaaaa'
        assert masked_seq.seq == exp_seq

    def test_trim_seq_by_qual_trimpoly(self):
        'It test trimpoly  but with trim low quality paramters'
        seq  = 'ATCGATCTGATCTAGTCGATGTCTAGCTGAGCTACATAGCTAACGATCTAGTCTAGTCTATG'
        seq += 'TCATGTCATGTCGATGTCTAGTCTAGTCTAGTGAGTCACTGACTAGATCATGACATCGANNN'
        seq += 'NNNNNNNNNNNNNNNNNNTACTAGTC'
        qual = [10] * 150
        seq1 = SeqWithQuality(seq=seq, qual=qual)
        trimmed_seq = strip_seq_by_quality_trimpoly(seq1)
        # It does not mask anything with this example, but it check that
        # if works
        assert trimmed_seq.seq.endswith('ATGACATCGA')

        try:
            seq1 = SeqWithQuality(seq='ATCTATCATATCAT')
            strip_seq_by_quality_trimpoly(seq1)
            self.fail('ValueError expected')
            #pylint: disable-msg=W0704
        except ValueError:
            pass

    def test_strip_seq_by_qual_lucy(self):
        'It tests lucy with a ga runner class for lucy'

        seq =  'ATCGATCAGTCAGACTGACAGACTCAGATCAGATCAGCATCAGCATACGATACGCATCAGACT'
        seq += 'ACGATCGATCGATCGACAGATCATCGATCATCGACGACTAGACGATCATCGATACGCAGACTC'
        seq += 'AGCAGACTACGAGATCAGCAGCATCAGCAGCA'
        qual =  '00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 00 '
        qual += '00 00 00 00 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 '
        qual += '60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 '
        qual += '60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 '
        qual += '60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 '
        qual += '60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 '
        qual += '60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 60 '
        qual += '60 60 60 60 60 60 60 60 60 60 60 60 60 60 00 00 00 00'

        qual = qual.split()
        seqrec = SeqWithQuality(name='hola', seq=seq, qual=qual)
        striped_seq = strip_seq_by_quality_lucy(seqrec)

        seq = striped_seq.seq
        assert seq.startswith('CAGATCAGATCAGCATCAGCAT')
        assert seq.endswith('CGAGATCAGCAGCATCAGC')

        #a sequence too short or with no name will raise an error
        try:
            #short sequence
            seqrec = SeqWithQuality(name='hola', seq='ATCT', qual=[1, 2, 3, 4])
            strip_seq_by_quality_lucy(seqrec)
            self.fail()
            #pylint: disable-msg=W0704
        except ValueError:
            pass
        try:
            #seq with no name
            seqrec = SeqWithQuality(seq=seq, qual=qual)
            strip_seq_by_quality_lucy(seqrec)
            self.fail()
            #pylint: disable-msg=W0704
        except ValueError:
            pass
    @staticmethod
    def test_strip_vector_align_exonera():
        'It tests strip_vector_by_alignment'

        seq  = 'ATGCATCAGATGCATGCATGACTACGACTACGATCAGCATCAGCGATCAGCATCGATACGATC'
        seq2 = SeqWithQuality(name='seq', seq=seq)
        vec1 = SeqWithQuality(name='vec1', seq='atcgatcgatagcatacgat')
        vec2 = SeqWithQuality(name='vec2', seq='atgcatcagatcgataaaga')
        fhand_vectors = temp_multi_fasta_file([vec1, vec2])
        seq1  = SeqWithQuality(name=seq2.name, seq=vec1+seq2+vec2 )
        seq3 = strip_vector_by_alignment(seq1, fhand_vectors, 'exonerate')
        assert str(seq2.seq) == str(seq3.seq)


        fhand_vectors.seek(0)
        seq1  = SeqWithQuality(name=seq2.name, seq=vec1+vec2+seq2 )
        seq3 = strip_vector_by_alignment(seq1, fhand_vectors, 'exonerate')
        assert str(seq2.seq) == str(seq3.seq)

        # overlaping vectors
        fhand_vectors.seek(0)
        seq1  = SeqWithQuality(name=seq2.name, seq=vec1[:-2]+vec2+seq2+vec2)
        seq3 = strip_vector_by_alignment(seq1, fhand_vectors, 'exonerate')
        assert str(seq2.seq) == str(seq3.seq)

        # Now only vectors
        fhand_vectors.seek(0)
        seq1 = SeqWithQuality(name=seq2.name, seq=vec1+vec2+vec2)
        seq3 = strip_vector_by_alignment(seq1, fhand_vectors, 'exonerate')
        assert seq3 is None

        # with some extra seq at the begining and end
        fhand_vectors.seek(0)
        seq1 = SeqWithQuality(name=seq2.name,
                              seq=seq2[:20]+vec1+seq2+vec2+seq2[:20] )
        seq3 = strip_vector_by_alignment(seq1, fhand_vectors, 'exonerate')
        assert str(seq2.seq) == str(seq3.seq)

        # Now without vectors
        fhand_vectors.seek(0)
        seq1 = seq2
        seq3 = strip_vector_by_alignment(seq1, fhand_vectors, 'exonerate')
        assert str(seq2.seq) == str(seq3.seq)

        fhand_vectors.seek(0)
        seq1  = SeqWithQuality(name=seq2.name, seq=vec1[::-1]+vec2+seq2 )
        seq3 = strip_vector_by_alignment(seq1, fhand_vectors, 'exonerate')
        assert str(seq2.seq) == str(seq3.seq)

    @staticmethod
    def test_strip_vector_align_blast():
        'It tests strip_vector_by_alignment using blast and UniVec'
        vector = 'UniVec'
        vec1  = 'CTCGGGCCGTCTCTTGGGCTTGATCGGCCTTCTTGCGCATCTCACGCGCTCCTGCGGCGGCC'
        vec1 += 'TGTAGGGCAGGCTCATACCCCTGCCGAACCGCTTTTGTCAGCCGGTCGGCCACGGCTTCCGG'
        vec1 += 'CGTCTCAACGCGCTTT'
        seq1 = 'ATGCATCAGATGCATGCATGACTACGACTACGATCAGCATCAGCGATCAGCATCGATACGATC'
        seq  = SeqWithQuality(name='seq', seq=seq1+vec1)
        striped_seq = strip_vector_by_alignment(seq, vector, 'blast')
        striped_seq = str(striped_seq.seq)

        assert seq1[4:14]  in striped_seq
        assert seq1[-14:-4]  in striped_seq
        assert vec1[4:14]  not in striped_seq
        assert vec1[-14:-4] not  in striped_seq

    @staticmethod
    def test_repeatmasking():
        'It test that we can mask a repeat element using repeat masker'
        seq  = 'GGTGATGCTGCCAACTTACTGATTTAGTGTATGATGGTGTTTTTGAGGTGCTCCAGTGGCT'
        seq += 'TCTGTTTCTATCAGCTGTCCCTCCTGTTCAGCTACTGACGGGGTGGTGCGTAACGGCAAAA'
        seq += 'GCACCGCCGGACATCAGCGCTATCTCTGCTCTCACTGCCGTAAAACATGGCAACTGCAGTT'
        seq += 'CACTTACACCGCTTCTCAACCCGGTACGCACCAGAAAATCATTGATATGGCCATGAATGGC'
        seq += 'GTTGGATGCCGGGCAACAGCCCGCATTATGGGCGTTGGCCTCAACACGATTTTACGAACCG'
        seq += 'TTTGACTTACGTATTTGCCCATTGTGATTCTAGTCGATTTGCATAACGTGTACGTATCGGT'
        seq += 'ATTGTGACTGATTCGATGCTATTGCAAACGTTTTGATTGTGTGATCGTGATGCATGCTAGT'
        seq += 'CTGATCGAGTCTGATCGTAGTCTAGTCGTAGTCGATGTCGATTTATCGTAGTCGATGCTAG'
        seq += 'TCTAGTCTAGTCTACTAGTCTAGTCATGCTAGTCGAGTCGAT'
        seqrec  = SeqWithQuality(name='seq', seq=seq)
        masked_seq = mask_repeats_by_repeatmasker(seqrec)
        masked_str = str(masked_seq.seq)
        assert seq[0:10].lower() in masked_str
        assert 'tggcctcaacacgat' in masked_str
        assert 'CGTTTGACTT'      in masked_str

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
