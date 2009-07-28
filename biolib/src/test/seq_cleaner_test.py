'''
Created on 2009 uzt 6

@author: peio
'''
import unittest, os
import biolib
from tempfile import NamedTemporaryFile

from biolib.seqs import SeqWithQuality
from biolib.biolib_utils import NamedTemporaryDir
from biolib.biolib_seqio_utils import temp_multi_fasta_file
from biolib.seq_cleaner import (create_vector_striper_by_alignment,
                                create_masker_for_polia,
                                create_masker_for_low_complexity,
                                create_striper_by_quality_trimpoly,
                                create_striper_by_quality,
                                create_striper_by_quality_lucy,
                                create_striper_by_quality_lucy2,
                                create_masker_repeats_by_repeatmasker,
                                configure_pipeline, pipeline_runner,
                                checkpoint)


DATA_DIR = os.path.join(os.path.split(biolib.__path__[0])[0], 'data')


class SeqCleanerTest(unittest.TestCase):
    'It tests cleaner function from seq_cleaner'

    @staticmethod
    def test_strip_seq_by_quality():
        'test trim_seq_by_quality '
        qual = [20, 20, 20, 60, 60, 60, 60, 60, 20, 20, 20, 20]
        seq  = 'ataataataata'
        seq1 = SeqWithQuality(qual=qual, seq=seq)
        strip_seq_by_quality = create_striper_by_quality(quality_treshold=40,
                                                         min_seq_length=2,
                                                         min_quality_bases=3)
        new_seq = strip_seq_by_quality(seq1)
        assert new_seq.seq == 'ataat'

        qual = [60, 60, 60, 60, 60, 60, 60]
        seq  = 'ataataa'
        new_seq = strip_seq_by_quality(SeqWithQuality(qual=qual, seq=seq))
        assert  new_seq.seq == 'ataataa'

        qual = [60, 60, 60, 60, 60, 60, 0]
        seq  = 'ataataa'
        new_seq = strip_seq_by_quality(SeqWithQuality(qual=qual, seq=seq))
        assert new_seq.seq == 'ataata'
        assert new_seq.qual == [60, 60, 60, 60, 60, 60]

    @staticmethod
    def test_mask_low_complexity():
        'It test mask_low_complexity function'
        seq = 'TCGCATCGATCATCGCAGATCGACTGATCGATCGATCGGGGGGGGGGGGGGGGGGGGGGGG'
        qual = [30] * 61
        seq1 = SeqWithQuality(seq=seq, qual=qual)
        mask_low_complexity = create_masker_for_low_complexity()
        masked_seq = mask_low_complexity(seq1)
        assert masked_seq.seq[-10:] == 'gggggggggg'
        assert len(masked_seq.qual) == 61

        seq  = 'GGGGGTTTCTTAAATTCGCCTGGAGATTTCATTCGGGGGGGGGGTTCTCCCCAGGGGGGGGTG'
        seq += 'GGGAAACCCCCCGTTTCCCCCCCCGCGCGCCTTTTCGGGGAAAATTTTTTTTTGTTCCCCCCG'
        seq += 'GAAAAAAAAATATTTCTCCTGCGGGGCCCCCGCGAAGAAAAAAGAAAAAAAAAAAGAGGAGGA'
        seq += 'GGGGGGGGGGGGCGAAAATATAGTTTGG'
        seq1 = SeqWithQuality(seq=seq)
        masked_seq = mask_low_complexity(seq1)
        expected =  'GGGGGTTTCTTAAATTCGCCTGGAGATTTCATtcggggggggggttctccccaggggg'
        expected += 'gggtggggAAaccccccgtttccccccccgcgcgccttttcggggaaaattttttttt'
        expected += 'gttccccccGGAAAAAAAAATATTTCTCCTGCGGGGCCCCCGCGaagaaaaaagaaaa'
        expected += 'aaaaaaaGAGGAGGAGGGgggggggggCGAAAATATAGTTTGG'

        assert  str(masked_seq) == expected


    @staticmethod
    def test_mask_polya():
        'It test mask_polyA function'
        seq = 'TCGCATCGATCATCGCAGATCGACTGATCGATCGATCAAAAAAAAAAAAAAAAAAAAAAA'
        seq1 = SeqWithQuality(seq=seq)
        mask_polya = create_masker_for_polia()
        masked_seq = mask_polya(seq1)
        exp_seq = 'TCGCATCGATCATCGCAGATCGACTGATCGATCGATCaaaaaaaaaaaaaaaaaaaaaaa'
        assert masked_seq.seq == exp_seq

    def test_trim_seq_by_qual_trimpoly(self):
        'It test trimpoly  but with trim low quality parameters'
        seq  = 'ATCGATCTGATCTAGTCGATGTCTAGCTGAGCTACATAGCTAACGATCTAGTCTAGTCTATG'
        seq += 'TCATGTCATGTCGATGTCTAGTCTAGTCTAGTGAGTCACTGACTAGATCATGACATCGANNN'
        seq += 'NNNNNNNNNNNNNNNNNNTACTAGTC'
        qual = [10] * 150
        seq1 = SeqWithQuality(seq=seq, qual=qual)
        strip_seq_by_quality_trimpoly = create_striper_by_quality_trimpoly()
        trimmed_seq = strip_seq_by_quality_trimpoly(seq1)
        # It does not mask anything with this example, but it check that
        # if works
        assert trimmed_seq.seq.endswith('ATGACATCGA')



    def xtest_strip_seq_by_qual_lucy(self):
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
        strip_seq_by_quality_lucy = create_striper_by_quality_lucy()
        striped_seq = strip_seq_by_quality_lucy(seqrec)

        seq = striped_seq.seq
        assert seq.startswith('CAGATCAGATCAGCATCAGCAT')
        assert seq.endswith('CGAGATCAGCAGCATCAGC')

        try:
            #seq with no name
            seqrec = SeqWithQuality(seq=seq, qual=qual)
            strip_seq_by_quality_lucy(seqrec)
            self.fail("Value Error expected")
            #pylint: disable-msg=W0704
        except ValueError:
            pass

        # this is to tune the lucy configuration, It should cut until the
        # quality starts to
        qual  = '40 40 40 37 40 40 37 37 37 37 37 37 37 37 40 42 42 42 44 44 '
        qual += '56 56 42 40 40 40 40 36 36 28 35 32 35 35 40 42 37 37 35 37 '
        qual += '32 35 35 35 35 35 35 38 33 33 24 33 33 42 33 35 35 35 35 33 '
        qual += '36 30 30 24 29 29 35 35 35 35 29 29 29 35 38 38 38 37 35 33 '
        qual += '29 35 35 34 30 30 30 33 29 31 31 29 29 29 28 28 24 21 16 16 '
        qual += '21 24 29 29 32 40 27 27 25 25 21 30 27 28 28 32 23 23 21 24 '
        qual += '24 17 18 19 21 15 19 11 9 9 11 23 17 15 10 10 10 20 27 25 23 '
        qual += '18 22 23 24 18 10 10 13 13 18 19 10 12 12 18 16 14 10 10 11 '
        qual += '16 13 21 19 31 19 27 27 28 26 29 25 25 20 19 23 28 28 19 20 '
        qual += '13 9 9 9 9 9 17 15 21 17 14 12 21 17 19 24 28 24 23 '
        quality = qual.split()
        seq =  'ATCGATCAGTCAGACTGACAGACTCAGATCAGATCAGCATCAGCATACGATACGCATCAGACT'
        seq += 'ACGATCGATCGATCGACAGATCATCGATCATCGACGACTAGACGATCATCGATACGCAGACTC'
        seq += 'AGCAGACTACGAGATCAGCAGCATCAGCAGCAAGCAGACTACGAGATCAGCAGCATCAGCAGC'
        seq += 'ATTACGATGAT'
        seq1 = SeqWithQuality(seq=seq, qual=quality, name='seq1')
        striped_seq = strip_seq_by_quality_lucy(seq1)
        assert len(striped_seq.qual) > 170
        assert len(striped_seq.qual) < 185

    @staticmethod
    def test_strip_seq_by_quality_lucy2():
        'It tests strip_seq_by_quality_liucy2'
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
        seqrec1 = SeqWithQuality(name='seq1', seq=seq, qual=qual)

        qual  = '40 40 40 37 40 40 37 37 37 37 37 37 37 37 40 42 42 42 44 44 '
        qual += '56 56 42 40 40 40 40 36 36 28 35 32 35 35 40 42 37 37 35 37 '
        qual += '32 35 35 35 35 35 35 38 33 33 24 33 33 42 33 35 35 35 35 33 '
        qual += '36 30 30 24 29 29 35 35 35 35 29 29 29 35 38 38 38 37 35 33 '
        qual += '29 35 35 34 30 30 30 33 29 31 31 29 29 29 28 28 24 21 16 16 '
        qual += '21 24 29 29 32 40 27 27 25 25 21 30 27 28 28 32 23 23 21 24 '
        qual += '24 17 18 19 21 15 19 11 9 9 11 23 17 15 10 10 10 20 27 25 23 '
        qual += '18 22 23 24 18 10 10 13 13 18 19 10 12 12 18 16 14 10 10 11 '
        qual += '16 13 21 19 31 19 27 27 28 26 29 25 25 20 19 23 28 28 19 20 '
        qual += '13 9 9 9 9 9 17 15 21 17 14 12 21 17 19 24 28 24 23 '
        quality = qual.split()
        seq =  'ATCGATCAGTCAGACTGACAGACTCAGATCAGATCAGCATCAGCATACGATACGCATCAGACT'
        seq += 'ACGATCGATCGATCGACAGATCATCGATCATCGACGACTAGACGATCATCGATACGCAGACTC'
        seq += 'AGCAGACTACGAGATCAGCAGCATCAGCAGCAAGCAGACTACGAGATCAGCAGCATCAGCAGC'
        seq += 'ATTACGATGAT'
        seqrec2 = SeqWithQuality(seq=seq, qual=quality, name='seq2')
        seq_iter = iter([seqrec1, seqrec2])
        lucy_striper = create_striper_by_quality_lucy2()
        #pylint:disable-msg=W0612
        seq_iter, output_files = lucy_striper(seq_iter)
        seqs = list(seq_iter)
        seq = seqs[0].seq
        assert seq.startswith('CAGATCAGATCAGCATCAGCAT')
        assert seq.endswith('CGAGATCAGCAGCATCAGC')
        assert len(seqs) == 2

    @staticmethod
    def test_strip_vector_align_exonera():
        'It tests strip_vector_by_alignment'

        vec1 = SeqWithQuality(name='vec1', seq='atcgatcgatagcatacgat')
        vec2 = SeqWithQuality(name='vec2', seq='atgcatcagatcgataaaga')
        fhand_vectors = temp_multi_fasta_file([vec1, vec2])

        strip_vector_by_alignment = \
                create_vector_striper_by_alignment(fhand_vectors, 'exonerate')

        seq  = 'ATGCATCAGATGCATGCATGACTACGACTACGATCAGCATCAGCGATCAGCATCGATACGATC'
        seq2 = SeqWithQuality(name='seq', seq=seq)
        seq1  = SeqWithQuality(name=seq2.name, seq=vec1.seq+seq2.seq+vec2.seq)

        seq3 = strip_vector_by_alignment(seq1)
        assert str(seq2.seq) == str(seq3.seq)


        fhand_vectors.seek(0)
        seq1  = SeqWithQuality(name=seq2.name, seq=vec1.seq+vec2.seq+seq2.seq )
        seq3 = strip_vector_by_alignment(seq1)
        assert str(seq2.seq) == str(seq3.seq)

        # overlaping vectors
        fhand_vectors.seek(0)
        seq1  = SeqWithQuality(name=seq2.name, seq=vec1[:-2]+vec2+seq2+vec2)
        seq3 = strip_vector_by_alignment(seq1)
        assert str(seq2.seq) == str(seq3.seq)

        # Now only vectors
        fhand_vectors.seek(0)
        seq1 = SeqWithQuality(name=seq2.name, seq=vec1+vec2+vec2)
        seq3 = strip_vector_by_alignment(seq1)
        assert seq3 is None

        # with some extra seq at the begining and end
        fhand_vectors.seek(0)
        seq1 = SeqWithQuality(name=seq2.name,
                              seq=seq2[:20]+vec1+seq2+vec2+seq2[:20] )
        seq3 = strip_vector_by_alignment(seq1)
        assert str(seq2.seq) == str(seq3.seq)

        # Now without vectors
        fhand_vectors.seek(0)
        seq1 = seq2
        seq3 = strip_vector_by_alignment(seq1)
        assert str(seq2.seq) == str(seq3.seq)

        fhand_vectors.seek(0)
        seq1  = SeqWithQuality(name=seq2.name, seq=vec1[::-1]+vec2+seq2 )
        seq3 = strip_vector_by_alignment(seq1)
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
        strip_vector_by_alignment = \
                            create_vector_striper_by_alignment(vector, 'blast')
        striped_seq = strip_vector_by_alignment(seq)
        striped_seq = str(striped_seq.seq)

        assert seq1[4:14]  in striped_seq
        assert seq1[-14:-4]  in striped_seq
        assert vec1[4:14]  not in striped_seq
        assert vec1[-14:-4] not  in striped_seq

    @staticmethod
    def xtest_repeatmasking():
        'It test that we can mask a repeat element using repeat masker'
        mask_repeats_by_repeatmasker = \
                 create_masker_repeats_by_repeatmasker(species='eudicotyledons')

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


        #This test with no repetitive regions
        seq  = 'ATCGATCTGATCTAGTCGATGTCTAGCTGAGCTACATAGCTAACGATCTAGTCTAGTCTATGA'
        seq += 'TGCATCAGATGCATGAAATCGATCTGATCTAGTCGATGTCTAGCTGAGCTACATAGCTAACGA'
        seq += 'TCTAGTCTAGTCTATGATGCATCAGCTACGATGATCATGTCATGTCGATGTCTAGTCTAGTCT'
        seq += 'AGTGAGTCACTGACTAGATCATGACATCGATACTAGTC'
        seqrec  = SeqWithQuality(name='seq', seq=seq)
        masked_seq = mask_repeats_by_repeatmasker(seqrec)

        masked_str = str(masked_seq.seq)
        assert  masked_str == seq

ADAPTORS = '''>adaptor1
atcgatcgatagcatacgat
>adaptor2
atgcatcagatcgataaaga'''

EXPECTED = '''>seq1
ATCGATCTGATCTAGTCGATGTCTAGCTGAGCTACATAGCTAACGATCTAGTCTAGTCTATGATGCATCAGATGCATGAATCGCATCGATCATCGCAGATCGACTGATCGATATGCATCAGATCGCGATCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAATCGATCTGATCTAGTCGATGTCTAGCTGAGCTACATAGCTAACGATCTAGTCTAGTCTATGATGCATCAGATGCATGAA
>seq2
ATCGATCTGATCTAGTCGATGTCTAGCTGAGCTACATAGCTAACGATCTAGTCTAGTCTATGATGCATCAGATGCATGAAATCGATCTGATCTAGTCGATGTCTAGCTGAGCTACATAGCTAACGATCTAGTCTAGTCTATGATGCATCAGATGCATGAAATCAGCATGACTCATCGCATCGATCATCGCAGATCGACTGATCGATCGATCAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>seq3
ATCGATCTGATCTAGTCGATGTCTAGCTGAGCTACATAGCTAACGATCTAGTCTAGTCTATGATGCATCAGATGCATGAAATCGATCTGATCTAGTCGATGTCTAGCTGAGCTACATAGCTAACGATCTAGTCTAGTCTATGATGCATCAGCTACGATGATCATGTCATGTCGATGTCTAGTCTAGTCTAGTGAGTCACTGACTAGATCATGACATCGANNNNNNNNNNNNNNNNNNNNNNTACTAGTC
>seq4
ATCGATCAGTCAGACTGACAGACTCAGATCAGATCAGCATCAGCATACGATACGCATCAGACTTCAGCATCGATCGACTAACGATCGATCGATCGACAGATCATCGATCATCGACGACTAGACGATCATCGATACGCAGACTCCGACTACGACTACGATAAGCAGACTACGAGATCAGCAGCATCAGCAGCANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
>seq5
ATGCATCAGATGCATGCATGACTACGACTACGATCAGCATCAGCGATCAGCATCGATACGATCATCGACTGCATCGATGAATCGATCAGTCAGACTGACAGACTCAGATCAGATCAGCATCAGCATACGATACGCATCAGACTTCAGCATCGATCGACTAATCGATCAGTCAGACTGACAGACTCAGATCAGATCAGCATCAGCATACGATACGCATCAGACTTCAGCATCGATCGACTA
>seq6
AACCGTTTGACTTACGATATTTGCCCATTGTGATTCTAGTCGATTTGCATAACGTGTACGTATCGGTATTGTGACTGATTCGATGCTATTGCAAACAGTTTTGATTGTGTGATCGTGATGCATGCTAGTCTGATCGAGTCTGATCGTAGTCTAGTCGTAGTCGATGTCGATTTATCAGTAGTCGATGCTAGTCTAGTCTAGTCTACTAGTCTAGTCATGCTAGTCGAGTCGAT
'''


class PipelineTests(unittest.TestCase):
    'It test pipeline related functions'

    def test_configure_pipeline(self):
        'It tests configure pipeline'
        pipeline      = 'sanger_with_quality_clean'
        configuration = {'remove_vectors': {'vectors':'Univec'},
                         'remove_adaptors':{'vectors':'hola'}}
        pipeline      = configure_pipeline(pipeline, configuration)

        assert pipeline[0]['arguments']['vectors'] == 'Univec'

        # Now it should fail because one of the arguments is Not set
        configuration = {'remove_vectors': {'vectors':'Univec'}}
        try:
            pipeline = configure_pipeline(pipeline, configuration)
            self.fail()
            #pylint: disable-msg=W0704
        except Exception:
            pass

    @staticmethod
    def test_pipeline_run():
        'It tests that the pipeline runs ok'
        pipeline = 'sanger_with_quality_clean'

        fhand_adaptors = NamedTemporaryFile()
        fhand_adaptors.write(ADAPTORS)
        fhand_adaptors.flush()

        configuration = {'remove_vectors': {'vectors':'UniVec'},
                         'remove_adaptors':{'vectors':fhand_adaptors.name}}

        io_fhands = {}
        io_fhands['in_seq']   = open(os.path.join(DATA_DIR, 'seq.fasta'), 'r')
        io_fhands['in_qual']  = open(os.path.join(DATA_DIR, 'qual.fasta'), 'r')
        io_fhands['out_seq']  = NamedTemporaryFile()
        io_fhands['out_qual'] = NamedTemporaryFile()

        dir_ = NamedTemporaryDir()
        work_dir = dir_.get_name()

        pipeline_runner(pipeline, configuration, io_fhands, work_dir)
        io_fhands['out_seq'].seek(0)
        result_seq = io_fhands['out_seq'].read()
        assert result_seq.count('>') == 6


class CheckPointTest(unittest.TestCase):
    'It checks checkpoint function'
    @staticmethod
    def test_checkpoint():
        'It checks checkpoint function'
        seq  = 'atgatgatgt'
        qual = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        seq1 = SeqWithQuality(seq=seq, qual=qual, name='seq1')
        seq2 = SeqWithQuality(seq=seq, qual=qual, name='seq2')
        seqs_iter = iter([seq1, seq2])


        fhand_seq  = NamedTemporaryFile()
        fhand_qual = NamedTemporaryFile()
        seq_iter2 = checkpoint(seqs_iter, fhand_seq, fhand_qual)
        assert str(seq_iter2.next().seq) == seq

        fhand_seqs_out = open(fhand_seq.name)
        fhand_qual_out = open(fhand_qual.name)

        assert fhand_seqs_out.readline()[0] == '>'
        assert fhand_seqs_out.readline()[0] == 'a'
        assert fhand_qual_out.readline()[0] == '>'
        assert fhand_qual_out.readline()[0] == '0'

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
