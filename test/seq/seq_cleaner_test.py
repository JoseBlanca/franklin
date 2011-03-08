'''
Created on 2009 uzt 6

@author: peio
'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of project.
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

import unittest, os

from franklin.seq.seqs import SeqWithQuality, Seq
from franklin.seq.writers import temp_fasta_file
from franklin.seq.readers import seqs_in_file
from franklin.seq.seq_cleaner import (create_vector_striper_by_alignment,
                                create_masker_for_polia,
                                create_masker_for_low_complexity,
                                create_striper_by_quality_trimpoly,
                                create_striper_by_quality,
                                create_striper_by_quality_lucy,
                                _get_non_matched_locations,
                                _get_unmasked_locations,
                                _get_matched_locations,
                                split_seq_by_masked_regions,
                                create_word_masker,
                                create_edge_stripper,
                                create_word_striper_by_alignment,
                                create_upper_mapper,
                                create_seq_trim_and_masker,
                                _mask_sequence, TRIMMING_RECOMENDATIONS,
                                _get_all_segments,
                                _get_non_matched_from_matched_locations)

from franklin.utils.misc_utils import DATA_DIR
from franklin.utils.cmd_utils import BLAST_TOOL

class SeqCleanerTest(unittest.TestCase):
    'It tests cleaner function from seq_cleaner'

    @staticmethod
    def test_upper():
        'It tests the upper mapper'
        seq = SeqWithQuality(seq=Seq('actg'))
        upper = create_upper_mapper()
        seq = upper(seq)
        assert seq.seq == 'ACTG'

    @staticmethod
    def test_strip_seq_by_quality():
        'test trim_seq_by_quality '
        qual = [20, 20, 20, 60, 60, 60, 60, 60, 20, 20, 20, 20]
        seq  = 'ataataataata'
        desc = 'hola'
        seq1 = SeqWithQuality(qual=qual, seq=Seq(seq), description=desc)
        strip_seq_by_quality = create_striper_by_quality(quality_treshold=40,
                                                         min_seq_length=2,
                                                         min_quality_bases=3)
        seq_trimmer = create_seq_trim_and_masker()
        new_seq = strip_seq_by_quality(seq1)
        new_seq = seq_trimmer(new_seq)
        assert new_seq.seq == 'ataat'
        assert new_seq.description == desc

        qual = [60, 60, 60, 60, 60, 60, 60]
        seq  = 'ataataa'
        new_seq = strip_seq_by_quality(SeqWithQuality(qual=qual, seq=Seq(seq)))
        new_seq = seq_trimmer(new_seq)
        assert  new_seq.seq == 'ataataa'

        qual = [60, 60, 60, 60, 60, 60, 0]
        seq  = 'ataataa'
        new_seq = strip_seq_by_quality(SeqWithQuality(qual=qual, seq=Seq(seq)))
        new_seq = seq_trimmer(new_seq)
        assert new_seq.seq == 'ataata'
        assert new_seq.qual == [60, 60, 60, 60, 60, 60]

        qual = [40, 18, 10, 40, 40, 5, 8, 30, 14, 3, 40, 40, 40, 11, 6, 5, 3,
               20, 10, 12, 8, 5, 4, 7, 1]
        seq = 'atatatatagatagatagatagatg'
        strip_seq_by_quality = create_striper_by_quality(quality_treshold=20,
                                                         min_seq_length=2,
                                                         min_quality_bases=3,
                                                         quality_window_width=2)
        new_seq = strip_seq_by_quality(SeqWithQuality(qual=qual, seq=Seq(seq)))
        new_seq = seq_trimmer(new_seq)
        assert new_seq.qual == [40, 18, 10, 40, 40, 5, 8, 30, 14, 3, 40, 40, 40,
                                11]


        qual = [40, 40, 13, 11, 40, 9, 40, 4, 27, 38, 40, 4, 11, 40, 40, 10, 10,
                21, 3, 40, 9, 9, 12, 10, 9]
        seq  = 'atatatatatatatatatatatata'
        strip_seq_by_quality = create_striper_by_quality(quality_treshold=20,
                                                         min_seq_length=2,
                                                         min_quality_bases=3,
                                                         quality_window_width=1)
        new_seq = strip_seq_by_quality(SeqWithQuality(qual=qual, seq=Seq(seq)))
        new_seq = seq_trimmer(new_seq)
        assert new_seq.qual == [40, 40, 13, 11, 40, 9, 40, 4, 27, 38, 40]

    @staticmethod
    def test_mask_low_complexity():
        'It test mask_low_complexity function'
        seq = 'TCGCATCGATCATCGCAGATCGACTGATCGATCGATCGGGGGGGGGGGGGGGGGGGGGGGG'
        seq = Seq(seq)
        qual = [30] * 61
        desc = 'hola'
        seq1 = SeqWithQuality(seq=seq, qual=qual, description=desc)
        mask_low_complexity = create_masker_for_low_complexity()
        masked_seq = mask_low_complexity(seq1)
        sequence_trimmer = create_seq_trim_and_masker()
        masked_seq = sequence_trimmer(masked_seq)
        assert masked_seq.seq[-10:] == 'gggggggggg'
        assert len(masked_seq.qual) == 61
        assert masked_seq.description == desc

        seq  = 'GGGGGTTTCTTAAATTCGCCTGGAGATTTCATTCGGGGGGGGGGTTCTCCCCAGGGGGGGGTG'
        seq += 'GGGAAACCCCCCGTTTCCCCCCCCGCGCGCCTTTTCGGGGAAAATTTTTTTTTGTTCCCCCCG'
        seq += 'GAAAAAAAAATATTTCTCCTGCGGGGCCCCCGCGAAGAAAAAAGAAAAAAAAAAAGAGGAGGA'
        seq += 'GGGGGGGGGGGGCGAAAATATAGTTTGG'
        seq1 = SeqWithQuality(seq=Seq(seq))
        masked_seq = mask_low_complexity(seq1)
        masked_seq = sequence_trimmer(masked_seq)
        expected =  'GGGGGTTTCTTAAATTCGCCTGGAGATTTCATtcggggggggggttctccccaggggg'
        expected += 'gggtggggAAaccccccgtttccccccccgcgcgccttttcggggaaaattttttttt'
        expected += 'gttccccccGGAAAAAAAAATATTTCTCCTGCGGGGCCCCCGCGaagaaaaaagaaaa'
        expected += 'aaaaaaaGAGGAGGAGGGgggggggggCGAAAATATAGTTTGG'

        assert  str(masked_seq.seq) == expected


    @staticmethod
    def test_mask_polya():
        'It test mask_polyA function'
        seq = 'TCGCATCGATCATCGCAGATCGACTGATCGATCGATCAAAAAAAAAAAAAAAAAAAAAAA'
        seq = Seq(seq)
        seq1 = SeqWithQuality(seq=seq, description='hola')
        mask_polya = create_masker_for_polia()
        masked_seq = mask_polya(seq1)
        sequence_trimmer = create_seq_trim_and_masker()
        masked_seq = sequence_trimmer(masked_seq)
        exp_seq = 'TCGCATCGATCATCGCAGATCGACTGATCGATCGATCaaaaaaaaaaaaaaaaaaaaaaa'
        assert masked_seq.seq == exp_seq
        assert masked_seq.description == 'hola'

    def test_trim_seq_by_qual_trimpoly(self):
        'It test trimpoly  but with trim low quality parameters'
        seq  = 'ATCGATCTGATCTAGTCGATGTCTAGCTGAGCTACATAGCTAACGATCTAGTCTAGTCTATG'
        seq += 'TCATGTCATGTCGATGTCTAGTCTAGTCTAGTGAGTCACTGACTAGATCATGACATCGANNN'
        seq += 'NNNNNNNNNNNNNNNNNNTACTAGTC'
        seq = Seq(seq)
        qual = [10] * 150
        desc = 'hola'
        seq1 = SeqWithQuality(seq=seq, qual=qual, description=desc)
        strip_seq_by_quality_trimpoly = create_striper_by_quality_trimpoly()
        trimmed_seq = strip_seq_by_quality_trimpoly(seq1)
        sequence_trimmer = create_seq_trim_and_masker()
        trimmed_seq = sequence_trimmer(trimmed_seq)
        # It does not mask anything with this example, but it check that
        # if works
        assert trimmed_seq.seq.endswith('ATGACATCGA')

        #another seq with the problem at the begining
        seq  = 'TNNNNNNAGGGCTTTCCTGACAGCTANNNNNTTTGCGGGCAACATCCAGAACAAGCACCG'
        seq += 'GCAGATTGGCAATGCCGTGCCCCCGCCTCTTGCCTATGCACTTGGGAGGAAGCTGAAGGA'
        seq += 'AGCCGTTGACAAGCGTCAGGAAGCCAGCGCAGGCGTGCCTGCACCATGAGAAGTTTTCCT'
        seq = Seq(seq)
        seq1 = SeqWithQuality(seq=seq)
        strip_seq_by_quality_trimpoly = create_striper_by_quality_trimpoly()
        trimmed_seq = strip_seq_by_quality_trimpoly(seq1)
        trimmed_seq = sequence_trimmer(trimmed_seq)
        str_seq = str(trimmed_seq.seq)
        assert str_seq.startswith('TTTGCGGGCAACA')
        assert str_seq.endswith('GAGAAGTTTTCCT')

        seq  = 'TGACATCGAACCTCGGCGCCGAGCACCTCCTCGCTGGGATGGTGGGCAAGAACTCCATGA'
        seq += 'AGGTCGCTCGCGATCTGGTCATGCAGGAGGTGAGGAGGCACTTCCGCCCTGAGCTGCTGA'
        seq += 'ACCGTCTCGACGAGATCGTGATCTTCGATCCTCTGTCCCACGAGCAGCTGAGGAAGGTCG'
        seq += 'CTCGCCTTCAGATGAAGGATGTGGCCGTCCGTCTTGCCGAANNNNNCATCGCTCTGGCTG'
        seq += 'TGACCGANNNNNCATTGGACATCATCTTGTCTCTCTCTNNNNNNTCNNNNT'
        seq1 = SeqWithQuality(seq=Seq(seq))
        strip_seq_by_quality_trimpoly = create_striper_by_quality_trimpoly()
        trimmed_seq = strip_seq_by_quality_trimpoly(seq1)
        trimmed_seq = sequence_trimmer(trimmed_seq)
        str_seq = str(trimmed_seq.seq)
        assert str_seq.startswith('TGACATCGAA')
        assert str_seq.endswith('TCTTGCCGAA')

        seq  = 'TACGGCCGGGGTNNCNNANNNNGCATTCTCGCAGGGTCTTTCTACACTATTAGATAAGAT'
        seq += 'GGATCCTTCTCAGAGAGTGAAGTTTGTTCAGGAAGTCAAGAAGGTTCTTGGATGATGATA'
        seq += 'TGATACCAACACATCCAACACAATATGCGCATGCTACATGTTATTTTTCAAGTACATACA'
        seq += 'TAGAAGGATATTGCTTGGCCTTGATTGATCATGTCTGATCTAAGTCGATCATTATTTTCT'
        seq += 'TGAAACTTCCTTTCGGACGTGGTGCTATGGTTGATGAATTTGGATGTGTGCGTTCTGCCA'
        seq += 'GGTGTAAGCCCAAAGGTTTATACAGACCGAGTTAAGGTTAGGAAGAGCACGAGTGAACTT'
        seq1 = SeqWithQuality(seq=Seq(seq))
        strip_seq_by_quality_trimpoly = create_striper_by_quality_trimpoly()
        trimmed_seq = strip_seq_by_quality_trimpoly(seq1)
        trimmed_seq = sequence_trimmer(trimmed_seq)
        str_seq = str(trimmed_seq.seq)
        assert str_seq.startswith('GCATTCTCGCAG')
        assert str_seq.endswith('GTGAACTT')


    @staticmethod
    def test_strip_seq_by_quality_lucy():
        'It tests strip_seq_by_quality_lucy2'
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
        seq = Seq(seq)
        seqrec1 = SeqWithQuality(name='seq1', seq=seq, qual=qual,
                                  description ='desc1')

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
        seq = Seq(seq)
        seqrec2 = SeqWithQuality(seq=seq, qual=quality, name='seq2',
                                 description ='desc2')
        seq_iter = iter([seqrec1, seqrec2])
        seq_trimmer = create_seq_trim_and_masker()
        lucy_striper = create_striper_by_quality_lucy()
        #pylint:disable-msg=W0612
        seq_iter = lucy_striper(seq_iter)
        new_seqs = []
        for seq in seq_iter:
            new_seqs.append(seq_trimmer(seq))
        seq = new_seqs[0].seq
        assert seqrec1.description == new_seqs[0].description
        assert seq.startswith('CAGATCAGATCAGCATCAGCAT')
        assert seq.endswith('CGAGATCAGCAGCATCAGC')
        assert len(new_seqs) == 2
        assert new_seqs[1].description == 'desc2'

        # now we test the sequence with adaptors
        vector_fpath = os.path.join(DATA_DIR, 'lucy', 'icugi_vector.fasta')
        splice_fpath = os.path.join(DATA_DIR, 'lucy', 'icugi_splice.fasta')
        parameters = {'vector':[vector_fpath, splice_fpath],
                      'bracket':[10, 0.02]}
        lucy_striper = create_striper_by_quality_lucy(parameters)
        seq_fhand = open(os.path.join(DATA_DIR, 'lucy',
                                      'seq_with_adaptor1.fastq'))
        seq_iter = lucy_striper(seqs_in_file(seq_fhand, format='fastq'))
        new_seqs = []
        for seq in seq_iter:
            new_seqs.append(seq_trimmer(seq))

    @staticmethod
    def test_strip_vector_exonerate():
        'It tests strip_vector_by_alignment'

        vec1 = SeqWithQuality(name='vec1', seq=Seq('atcgatcgatagcatacgat'))
        vec2 = SeqWithQuality(name='vec2', seq=Seq('atgcatcagatcgataaaga'))
        fhand_vectors = temp_fasta_file([vec1, vec2])
        seq_trimmer = create_seq_trim_and_masker()
        strip_vector_by_alignment = \
                create_vector_striper_by_alignment(fhand_vectors, 'exonerate')

        seq  = 'ATGCATCAGATGCATGCATGACTACGACTACGATCAGCATCAGCGATCAGCATCGATACGATC'
        seq = Seq(seq)
        seq2 = SeqWithQuality(name='seq', seq=seq)
        seq1 = SeqWithQuality(name=seq2.name,
                              seq=vec1.seq + seq2.seq + vec2.seq,
                              description='hola')

        seq3 = strip_vector_by_alignment(seq1)
        seq3 = seq_trimmer(seq3)

        assert str(seq2.seq) == str(seq3.seq)
        assert seq3.description == 'hola'

        fhand_vectors.seek(0)
        seq1  = SeqWithQuality(name=seq2.name, seq=vec1.seq+vec2.seq+seq2.seq)
        seq3 = strip_vector_by_alignment(seq1)
        seq3 = seq_trimmer(seq3)
        assert str(seq2.seq) == str(seq3.seq)

        # overlaping vectors
        fhand_vectors.seek(0)
        new_seq = vec1.seq[:-2]+vec2.seq+seq2.seq+vec2.seq
        seq1  = SeqWithQuality(name=seq2.name, seq=new_seq)
        seq3 = strip_vector_by_alignment(seq1)
        seq3 = seq_trimmer(seq3)
        assert str(seq2.seq) == str(seq3.seq)

        # Now only vectors
        fhand_vectors.seek(0)
        new_seq = vec1.seq+vec2.seq+vec2.seq
        seq1 = SeqWithQuality(name=seq2.name, seq=new_seq)
        seq3 = strip_vector_by_alignment(seq1)
        seq3 = seq_trimmer(seq3)
        assert seq3 is None

        # with some extra seq at the begining and end
        fhand_vectors.seek(0)
        seq1 = SeqWithQuality(name=seq2.name,
                    seq=seq2.seq[:20]+vec1.seq+seq2.seq+vec2.seq+seq2.seq[:20])
        seq3 = strip_vector_by_alignment(seq1)
        seq3 = seq_trimmer(seq3)
        assert str(seq2.seq) == str(seq3.seq)

        # Now without vectors
        fhand_vectors.seek(0)
        seq1 = seq2
        seq3 = strip_vector_by_alignment(seq1)
        seq3 = seq_trimmer(seq3)
        assert str(seq2.seq) == str(seq3.seq)

        fhand_vectors.seek(0)
        seq1  = SeqWithQuality(name=seq2.name,
                               seq=vec1.seq[::-1]+vec2.seq+seq2.seq)
        seq3 = strip_vector_by_alignment(seq1)
        seq3 = seq_trimmer(seq3)
        assert str(seq2.seq) == str(seq3.seq)

    @staticmethod
    def test_strip_short_adaptors():
        'It tests the short adaptor removal with J. Forment sequences'

        seq  = 'CgCGTGTCTCTAgATAGGGACAGTAGGAATCTCGTTAATCCATTCATGCGCGTCACTAATTAG'
        seq += 'ATGACGAGGCATTTGGCTACCTTAAGAGAGTCATAGTTACTCCCGCCGTTTA'
        seq  = Seq(seq)
        seq  = SeqWithQuality(name='seq', seq=seq)
        seq_trimmer = create_seq_trim_and_masker()
        strip_adap = create_word_striper_by_alignment(words=['CGTGTCTCTA',
                                                             'TATATATA'])

        clean_seq = strip_adap(seq)
        clean_seq = seq_trimmer(clean_seq)
        assert str(clean_seq.seq).startswith('gATAGGGACAGTAGGAATCTCGTTAATC')

        #two words
        #       000000000011111111112222222222333333333344444444
        #       012345678901234567890123456789012345678901234567
        seq  = 'CgCGTGTCTCTAgATAGGGACAGTAGGAATTTTTTTcCGTGTCTCTAc'
        seq  = SeqWithQuality(name='seq', seq=Seq(seq))
        clean_seq = strip_adap(seq)
        clean_seq = seq_trimmer(clean_seq)
        assert str(clean_seq.seq) == 'gATAGGGACAGTAGGAATTTTTTTc'

        #We don't want the first region even if it's the longest one
        seq  = 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCgCGTGTCTCTAgC'
        seq  = SeqWithQuality(name='seq', seq=Seq(seq))
        clean_seq = strip_adap(seq)
        clean_seq = seq_trimmer(clean_seq)
        assert str(clean_seq.seq) == 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCg'

        #Everything works if the adaptors are at the beginning and ends
        seq  = 'CGTGTCTCTAcccccccccccccccTATATATAggggggggggCGTGTCTCTA'
        seq  = SeqWithQuality(name='seq', seq=Seq(seq))
        clean_seq = strip_adap(seq)
        clean_seq = seq_trimmer(clean_seq)
        assert str(clean_seq.seq) == 'ccccccccccccccc'

        #It test if we remove words from the beginning of the seq
        word = 'ATAT'
        seq2 = 'tctcatcatca'
        seq1 = word + seq2 + word
        word = '^' + word
        seq1 = Seq(seq1)
        seq  = SeqWithQuality(seq1, qual=[30] * len(seq1))

        remover = create_word_striper_by_alignment([word])
        seq = remover(seq)
        seq = seq_trimmer(seq)
        assert seq.seq == seq2

    @staticmethod
    def test_strip_vector_align_blast():
        'It tests strip_vector_by_alignment using blast and UniVec'
        vector = os.path.join(DATA_DIR, 'blast', 'univec+')
        vec1  = 'CTCGGGCCGTCTCTTGGGCTTGATCGGCCTTCTTGCGCATCTCACGCGCTCCTGCGGCGGCC'
        vec1 += 'TGTAGGGCAGGCTCATACCCCTGCCGAACCGCTTTTGTCAGCCGGTCGGCCACGGCTTCCGG'
        vec1 += 'CGTCTCAACGCGCTTT'
        seq1 = 'ATGCATCAGATGCATGCATGACTACGACTACGATCAGCATCAGCGATCAGCATCGATACGATC'
        seq  = SeqWithQuality(name='seq', seq=Seq(seq1+vec1))
        seq_trimmer = create_seq_trim_and_masker()
        strip_vector_by_alignment = \
                            create_vector_striper_by_alignment(vector,
                                                               BLAST_TOOL)
        striped_seq = strip_vector_by_alignment(seq)
        striped_seq = seq_trimmer(striped_seq)
        striped_seq = str(striped_seq.seq)

        assert seq1[4:14]  in striped_seq
        assert seq1[-14:-4]  in striped_seq
        assert vec1[4:14]  not in striped_seq
        assert vec1[-14:-4] not  in striped_seq

class SeqSplitterTests(unittest.TestCase):
    'Here we test seq splitter functions'
    @staticmethod
    def test_non_matches_fragment_detector():
        'it test  that the functions detects all the fragments of an alignments'
        alignments = [{'matches':[
                             {'match_parts':[{'query_start':1, 'query_end':10},
                                             {'query_start':6, 'query_end':17},
                                             {'query_start':3, 'query_end':18}]
                              }]}]
        locations = _get_non_matched_locations(alignments)
        assert locations[0] == (1, 10)
        assert locations[1] == (6, 17)
        assert locations[2] == (3, 18)

    @staticmethod
    def test_unmasked_seq_location_detector():
        'it detects the location of the unmasked seq regions'
        seq = SeqWithQuality(seq=Seq('AATTaaTTaaTTT'), name='seq')
        locations = _get_unmasked_locations(seq)
        assert locations[0] == (0, 3)
        assert locations[1] == (6, 7)
        assert locations[2] == (10, 12)

        seq = SeqWithQuality(seq=Seq('AATT'), name='seq')
        locations = _get_unmasked_locations(seq)
        assert locations[0] == (0, 3)

        seq = SeqWithQuality(seq=Seq('aatt'), name='seq')
        locations = _get_unmasked_locations(seq)
        assert not locations

    @staticmethod
    def test_get_matched_regions():
        'it tests get_matched_region function'
        seq1 = 'AATTaatAATTAATtctcTTCtctctctctctcGCGCGCGCGCCC'
        seq = SeqWithQuality(seq=Seq(seq1), name='seq')
        locations =  _get_unmasked_locations(seq)

        seqs = list(_get_matched_locations(seq, locations, 3))
        assert str(seqs[0].seq) == 'AATT'
        assert str(seqs[1].seq) == 'AATTAAT'
        assert str(seqs[2].seq) == 'TTC'
        assert str(seqs[3].seq) == 'GCGCGCGCGCCC'

        seqs = list(_get_matched_locations(seq, locations, 5))
        assert str(seqs[0].seq) == 'AATTAAT'
        assert str(seqs[1].seq) == 'GCGCGCGCGCCC'

    @staticmethod
    def test_strip_masked_section():
        'It tests the strip_masked_functions'
        seq1 = 'AATTaatAATTAATtctcTTCtctctctctctcGCGCGCGCGCCC'
        seq = SeqWithQuality(seq=Seq(seq1), name='seq1')
        seq2 = 'AATTaatAATTAATtctcTTCtctctctctctcGCGCGCGCGCCC'
        seq_ = SeqWithQuality(seq=Seq(seq2), name='seq2')

        seq_iter = iter([seq, seq_])

        new_seqs = list(split_seq_by_masked_regions(seq_iter, 5))
        assert str(new_seqs[0].seq) == 'AATTAAT'
        assert str(new_seqs[1].seq) == 'GCGCGCGCGCCC'
        assert str(new_seqs[2].seq) == 'AATTAAT'
        assert str(new_seqs[3].seq) == 'GCGCGCGCGCCC'

    @staticmethod
    def test_word_masker():
        'It test if we remove words from the beginning of the seq'
        word = 'ATAT'
        seq1 = word + 'tctcatcatca'.upper()
        seq1 = Seq(seq1)
        seq  = SeqWithQuality(seq1, qual=[30] * len(seq1))

        remover = create_word_masker([word])
        seq = remover(seq)
        sequence_trimmer = create_seq_trim_and_masker()
        seq = sequence_trimmer(seq)
        assert seq.seq[0] == 'a'

        seq1 = 'atactctcatcatca'.upper()
        seq  = SeqWithQuality(Seq(seq1), qual=[30] * len(seq1))
        seq = remover(seq)
        seq = sequence_trimmer(seq)
        assert seq.seq == seq1

        word = 'CA'
        seq1 = 'ATCATCATCATCA'
        seq  = SeqWithQuality(Seq(seq1), qual=[30] * len(seq1))
        remover = create_word_masker([word], False)
        seq = remover(seq)
        seq = sequence_trimmer(seq)
        assert seq.seq == 'ATcaTcaTcaTca'

    @staticmethod
    def test_edge_stripper():
        'It test if we remove edges of the seq'
        seq1 = 'gggtctcatcatcaggg'
        seq  = SeqWithQuality(Seq(seq1), qual=[30] * len(seq1))

        edge_stripperr = create_edge_stripper(left_length=3, right_length=3)
        sequence_trimmer = create_seq_trim_and_masker()
        seq = edge_stripperr(seq)
        seq = sequence_trimmer(seq)
        assert seq.seq == 'tctcatcatca'

    @staticmethod
    def test_sequence_stripper():
        'It can cut using trimming recomendations'
        seq1 = 'gggtctcatcatcaggg'.upper()
        seq = SeqWithQuality(Seq(seq1), qual=[30] * len(seq1),
                             annotations={TRIMMING_RECOMENDATIONS:{}})

        trim_rec = seq.annotations[TRIMMING_RECOMENDATIONS]
        seq_trimmer = create_seq_trim_and_masker()

        trim_rec['vector']  = [(0,3), (8, 12)]
        seq2 = seq_trimmer(seq)
        assert str(seq2.seq) == 'CTCA'

        trim_rec['vector']  = [(0, 0), (8, 12)]
        seq2 = seq_trimmer(seq)
        assert str(seq2.seq) == 'GGTCTCA'


        trim_rec['vector']  = [(0, 0), (8, 12)]
        trim_rec['quality'] = []
        trim_rec['mask']    = [(0, 3), (5, 6)]
        seq2 = seq_trimmer(seq)
        assert seq2.seq == 'ggtCtcA'
        assert 'vector' not in trim_rec
        assert 'quality' not in trim_rec
        assert 'mask' not in trim_rec

        seq_trimmer = create_seq_trim_and_masker(mask=False)
        trim_rec['vector']  = [(0, 0), (8, 12)]
        trim_rec['quality'] = []
        trim_rec['mask']    = [(0, 3), (5, 6)]
        seq2 = seq_trimmer(seq)
        assert seq2.seq == 'GGTCTCA'
        assert seq2.annotations == {'trimming_recommendations':
                                                    {'mask': [(0, 2), (4, 5)]}}

        seq_trimmer = create_seq_trim_and_masker()
        seq3 = seq_trimmer(seq2)
        assert seq3.seq == 'ggtCtcA'

    @staticmethod
    def test_sequence_masker():
        'It test the sequence masker'
        seq1 = 'ATGTGATGATGATA'
        seq = SeqWithQuality(Seq(seq1), qual=[30] * len(seq1),
                             annotations={'trimming-recomendations':{}})
        segments = [(0, 5) , (9, len(seq1))]
        seq = _mask_sequence(seq, segments)
        assert str(seq.seq) == 'atgtgaTGAtgata'

class SegmentsTests(unittest.TestCase):
    'Here we test seq segments functions'
    @staticmethod
    def test_get_all_segments():
        'Give a list of discontinious segments we get all segments'
        segments =  _get_all_segments([(0, 10), (15, 20)], 30)
        assert segments == [((0, 10), True), ((11, 14), False),
                            ((15, 20), True), ((21, 29), False)]
        segments =  _get_all_segments([(15, 20)], 30)
        assert segments == [((0, 14), False),
                            ((15, 20), True), ((21, 29), False)]
        segments =  _get_all_segments([(15, 29)], 30)
        assert segments == [((0, 14), False), ((15, 29), True)]

    @staticmethod
    def test_non_matched():
        'Given a list of segments we get the complementary matches'
        segments = _get_non_matched_from_matched_locations([(0, 10), (15, 20)], 30)
        assert segments == [(11, 14), (21, 29)]


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'SeqCleanerTest.test_strip_vector_exonerate']
    unittest.main()
