'''
Created on 11/05/2010

@author: peio
'''
from franklin.utils.misc_utils import DATA_DIR
import unittest, StringIO, tempfile
from franklin.seq.readers import (seqs_in_file, guess_seq_file_format,
                                  guess_seq_type, num_seqs_in_file)
from os.path import join

class GuessFormatSeqFileTest(unittest.TestCase):
    'It tests that we can guess the format of a sequence file'
    @staticmethod
    def test_guess_format():
        'It test that we can guess the format for the sequence files'
        fhand = StringIO.StringIO('>fasta\nACTAG\n')
        assert guess_seq_file_format(fhand) == 'fasta'

        fhand = StringIO.StringIO('LOCUS AX0809\n')
        assert guess_seq_file_format(fhand) == 'genbank'

        fhand = StringIO.StringIO('@fastq\nACTAG\n')
        fhand.name = 'hola.sfastq'
        assert guess_seq_file_format(fhand) == 'fastq'


    @staticmethod
    def test_staticmethod():
        'If an empty file is given it should not fail'
        fhand = StringIO.StringIO()
        assert guess_seq_file_format(fhand) is None

class GuessSeqFileTypeTest(unittest.TestCase):
    'It tests that we can guess the type of a sequence file(short|long)'
    @staticmethod
    def test_guess_format():
        'It test that we can guess the format for the sequence files'
        seqs = '>fasta\nACTAG\n>fasta\nACTAG\n>fasta\nACTAG\n>fasta\nACTAG\n'
        fhand = StringIO.StringIO(seqs)
        assert guess_seq_type(fhand, 'fasta', 6) == 'short_seqs'
        fhand = StringIO.StringIO(seqs)
        assert guess_seq_type(fhand, 'fasta', 3) == 'long_seqs'

class SeqsInFileTests(unittest.TestCase):
    'It test that we can get seqrecords out of a seq file.'
    @staticmethod
    def test_seqs_in_file():
        'It test that we get seqs without quality from a sequence file'
        fcontent = '>hola\nACGATCTAGTCATCA\n>caracola\nATCGTAGCTGATGT'
        fhand = StringIO.StringIO(fcontent)
        expected = [('hola', 'ACGATCTAGTCATCA'), ('caracola', 'ATCGTAGCTGATGT')]
        for index, seq in enumerate(seqs_in_file(fhand)):
            assert seq.name == expected[index][0]
            assert str(seq.seq) == expected[index][1]

    def test_seqquals_in_file(self):
        'It test that we get seqs with quality from two sequence files'
        fcontent = '>hola\nACGA\n>caracola\nATCG'
        fhand = StringIO.StringIO(fcontent)
        fcontent_qual = '>hola\n1 2 3 4\n>caracola\n5 6 7 8'
        fhand_qual = StringIO.StringIO(fcontent_qual)
        expected = [('hola', 'ACGA', [1, 2, 3, 4]),
                    ('caracola', 'ATCG', [5, 6, 7, 8])]
        for index, seq in enumerate(seqs_in_file(fhand, fhand_qual)):
            assert seq.name == expected[index][0]
            assert str(seq.seq) == expected[index][1]
            assert seq.qual == expected[index][2]

        #when the seq and qual names do not match we get an error
        fcontent = '>hola\nACGA\n>caracola\nATCG'
        fhand = StringIO.StringIO(fcontent)
        fcontent_qual = '>caracola\n1 2 3 4\n>hola\n5 6 7 8'
        fhand_qual = StringIO.StringIO(fcontent_qual)
        try:
            for seq in seqs_in_file(fhand, fhand_qual):
                #pylint: disable-msg=W0104
                seq.name
            self.fail()
            #pylint: disable-msg=W0704
        except ValueError:
            pass

    @staticmethod
    def test_emptyfile():
        'It should work with empty files'
        fhand = tempfile.NamedTemporaryFile()
        assert not list(seqs_in_file(fhand))

    @staticmethod
    def test_fasq():
        'It test that we can get a seq iter from a fasq file'

        #fasq
        fcontent  = '@seq1\n'
        fcontent += 'CCCT\n'
        fcontent += '+\n'
        fcontent += ';;3;\n'
        fcontent += '@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36\n'
        fcontent += 'GTTGC\n'
        fcontent += '+\n'
        fcontent += ';;;;;\n'
        fhand = StringIO.StringIO(fcontent)

        expected = [('seq1', 'CCCT', [26, 26, 18, 26]),
                    ('SRR001666.1', 'GTTGC', [26, 26, 26, 26, 26])]
        for index, seq in enumerate(seqs_in_file(fhand, format='fastq')):
            assert seq.name == expected[index][0]
            assert str(seq.seq) == expected[index][1]
            assert seq.qual == expected[index][2]

        #fastq-illumina
        fcontent  = '@seq1\n'
        fcontent += 'CCCT\n'
        fcontent += '+\n'
        fcontent += 'AAAA\n'
        fcontent += '@seq2\n'
        fcontent += 'GTTGC\n'
        fcontent += '+\n'
        fcontent += 'AAAAA\n'
        fhand = StringIO.StringIO(fcontent)

        expected = [('seq1', 'CCCT', [1, 1, 1, 1]),
                    ('seq2', 'GTTGC', [1, 1, 1, 1, 1])]
        for index, seq in enumerate(seqs_in_file(fhand,
                                                 format='fastq-illumina')):
            assert seq.name == expected[index][0]
            assert str(seq.seq) == expected[index][1]
            assert seq.qual == expected[index][2]

        #fastq-solexa
        fcontent  = '@seq1\n'
        fcontent += 'CCCT\n'
        fcontent += '+\n'
        fcontent += 'BBBB\n'
        fcontent += '@seq2\n'
        fcontent += 'GTTGC\n'
        fcontent += '+\n'
        fcontent += 'BBBBB\n'
        fhand = StringIO.StringIO(fcontent)

        expected = [('seq1', 'CCCT', [4, 4, 4, 4]),
                    ('seq2', 'GTTGC', [4, 4, 4, 4, 4])]
        for index, seq in enumerate(seqs_in_file(fhand,
                                                 format='fastq-solexa')):
            assert seq.name == expected[index][0]
            assert str(seq.seq) == expected[index][1]
            assert seq.qual == expected[index][2]

class TestNumSeqsInFile(unittest.TestCase):
    'tests num_seqs_in_file'

    @staticmethod
    def test_num_seqs_in_file():
        'tests num_seqs_in_file'
        fasta_fhand = open(join(DATA_DIR, 'seq.fasta'))
        assert num_seqs_in_file(fasta_fhand) == 6

        fastq_fhand = open(join(DATA_DIR, 'solexa.fastq'))
        assert num_seqs_in_file(fastq_fhand, format='sfastq') == 3

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
