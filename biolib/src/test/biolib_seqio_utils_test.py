'''
Created on 2009 uzt 28

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
import StringIO

from biolib.seqs import SeqWithQuality
from biolib.biolib_seqio_utils import (seqs_in_file, guess_seq_file_format,
                                       temp_fasta_file, FileSequenceIndex)

class GuessFormatSeqFileTest(unittest.TestCase):
    'It tests that we can guess the format of a sequence file'
    @staticmethod
    def test_guess_format():
        'It test that we can guess the format for the sequence files'
        fhand = StringIO.StringIO('>fasta\nACTAG\n')
        assert guess_seq_file_format(fhand) == 'fasta'

        fhand = StringIO.StringIO('LOCUS AX0809\n')
        assert guess_seq_file_format(fhand) == 'genbank'

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
        except RuntimeError:
            pass

class TestFastaFileUtils(unittest.TestCase):
    'here we test our utils related to fast format'

    @staticmethod
    def test_temp_fasta_file_one_seq():
        'It test temp_fasta_file'
        seqrec1 = SeqWithQuality(seq='ATGATAGATAGATGF', name='seq1')
        fhand = temp_fasta_file(seqrec1)
        content = open(fhand.name).read()
        assert content == ">seq1\nATGATAGATAGATGF\n"

    @staticmethod
    def test_temp_fasta_file_seq_iter():
        'It test temp_fasta_file'
        seqrec1 = SeqWithQuality(seq='ATGATAGATAGATGF', name='seq1')
        seqrec2 = SeqWithQuality(seq='ATGATAGATAGA', name='seq2')
        seq_iter = iter([seqrec1, seqrec2])
        fhand = temp_fasta_file(seq_iter)
        content = open(fhand.name).read()
        assert content == ">seq1\nATGATAGATAGATGF\n>seq2\nATGATAGATAGA\n"

class TestSequenceFileIndexer(unittest.TestCase):
    'It test the FileSequenceIndex class'

    @staticmethod
    def test_sequence_fasta_index():
        'It test the file index class basic functionality'
        fhand = StringIO.StringIO('>seq1 una seq\nACTG\n>seq2 otra seq\nGTAC\n')
        index = FileSequenceIndex(fhand)
        seqrec1 = index['seq1']
        seqrec2 = index['seq2']
        assert seqrec1.name == 'seq1'
        assert seqrec1.description == 'una seq'
        assert seqrec2.seq == 'GTAC'

    @staticmethod
    def test_sequence_fasta_qual_index():
        'It test the file index class basic functionality'
        fhand = StringIO.StringIO('>seq1 una seq\n1 2 3\n>seq2 otra seq\n10\n')
        index = FileSequenceIndex(fhand)
        seqrec1 = index['seq1']
        seqrec2 = index['seq2']
        assert seqrec1.name == 'seq1'
        assert seqrec1.description == 'una seq'
        assert seqrec1.qual == [1, 2, 3]
        assert seqrec2.qual == [10]
