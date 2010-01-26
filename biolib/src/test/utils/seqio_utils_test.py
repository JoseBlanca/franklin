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
import StringIO, tempfile

from biolib.seq.seqs import SeqWithQuality
from biolib.utils.seqio_utils import (seqs_in_file, guess_seq_file_format,
                                      temp_fasta_file, FileSequenceIndex,
                                      quess_seq_type, cat, seqio,
    write_seqs_in_file)
from Bio.Seq import  Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

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
        assert quess_seq_type(fhand, 'fasta', 6) == 'short_seqs'
        fhand = StringIO.StringIO(seqs)
        assert quess_seq_type(fhand, 'fasta', 3) == 'long_seqs'



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

class TestSeqio(unittest.TestCase):
    'It test the converter'
    @staticmethod
    def test_fastq_to_fasta_qual():
        'It tests the conversion from fastq to fasta'
        fcontent  = '@seq1\n'
        fcontent += 'CCCT\n'
        fcontent += '+\n'
        fcontent += ';;3;\n'
        fcontent += '@SRR001666.1\n'
        fcontent += 'GTTGC\n'
        fcontent += '+\n'
        fcontent += ';;;;;\n'
        fhand = StringIO.StringIO(fcontent)

        out_seq_fhand = StringIO.StringIO()
        out_qual_fhand = StringIO.StringIO()
        seqio(in_seq_fhand=fhand, in_format='fastq',
              out_seq_fhand=out_seq_fhand, out_qual_fhand=out_qual_fhand,
              out_format='fasta')
        assert out_seq_fhand.getvalue() == '>seq1\nCCCT\n>SRR001666.1\nGTTGC\n'
        qual = '>seq1\n26 26 18 26\n>SRR001666.1\n26 26 26 26 26\n'
        assert out_qual_fhand.getvalue() == qual

    @staticmethod
    def test_fastq_to_fastq_solexa():
        'It tests the conversion using the Biopython convert function'
        fcontent  = '@seq1\n'
        fcontent += 'CCCT\n'
        fcontent += '+\n'
        fcontent += ';;3;\n'
        fcontent += '@SRR001666.1\n'
        fcontent += 'GTTGC\n'
        fcontent += '+\n'
        fcontent += ';;;;;\n'
        fhand = StringIO.StringIO(fcontent)

        out_seq_fhand = StringIO.StringIO()
        seqio(in_seq_fhand=fhand, in_format='fastq',
              out_seq_fhand=out_seq_fhand, out_format='fastq-solexa')
        result = '@seq1\nCCCT\n+\nZZRZ\n@SRR001666.1\nGTTGC\n+\nZZZZZ\n'
        assert out_seq_fhand.getvalue() == result


class TestCat(unittest.TestCase):
    'It tests the sequence converter'
    @staticmethod
    def test_cat():
        'It tests the cat function'
        inh1 = StringIO.StringIO('>seq1\nACTG\n')
        inh2 = StringIO.StringIO('>seq2\nGTCA\n')
        outh = StringIO.StringIO()
        cat(infiles=[inh1, inh2], outfile=outh)
        assert outh.getvalue() == '>seq1\nACTG\n>seq2\nGTCA\n'

        #it works also with None Values
        outh = StringIO.StringIO()
        cat(infiles=[None, None], outfile=outh)
        assert outh.getvalue() == ''

class TestReprIn_Out(unittest.TestCase):
    'Tests repr write and read functions'
    @staticmethod
    def test_repr():
        'It test seqs in file read and write with repr'
        fhand = StringIO.StringIO()
        seq1 = SeqWithQuality(id='seqid', name='seqname',
                         description='seqdescription', seq=Seq('ATGAT'))
        seq1.letter_annotations["phred_quality"] = [40, 40, 38, 30, 30]
        seq1.annotations['source'] = 'ara1'

        seqfeature = SeqFeature(location=FeatureLocation(5, 8),
                                type='orthologs',
                                qualifiers={'arabidposys':['arab1', 'arab2']})
        seq1.features.append(seqfeature)
        write_seqs_in_file([seq1], fhand, format='repr')
        seqs = seqs_in_file(fhand)
        seq0 = seqs.next()
        assert repr(seq0) == repr(seq1)







if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
