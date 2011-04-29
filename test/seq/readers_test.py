'''
Created on 11/05/2010

@author: peio
'''
from franklin.utils.misc_utils import TEST_DATA_DIR
import unittest, StringIO, tempfile, os

from Bio.Alphabet import SingleLetterAlphabet, DNAAlphabet
from Bio.SeqFeature import ExactPosition, FeatureLocation

from franklin.seq.writers import write_seqs_in_file
from franklin.seq.readers import (seqs_in_file, guess_seq_file_format,
                                  guess_seq_type, num_seqs_in_file,
                                  _cast_to_class, fasta_contents_in_file)
from franklin.seq.seqs import Seq, SeqWithQuality, SeqFeature
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

        fhand = StringIO.StringIO('@fastq\nACTAG\n+\nAt+AA')
        fhand.name = 'hola.fastq'
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

    @staticmethod
    def test_repr():
        'It test the repr reader'

        assert u'a (h)' == _cast_to_class("u'a (h)'")

        assert ('adios',) == _cast_to_class("('adios',)")
        assert ['arab1', 'arab2'] == _cast_to_class("['arab1', 'arab2']")
        assert ('arab1', 'arab2') == _cast_to_class("('arab1', 'arab2')")
        result = _cast_to_class("{1: 2}")
        assert {1: 2}  == result

        assert {'al': {'c': 1}, 'T': 2} ==  _cast_to_class("{'al': {'c': 1}, 'T': 2}")

        seq1 = SeqWithQuality(seq=Seq('ATCT'))
        seq2 = SeqWithQuality(seq=Seq('AAAA'))
        fcontent = repr(seq1) + '\n' + repr(seq2) + '\n'
        fhand = StringIO.StringIO(fcontent)
        seqs = list(seqs_in_file(fhand, format='repr'))
        assert repr(seqs[0]) == repr(seq1)
        assert repr(seqs[1]) == repr(seq2)

        #with quality
        seq1 = SeqWithQuality(seq=Seq('ATCT'), qual=[10, 2, 3, 4])
        fcontent = repr(seq1) + '\n'
        fhand = StringIO.StringIO(fcontent)
        seqs = list(seqs_in_file(fhand, format='repr'))
        assert repr(seqs[0]) == repr(seq1)

        #a seq with features
        seq1 = SeqWithQuality(seq=Seq('GAAAAGATGTG', SingleLetterAlphabet()),
                             id='seq', name='seq', description='', dbxrefs=[],
                    features=[SeqFeature(FeatureLocation(ExactPosition(478),
                                                         ExactPosition(478)),
                                         type='intron',
                                         qualifiers={'db': 'tomato'} ),
                              SeqFeature(FeatureLocation(ExactPosition(572),
                                                         ExactPosition(572)),
                                         type='intron',
                                         qualifiers={'db': 'tomato'} )],
                             annotations={}, qual=None)
        fcontent = repr(seq1) + '\n'
        fhand = StringIO.StringIO(fcontent)
        seq0 = list(seqs_in_file(fhand, format='repr'))[0]

        assert seq0.seq == seq1.seq
        assert seq0.qual == seq1.qual
        assert seq0.description == seq1.description
        assert seq0.annotations == seq1.annotations
        feat0 = seq0.features[0]
        feat1 = seq1.features[0]
        assert feat0.type == feat1.type
        assert feat0.qualifiers == feat1.qualifiers
        assert str(feat0.location) == str(feat1.location)

        #some more qualifiers
        fhand = tempfile.NamedTemporaryFile(suffix='.repr')
        seq1 = SeqWithQuality(id='seqid', name='seqname',
                         description='seqdescription', seq=Seq('ATGAT'))
        seq1.letter_annotations["phred_quality"] = [40, 40, 38, 30, 30]
        seq1.annotations['source'] = 'ara1'

        seqfeature = SeqFeature(location=FeatureLocation(5, 8),
                                type='orthologs',
                                qualifiers={'arabidposis':['arab1', 'arab2']})
        seq1.features.append(seqfeature)

        fcontent = repr(seq1) + '\n'
        fhand = StringIO.StringIO(fcontent)
        seq0 = list(seqs_in_file(fhand, format='repr'))[0]

        assert seq0.seq == seq1.seq
        assert seq0.qual == seq1.qual
        assert seq0.description == seq1.description
        assert seq0.annotations == seq1.annotations
        feat0 = seq0.features[0]
        feat1 = seq1.features[0]
        assert feat0.type == feat1.type
        assert feat0.qualifiers == feat1.qualifiers
        assert str(feat0.location) == str(feat1.location)

        #with snps
        repr_ = "SeqWithQuality(seq=Seq('GGGGATTTG', Alphabet()), features=[SeqFeature(FeatureLocation(ExactPosition(213),ExactPosition(213)), type='snv', qualifiers={'alleles': {('C', 3): {'read_groups': ['group1+454', 'group1+454', 'group1+454'], 'qualities': [44.0, 44.0, 44.0], 'libraries': ['group1', 'group1', 'group1'], 'read_names': ['seq1', 'seq4', 'seq7'], 'orientations': [True, True, True], 'samples': ['group1+454', 'group1+454', 'group1+454'], 'quality': 66.0, 'mapping_qualities': [149, 149, 149]}, ('T', 0): {'read_groups': ['group1+454', 'group1+454', 'group1+454', 'group1+454', 'group1+454', 'group1+454'], 'qualities': [44.0, 44.0, 44.0, 44.0, 44.0, 44.0], 'libraries': ['group1', 'group1', 'group1', 'group1', 'group1', 'group1'], 'read_names': ['seq2', 'seq3', 'seq5', 'seq6', 'seq8', 'seq9'], 'orientations': [True, True, True, True, True, True], 'samples': ['group1+454', 'group1+454', 'group1+454', 'group1+454', 'group1+454', 'group1+454'], 'quality': 66.0, 'mapping_qualities': [28, 28, 28, 28, 28, 28]}}, 'reference_allele': 'C'} )])\n"
        alleles = {('C', 3):
                     {'read_groups': ['group1+454', 'group1+454', 'group1+454'],
                      'qualities': [44.0, 44.0, 44.0],
                      'libraries': ['group1', 'group1', 'group1'],
                      'read_names': ['seq1', 'seq4', 'seq7'],
                      'orientations': [True, True, True],
                      'samples': ['group1+454', 'group1+454', 'group1+454'],
                      'quality': 66.0, 'mapping_qualities': [149, 149, 149]},
                   ('T', 0):
                     {'read_groups': ['group1+454', 'group1+454', 'group1+454',
                                      'group1+454', 'group1+454', 'group1+454'],
                      'qualities': [44.0, 44.0, 44.0, 44.0, 44.0, 44.0],
                      'libraries': ['group1', 'group1', 'group1', 'group1',
                                    'group1', 'group1'],
                      'read_names': ['seq2', 'seq3', 'seq5', 'seq6', 'seq8',
                                     'seq9'],
                      'orientations': [True, True, True, True, True, True],
                      'samples': ['group1+454', 'group1+454', 'group1+454',
                                  'group1+454', 'group1+454', 'group1+454'],
                      'quality': 66.0,
                      'mapping_qualities': [28, 28, 28, 28, 28, 28]}
                  }
        fcontent = repr_
        fhand = StringIO.StringIO(fcontent)
        seq0 = list(seqs_in_file(fhand, format='repr'))[0]
        alleles0 = seq0.features[0].qualifiers['alleles']
        assert alleles == alleles0

    @staticmethod
    def test_json_reader():
        'It tests the json sequence writer'
        #first we write some files
        seq0 = SeqWithQuality(seq=Seq('ATGATAGATAGATGF'), name='seq1')
        seq1 = SeqWithQuality(seq=Seq('GATACCA', DNAAlphabet()), name='seq2')
        fhand = tempfile.NamedTemporaryFile(suffix='.json')
        write_seqs_in_file([seq0, seq1], fhand, format='json')
        fhand.flush()

        #now we read them
        seqs = list(seqs_in_file(open(fhand.name)))
        assert seqs[0].seq == seq0.seq
        assert seqs[1].seq == seq1.seq
        assert str(seqs[1].seq.alphabet) == str(seq1.seq.alphabet)

    @staticmethod
    def test_csfasta_reader():
        'It test a csfasta reader'
        seq_fhand = open(os.path.join(TEST_DATA_DIR, 'seq.csfasta'))
        qual_fhand = open(os.path.join(TEST_DATA_DIR, 'solid_qual.qual'))

        seqs = list(seqs_in_file(seq_fhand, qual_fhand, format='csfasta'))
        assert '121101332.0133.2221.23.2.21' in str(seqs[0].seq)
        assert len(seqs) == 3

        seqs = list(seqs_in_file(seq_fhand, qual_fhand, format='csfasta',
                                 double_encoding=True))
        assert 'TTGNACTTNGGGCNGTNGNGCA' in str(seqs[0].seq)


    @staticmethod
    def test_fasta_content_iterator():
        'it test fasta_content_iterator'
        fhand = open(os.path.join(TEST_DATA_DIR, 'seq.fasta'))
        fastas = list(fasta_contents_in_file(fhand))
        assert fastas[0][0] == 'seq1'
        assert fastas[1][1] == 'polya'
        assert 'TAGTCTATGATGCATCAGATGCATGA' in fastas[2][2]
class TestNumSeqsInFile(unittest.TestCase):
    'tests num_seqs_in_file'

    @staticmethod
    def test_num_seqs_in_file():
        'tests num_seqs_in_file'
        fasta_fhand = open(join(TEST_DATA_DIR, 'seq.fasta'))
        assert num_seqs_in_file(fasta_fhand) == 7

        fastq_fhand = open(join(TEST_DATA_DIR, 'solexa.fastq'))
        assert num_seqs_in_file(fastq_fhand, format='sfastq') == 3

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
