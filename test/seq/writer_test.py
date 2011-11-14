'''
Created on 25/03/2010

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

import unittest, json
from StringIO import StringIO

from Bio.Alphabet import ProteinAlphabet, DNAAlphabet

from franklin.seq.seqs import Seq, SeqWithQuality
from franklin.seq.writers import (SsrWriter, OrfWriter,
                                  OrthologWriter, temp_fasta_file,
                                  write_seqs_in_file, SamWriter)
from franklin.seq.readers import seqs_in_file

from tempfile import NamedTemporaryFile
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition

class WriterTest(unittest.TestCase):
    'It test all writers for seqs'

    @staticmethod
    def test_ortologth_writer():
        'It test ssr writer'
        annotations = {'arabidpsis-orthologs':['ara1', 'ara2'],
                       'melon-orthologs':['melo1', 'melo2']}

        seq = SeqWithQuality(seq=Seq('CTTCATCCATTCTCTCATCCGCCGNTGTGGCCTTTGN'),
                             id='seq1', name='seq1', description='Some desc',
                             dbxrefs=[], features=[], annotations=annotations)
        fhand = NamedTemporaryFile(mode='a')

        ortologwriter = OrthologWriter(fhand)
        ortologwriter.write(seq)
        seq_result = open(fhand.name).read()
        assert  'ara1,ara2' in seq_result
        assert  'melo1,melo2' in seq_result

    @staticmethod
    def test_orf_writer():
        'It test ssr writer'
        qualifiers = {'pep': Seq('LHPFSHPPXWPLX', ProteinAlphabet()),
                      'dna': Seq('ATGGCTTCATCCATTCTCTCATCCGCCG', DNAAlphabet()),
                      'strand': 'reverse'}
        orf_feature = SeqFeature(FeatureLocation(ExactPosition(61),
                                                 ExactPosition(471)),
                                                 type='orf',
                                                 qualifiers=qualifiers)

        features = [orf_feature]
        seq = SeqWithQuality(seq=Seq('CTTCATCCATTCTCTCATCCGCCGNTGTGGCCTTTGN'),
                             id='seq1', name='seq1', description='Some desc',
                             dbxrefs=[], features=features, annotations={})
        fhand = NamedTemporaryFile(mode='a')
        pep_fhand = NamedTemporaryFile(mode='a')
        orfwriter = OrfWriter(fhand, pep_fhand)
        orfwriter.write(seq)
        seq_result = open(fhand.name).read()
        pep_result = open(pep_fhand.name).read()
        assert 'ATGGCTTCATCCATTCTCTCATCCGCCG' in seq_result
        assert 'LHPFSHPPXWPLX' in pep_result

    @staticmethod
    def test_ssr_writer():
        'It test ssr writer'
        srr_feature = SeqFeature(FeatureLocation(ExactPosition(0),
                                                 ExactPosition(29)),
                                                 type='microsatellite',
                                                 qualifiers={'score': 27,
                                                       'type': 'trinucleotide',
                                                       'unit': 'ATC'})
        features = [srr_feature]
        seq = SeqWithQuality(seq=Seq('CTTCATCCATTCTCTCATCCGCCGNTGTGGCCTTTGN'),
                             id='seq1', name='seq1', description='Some desc',
                             dbxrefs=[], features=features, annotations={})
        fhand = NamedTemporaryFile(mode='a')
        ssrwriter = SsrWriter(fhand)
        ssrwriter.write(seq)
        result = open(fhand.name).read()
        assert "seq1\t1\t30\t30\t27\ttrinucleotide\tATC" in result


    @staticmethod
    def test_sam_writer():
        'It test sam writer'
        reference = '>ref1\natatagatagatagatagat'
        reads = '>read1\ngatgatagatgatagata'
        ref_fhand = NamedTemporaryFile()
        read_fhand = NamedTemporaryFile()
        out_fhand = NamedTemporaryFile()

        ref_fhand.write(reference)
        read_fhand.write(reads)
        ref_fhand.flush()
        read_fhand.flush()
        ref_fhand.seek(0)
        read_fhand.seek(0)


        alignment = {'reference_name':'ref1',
                     'query_name': 'read1',
                     'mapped': False,
                     'strand':'+',
                     'position': 1,
                     'query_position':1,
                     'cigar': '12M',
                     'mapq':255,
                     }
        sam_writer = SamWriter(ref_fhand, read_fhand, out_fhand)
        sam_writer.write(alignment)

        result = open(out_fhand.name).read()
        line = 'read1\t4\tref1\t1\t255\t12M6S\t*\t0\t0\tgatgatagatgatagata\t*'
        assert line in result


        sam_writer = SamWriter(ref_fhand, read_fhand, out_fhand,
                               keep_unmapped=False)
        sam_writer.write(alignment)

        result = open(out_fhand.name).read()
        line = 'read1\t4\tref1\t1\t255\t12M6S\t*\t0\t0\tgatgatagatgatagata\t*'
        assert line in result

        reference='>ref1\natatagatagatagatagat'
        reads='@read1\ngatgatagatgatagata\n+\ngatgatagatgatagata'
        ref_fhand = NamedTemporaryFile()
        read_fhand = NamedTemporaryFile(suffix='.sfastq')

        out_fhand = NamedTemporaryFile()

        ref_fhand.write(reference)
        read_fhand.write(reads)
        ref_fhand.flush()
        read_fhand.flush()
        ref_fhand.seek(0)
        read_fhand.seek(0)


        alignment = {'reference_name':'ref1',
                     'query_name': 'read1',
                     'mapped': False,
                     'position': 1,
                     'cigar': '12M',
                     'mapq':255,
                     'opt_files':[]}
        sam_writer = SamWriter(ref_fhand, read_fhand, out_fhand)
        sam_writer.write(alignment)

        result = open(out_fhand.name).read()
        line = 'read1\t4\tref1\t1\t255\t12M6S\t*\t0\t0\tgatgatagatgatagata'

        assert line in result

class SequenceWriter(unittest.TestCase):
    'It tests that we can write a stream of SeqWithQuality into a file'

    @staticmethod
    def test_json_writer():
        'It tests the json sequence writer'
        seq0 = SeqWithQuality(seq=Seq('ATGATAGATAGATGF'), name='seq1')
        alleles = {('G', 3): {}}
        filters = {'a_filter':{('param',):False}}
        snv_feature = SeqFeature(FeatureLocation(ExactPosition(3),
                                                 ExactPosition(3)),
                                                 type='snv',
                                        qualifiers={'alleles':alleles,
                                                    'filters':filters})
        seq1 = SeqWithQuality(seq=Seq('GATACCA'), name='seq2',
                              features=[snv_feature])
        fhand = StringIO()
        write_seqs_in_file([seq0, seq1], fhand, format='json')
        lines = fhand.getvalue().splitlines()
        struct1 = json.loads(lines[2])
        assert struct1['seq']['seq'] == 'GATACCA'
        assert struct1['features'][0]['qualifiers']['alleles'].keys()[0] == "('G', 3)"

        fhand.seek(0)
        seqs = list(seqs_in_file(fhand))
        assert seqs[1].features[0].qualifiers['alleles'] == alleles
        assert seqs[1].features[0].qualifiers['filters'] == filters

    @staticmethod
    def test_pickle_writer():
        'It tests the pickle sequence writer'
        seq0 = SeqWithQuality(seq=Seq('ATGATAGATAGATGF'), name='seq1')
        alleles = {('G', 3): {}}
        filters = {'a_filter':{('param',):False}}
        snv_feature = SeqFeature(FeatureLocation(ExactPosition(3),
                                                 ExactPosition(3)),
                                                 type='snv',
                                        qualifiers={'alleles':alleles,
                                                    'filters':filters})
        seq1 = SeqWithQuality(seq=Seq('GATACCA'), name='seq2',
                              features=[snv_feature])
        fhand = StringIO()
        write_seqs_in_file([seq0, seq1], fhand, format='pickle')
        #print fhand.getvalue()

        fhand.seek(0)
        seqs = list(seqs_in_file(fhand))
        assert seqs[1].features[0].qualifiers['alleles'] == alleles
        assert seqs[1].features[0].qualifiers['filters'] == filters


class TestFastaFileUtils(unittest.TestCase):
    'Here we test a couple of utilities related to fast format'

    @staticmethod
    def test_temp_fasta_file_one_seq():
        'It test temp_fasta_file'
        seqrec1 = SeqWithQuality(seq=Seq('ATGATAGATAGATGF'), name='seq1')
        fhand = temp_fasta_file(seqrec1)
        content = open(fhand.name).read()
        assert content == ">seq1\nATGATAGATAGATGF\n"

    @staticmethod
    def test_temp_fasta_file_seq_iter():
        'It test temp_fasta_file'
        seqrec1 = SeqWithQuality(seq=Seq('ATGATAGATAGATGF'), name='seq1')
        seqrec2 = SeqWithQuality(seq=Seq('ATGATAGATAGA'), name='seq2')
        seq_iter = iter([seqrec1, seqrec2])
        fhand = temp_fasta_file(seq_iter)
        content = open(fhand.name).read()
        assert content == ">seq1\nATGATAGATAGATGF\n>seq2\nATGATAGATAGA\n"

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_gff_writer']
    unittest.main()
