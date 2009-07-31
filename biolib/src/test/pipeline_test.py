'''
Created on 2009 uzt 30

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

import unittest, os
import biolib
from tempfile import NamedTemporaryFile

from biolib.pipelines import  (configure_pipeline, seq_pipeline_runner,
                               pipeline_runner)
from biolib.biolib_seqio_utils import seqs_in_file
from biolib.contig_io import get_parser



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



DATA_DIR = os.path.join(os.path.split(biolib.__path__[0])[0], 'data')

class PipelineTests(unittest.TestCase):
    'It test pipeline related functions'

    def test_configure_pipeline(self):
        'It tests configure pipeline'
        pipeline      = 'sanger_with_qual'
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
        pipeline = 'sanger_with_qual'

        fhand_adaptors = NamedTemporaryFile()
        fhand_adaptors.write(ADAPTORS)
        fhand_adaptors.flush()

        configuration = {'remove_vectors': {'vectors':'UniVec'},
                         'remove_adaptors':{'vectors':fhand_adaptors.name}}

        seq_fhand  = open(os.path.join(DATA_DIR, 'seq.fasta'), 'r')
        qual_fhand = open(os.path.join(DATA_DIR, 'qual.fasta'), 'r')

        seq_iter = seqs_in_file(seq_fhand, qual_fhand)

        filtered_seq_iter = pipeline_runner(pipeline, seq_iter, configuration)

        seq_list = list(filtered_seq_iter)
        assert 'ATCGCGAtcgggggg' in str(seq_list[0].seq)
        assert len(seq_list) == 6

    @staticmethod
    def test_seq_pipeline_run():
        'It tests that the pipeline runs ok'
        pipeline = 'sanger_with_qual'

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

        seq_pipeline_runner(pipeline, configuration, io_fhands)
        io_fhands['out_seq'].seek(0)
        result_seq = io_fhands['out_seq'].read()
        assert result_seq.count('>') == 6

    @staticmethod
    def test_contig_pipeline_run():
        'It test the contig clean pipeline'
        pipeline = 'contig_clean'
        fhand = open(os.path.join(DATA_DIR, 'example.caf'), 'r')
        caf_parser = get_parser(fhand, format='caf')
        contigs = caf_parser.contigs()
        contigs = pipeline_runner(pipeline, items=contigs)
        assert  contigs.next().consensus.name == 'Contig1'







if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
